
import click as ck
import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
import math
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_curve, auc
from utils import EarlyStopping
from step0_TrainTestSetting_local import *
from step0_TrainTestSetting_global import path_base
from tqdm import tqdm
import os

print('\n\n ----- starting step2 --------------')


@ck.command()
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='train data file')
@ck.option('--test-data-file', '-tedf', default='data/test_data.pkl', help='XX')
@ck.option('--aa-ss', '-aass', default=params_local['aa_ss'], help='aa/ss8/ss3')
@ck.option('--path-base', '-pb', default=params_local['path_base'], help='..')
@ck.option('--model-file', '-mf', default='data/model_checkpoint.pth', help='XX')
@ck.option('--predictions-file', '-pf', default='data/predictions.pkl', help='XX')
@ck.option('--logger-file', '-lf', default='data/training.csv', help='XX')  #
@ck.option('--terms-file', '-tf', default='data/terms_gominre_trxte.pkl', help='XX')
@ck.option('--kernels-list', '-kl', default=params_local['kernels'], type=list, help='XX')
@ck.option('--filters-list', '-fl', default=params_local['filters'], type=list, help='XX')
@ck.option('--prot-letter-aa', '-plaa', default=params_local['PROT_LETTER_aa'], type=list, help='XX')
@ck.option('--prot-letter-ss8', '-plss8', default=params_local['PROT_LETTER_ss8'], type=list, help='XX')
@ck.option('--prot-letter-ss3', '-plss3', default=params_local['PROT_LETTER_ss3'], type=list, help='XX')
@ck.option('--run-train-switch', '-rtrs', default='F', help='XX')
@ck.option('--run-test-switch', '-rtes', default='T', help='XX')


def main(go_file, train_data_file, test_data_file, aa_ss, path_base, model_file, predictions_file,
         kernels_list, filters_list, logger_file, terms_file, prot_letter_aa, prot_letter_ss8, prot_letter_ss3,
         run_train_switch, run_test_switch):
    print('\n################## a long, long time ago ... ##################\n')

    ##################### pre_setting
    PROT_LETTER = []
    PROT_LETTER_len = -1
    if aa_ss == 'aa':
        PROT_LETTER = prot_letter_aa
        PROT_LETTER_len = len(PROT_LETTER)  # =20
    elif aa_ss == 'ss8':
        PROT_LETTER = prot_letter_ss8
        PROT_LETTER_len = len(PROT_LETTER)  # =8
    elif aa_ss == 'ss3':
        PROT_LETTER = prot_letter_ss3
        PROT_LETTER_len = len(PROT_LETTER)  # =3

    PROT_INDEX = dict()
    for i in range(len(PROT_LETTER)):
        PROT_INDEX[PROT_LETTER[i]] = i + 1

    kernels_tuple = [(value, PROT_LETTER_len + 1) for value in kernels_list]
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}
    terms_classes = len(terms_dict)



    ########################### train + valid ####################################
    train_df, valid_df = split_train_valid(train_data_file, terms, params_local['TrainValidRatio'])
    test_df = pd.read_pickle(test_data_file)
    train_dataset = ProteinGODataset(train_df, terms_dict, PROT_LETTER_len, PROT_INDEX)
    valid_dataset = ProteinGODataset(valid_df, terms_dict, PROT_LETTER_len, PROT_INDEX)
    test_dataset = ProteinGODataset(test_df, terms_dict, PROT_LETTER_len, PROT_INDEX)

    train_dataloader = DataLoader(train_dataset, batch_size=params_local['batch_size'], shuffle=True)
    valid_dataloader = DataLoader(valid_dataset, batch_size=params_local['batch_size'], shuffle=False)
    test_dataloader = DataLoader(test_dataset, batch_size=params_local['batch_size'], shuffle=False)

    device_ids = params_local['device_ids']  # cuda:3"   or    [0, 1, 2, 3]

    if isinstance(device_ids, str):  #
        print('--- single GPU = ', device_ids)
        device = torch.device(device_ids if torch.cuda.is_available() else "cpu")  #
        model = Model1(terms_classes, params_local, kernels_tuple, filters_list)
        model.to(device)

    else:
        print('--- multi GPU = ', device_ids)
        device = torch.device("cuda:{}".format(device_ids[0]) if torch.cuda.is_available() else "cpu")

        model = Model1(terms_classes, params_local, kernels_tuple, filters_list)
        model.to(device)
        model = nn.DataParallel(model, device_ids=device_ids)
        model.to(device)

    print('---------------- This is model: ----------------')
    print('model shape by general model ::::')
    print(model)



    optimizer = optim.Adam(model.parameters(), lr=params_local['learning_rate'])  # lr=3e-4
    criterion = nn.BCELoss()

    if run_train_switch == 'T':
        run_train(model, train_dataloader, valid_dataloader, device, optimizer, criterion, logger_file)



    ###########################  test ############################################

    if params_local['load_pretrained_model'] == 1:
        model_file = params_local['load_pretrained_model_addr']

    model.load_state_dict(torch.load(model_file))
    model.to(device)

    if run_test_switch == 'T':
        run_test(model, test_dataloader, device, criterion, test_df, predictions_file)


    print('\n################## And they all lived happily ever after! ##################\n')




def to_onehot(seq, PROT_LETTER_len, PROT_INDEX, start=0):
    onehot = np.zeros((params_local['MAXLEN'], PROT_LETTER_len + 1), dtype=np.int32)
    l = min(params_local['MAXLEN'], len(seq))
    for i in range(start, start + l):
        onehot[i, PROT_INDEX.get(seq[i - start], 0)] = 1
    onehot[0:start, 0] = 1
    onehot[start + l:, 0] = 1
    return onehot


def split_train_valid(data_file, terms, trainvalidratio):
    df = pd.read_pickle(data_file)
    n = len(df)
    index = np.arange(n)
    train_n = int(n * trainvalidratio)
    np.random.seed(seed=0)
    np.random.shuffle(index)
    train_df = df.iloc[index[:train_n]]
    valid_df = df.iloc[index[train_n:]]
    return train_df, valid_df

def compute_roc(labels, preds):
    fpr, tpr, _ = roc_curve(labels.flatten(), preds.flatten())
    roc_auc = auc(fpr, tpr)
    return roc_auc

class ProteinGODataset(Dataset):
    def __init__(self, df, terms_dict, PROT_LETTER_len, PROT_INDEX):
        self.df = df
        self.terms_dict = terms_dict
        self.PROT_LETTER_len = PROT_LETTER_len
        self.PROT_INDEX = PROT_INDEX

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        terms_classes = len(self.terms_dict)
        row = self.df.iloc[idx]
        seq = row.sequences
        onehot = to_onehot(seq, self.PROT_LETTER_len, self.PROT_INDEX, start=0)
        data_onehot = torch.zeros((params_local['MAXLEN'], self.PROT_LETTER_len + 1), dtype=torch.float32)
        labels = torch.zeros(terms_classes,
                             dtype=torch.float)
        data_onehot[:onehot.shape[0], :] = torch.from_numpy(
            onehot)

        for t_id in row.prop_annotations:
            if t_id in self.terms_dict:
                labels[self.terms_dict[t_id]] = 1

        return data_onehot, labels

######## model  ########
class Model1(nn.Module):
    def __init__(self, terms_classes, params, kernels_tuple, filters_list):
        super(Model1, self).__init__()

        self.kernels_tuple = kernels_tuple
        self.filters = filters_list

        self.convs = nn.ModuleList()
        self.pools = nn.ModuleList()
        self.flats = nn.ModuleList()

        for i in range(len(self.kernels_tuple)):
            conv = nn.Conv2d(
                in_channels=1,
                out_channels=self.filters[i],
                kernel_size=self.kernels_tuple[i],
                padding=(0,),
                stride=(1, 1)
            )

            pool = nn.MaxPool2d(kernel_size=(params_local['MAXLEN'] - self.kernels_tuple[i][0] + 1, 1))

            flat = nn.Flatten()
            self.convs.append(conv)
            self.pools.append(pool)
            self.flats.append(flat)

        if params_local['FC_depth'] > 0:
            self.fc_layers = nn.ModuleList()
            for i in range(params_local['FC_depth']):
                if i == 0:
                    fc = nn.Linear(sum(self.filters), terms_classes)
                else:
                    fc = nn.Linear(terms_classes, terms_classes)
                self.fc_layers.append(fc)
            self.output_layer = nn.Linear(terms_classes, terms_classes)

        else:
            self.output_layer = nn.Linear(sum(self.filters), terms_classes)

    def forward(self, x):
        nets = []
        for i in range(len(self.kernels_tuple)):
            conv = self.convs[i](x.unsqueeze(1))
            pool = self.pools[i](conv)
            flat = self.flats[i](pool)
            nets.append(flat)

        if len(self.kernels_tuple) == 1:
            net = nets[0]
        elif len(self.kernels_tuple) > 1:
            net = torch.cat(nets, dim=1)
        else:
            print('error ................')


        if params_local['FC_depth'] > 0:
            for k in range(len(self.fc_layers)):
                net = self.fc_layers[k](net)  #

        net = self.output_layer(net)
        net = torch.sigmoid(net)

        return net


####################### train + valid  ####################


def run_train(model, train_dataloader, valid_dataloader, device, optimizer, criterion, logger_file):

    early_stopping = EarlyStopping(patience=params_local['EarlyStopping_patience'], verbose=True)
    logger = pd.DataFrame(columns=['epoch', 'train_loss', 'valid_loss'])

    for epoch in range(params_local['epochs']):
        model.train()
        train_loss = 0

        train_dataloader_with_progress = tqdm(train_dataloader, desc=f"Epoch {epoch + 1}/{params_local['epochs']}", leave=False)

        for inputs_onehot, target_labels in train_dataloader_with_progress:
            inputs_onehot = inputs_onehot.to(device)
            target_labels = target_labels.to(device)


            optimizer.zero_grad()
            outputs = model(inputs_onehot)

            loss = criterion(outputs, target_labels)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()
            train_dataloader_with_progress.set_postfix({'Loss': torch.tensor(train_loss).item()})

        train_loss /= len(train_dataloader)
        train_dataloader_with_progress.close()

        ################## valid  ###############
        model.eval()
        valid_loss = 0.0

        valid_dataloader_with_progress = tqdm(valid_dataloader, desc="Validation", leave=False)

        with torch.no_grad():  #
            for inputs_onehot, target_labels in valid_dataloader_with_progress:
                inputs_onehot = inputs_onehot.to(device)
                target_labels = target_labels.to(device)

                outputs = model(inputs_onehot)
                valid_loss += criterion(outputs, target_labels).item()

            valid_loss /= len(valid_dataloader)
            valid_dataloader_with_progress.set_postfix({'Loss': torch.tensor(valid_loss).item()})
            valid_dataloader_with_progress.close()

            logger = logger.append({
                'epoch': epoch + 1,
                'train_loss': train_loss,
                'valid_loss': valid_loss
            }, ignore_index=True)

            print(f"Epoch {epoch + 1}/{params_local['epochs']} --------------- train_loss: {train_loss:.4f} - Valid_loss: {valid_loss:.4f}")


            if early_stopping(valid_loss, model):
                print("Early Stopping!")
                break

        logger.to_csv(logger_file, index=False)

        if early_stopping.early_stop:
            print("Early Stopping!")
            break



####################### test #####################

def run_test(model, test_dataloader, device, criterion, test_df, predictions_file):
    model.eval()
    test_loss = 0
    preds_list = []
    labels_list = []

    with torch.no_grad():
        for inputs_onehot, target_labels in test_dataloader:
            inputs_onehot = inputs_onehot.to(device)
            target_labels = target_labels.to(device)

            outputs = model(inputs_onehot)
            test_loss += criterion(outputs, target_labels).item()

            preds_list.append(outputs.cpu().numpy())
            labels_list.append(target_labels.cpu().numpy())

    test_loss /= len(test_dataloader)
    print('\n-------------test_loss = ', test_loss)

    preds = np.concatenate(preds_list)
    labels = np.concatenate(labels_list)

    roc_auc = compute_roc(labels, preds)
    print('\n------------  ROC AUC = ', roc_auc)

    test_df['labels'] = list(labels)
    test_df['preds'] = list(preds)
    test_df.to_pickle(predictions_file)


if __name__ == '__main__':
    main()








