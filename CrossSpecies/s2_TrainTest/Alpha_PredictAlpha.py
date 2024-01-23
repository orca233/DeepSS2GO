#!/usr/bin/env python

import click as ck
import numpy as np
import pandas as pd
import torch.nn as nn
import time
from utils import Ontology, NAMESPACES
import torch
from torch.utils.data import Dataset, DataLoader
import json
from step0_TrainTestSetting_local import *


MAXLEN = params_local['MAXLEN']

@ck.command()
@ck.option('--in-file', '-if', default='data/test_data.fa', help='Input FASTA file', required=True)
@ck.option('--out-file-bp', '-of', default='data/results_bp.csv', help='Output result file')
@ck.option('--out-file-cc', '-of', default='data/results_cc.csv', help='Output result file')
@ck.option('--out-file-mf', '-of', default='data/results_mf.csv', help='Output result file')
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')
@ck.option('--model-file', '-mf', default='data/model_checkpoint.pth', help='Tensorflow model file')
@ck.option('--terms-file', '-tf', default='data/terms_gominre_trxte.pkl', help='List of predicted terms')
@ck.option('--annotations-file', '-tf', default='data/train_data.pkl', help='Experimental annotations')
@ck.option('--chunk-size', '-cs', default=100000, help='Number of sequences to read at a time')
@ck.option('--diamond-file', '-df', default='data/diamond_aa.res', help='Diamond Mapping file')
@ck.option('--threshold', '-t', default=0.1, help='Prediction threshold')
@ck.option('--batch-size', '-bs', default=32, help='Batch size for prediction model')
@ck.option('--test-data-file', '-tedf', default='data/test_data.pkl', help='XX')
@ck.option('--aa-ss', '-aass', default=params_local['aa_ss'], help='aa/ss8/ss3')
@ck.option('--prot-letter-aa', '-plaa', default=params_local['PROT_LETTER_aa'], type=list, help='XX')
@ck.option('--prot-letter-ss8', '-plss8', default=params_local['PROT_LETTER_ss8'], type=list, help='XX')
@ck.option('--prot-letter-ss3', '-plss3', default=params_local['PROT_LETTER_ss3'], type=list, help='XX')
@ck.option('--kernels-list', '-kl', default=params_local['kernels'], type=list, help='XX')
@ck.option('--filters-list', '-fl', default=params_local['filters'], type=list, help='XX')
@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--alpha', '-a', default='json', help='alpha = json(with quote) or 0-1(without quote, eg 0.3 float)')


def main(in_file, out_file_bp, out_file_cc, out_file_mf, go_file, model_file, terms_file,  # out_file
         annotations_file, chunk_size, diamond_file, threshold, batch_size,
         test_data_file, aa_ss, prot_letter_aa, prot_letter_ss8, prot_letter_ss3,
         kernels_list, filters_list, ont, alpha):  # maxlen


    ##################### pre_setting ######################

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

    go = Ontology(go_file, with_rels=True)
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}
    terms_classes = len(terms_dict)



    ############ diamondIC  ##################

    annotations = {}
    df = pd.read_pickle(annotations_file)
    for row in df.itertuples():
        annotations[row.proteins] = set(row.prop_annotations)

    go.calculate_ic(annotations.values())
    diamond_preds = {}
    mapping = {}
    with open(diamond_file, 'rt') as f:
        for line in f:
            it = line.strip().split()
            if it[0] not in mapping:
                mapping[it[0]] = {}
            mapping[it[0]][it[1]] = float(it[2])

    for prot_id, sim_prots in mapping.items():
        annots = {}
        allgos = set()
        total_score = 0.0

        for p_id, score in sim_prots.items():
            allgos |= annotations[p_id]
            total_score += score

        allgos = list(sorted(allgos))
        sim = np.zeros(len(allgos), dtype=np.float32)

        for j, go_id in enumerate(allgos):
            s = 0.0
            for p_id, score in sim_prots.items():
                if go_id in annotations[p_id]:
                    s += score
            sim[j] = s / total_score
        for go_id, score in zip(allgos, sim):
            annots[go_id] = score
        diamond_preds[prot_id] = annots





    ############ check_model.pth  ##################

    test_df = pd.read_pickle(test_data_file)
    test_dataset = ProteinGODataset(test_df, terms_dict, PROT_LETTER_len, PROT_INDEX)
    test_dataloader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)

    device_ids = params_local['device_ids']

    if isinstance(device_ids, str):
        print('--- single GPU = ', device_ids)
        device = torch.device(device_ids if torch.cuda.is_available() else "cpu")
        model = Model1(terms_classes, params_local)
        model.to(device)

    else:
        print('--- multi GPU = ', device_ids)
        device = torch.device("cuda:{}".format(device_ids[0]) if torch.cuda.is_available() else "cpu")
        model = Model1(terms_classes, params_local, kernels_tuple, filters_list)
        model.to(device)
        model = nn.DataParallel(model, device_ids=device_ids)
        model.to(device)


    model.load_state_dict(torch.load(model_file))
    model.to(device)
    model.eval()



    # (1 - alpha - beta) * diamond + alpha * preds_aa + beta * preds_ss8

    alphas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}

    start_time = time.time()
    total_seq = 0

    if ont == 'bp':
        out_file_ont = out_file_bp
    elif ont == 'cc':
        out_file_ont = out_file_cc
    elif ont == 'mf':
        out_file_ont = out_file_mf

    w = open(out_file_ont, 'wt')
    for prot_ids, sequences in read_fasta(in_file, chunk_size):
        total_seq += len(prot_ids)
        deep_preds = {}
        ids, data = get_data(sequences, PROT_LETTER_len, PROT_INDEX)

        preds_list = []
        labels_list = []
        with torch.no_grad():
            for inputs_onehot, target_labels in test_dataloader:
                inputs_onehot = inputs_onehot.to(device)
                target_labels = target_labels.to(device)

                outputs = model(inputs_onehot)
                preds_list.append(outputs.cpu().numpy())
                labels_list.append(target_labels.cpu().numpy())
        preds = np.concatenate(preds_list)
        labels = np.concatenate(labels_list)  #

        assert preds.shape[1] == len(terms)
        for i, j in enumerate(ids):
            prot_id = prot_ids[j]
            if prot_id not in deep_preds:
                deep_preds[prot_id] = {}
            for l in range(len(terms)):
                if preds[i, l] >= 0.01:
                    if terms[l] not in deep_preds[prot_id]:
                        deep_preds[prot_id][terms[l]] = preds[i, l]
                    else:
                        deep_preds[prot_id][terms[l]] = max(deep_preds[prot_id][terms[l]], preds[i, l])


        if alpha == 'json':
            print('alpha is from json, alpha = ', alpha)
            last_release_metadata = 'Alpha_last_release.json'
            with open(last_release_metadata, 'r') as f:
                print('Reading file from json')
                last_release_data = json.load(f)
                alpha = last_release_data["alphas"][ont]
                print('111111111', type(alpha))
                alphas[NAMESPACES[ont]] = alpha
                print('alpha, alphas = ')
                print(alpha, alphas)

        else:
            print('alpha is from click, alpha = ', alpha)
            print('type_alpha = ', type(alpha))
            alpha = float(alpha)
            alphas[NAMESPACES[ont]] = alpha
            print('updated_type_alpha = ', type(alphas[NAMESPACES[ont]]))



        for prot_id in prot_ids:
            annots = {}
            if prot_id in diamond_preds:
                for go_id, score in diamond_preds[prot_id].items():
                    annots[go_id] = score * (1 - alphas[go.get_namespace(go_id)])
            for go_id, score in deep_preds[prot_id].items():
                if go_id in annots:
                    annots[go_id] += alphas[go.get_namespace(go_id)] * score
                else:
                    annots[go_id] = alphas[go.get_namespace(go_id)] * score

            gos = list(annots.keys())
            for go_id in gos:
                for g_id in go.get_ancestors(go_id):
                    if g_id in annots:
                        annots[g_id] = max(annots[g_id], annots[go_id])
                    else:
                        annots[g_id] = annots[go_id]

            sannots = sorted(annots.items(), key=lambda x: x[1], reverse=True)
            for go_id, score in sannots:
                if score >= threshold:

                    if go.get_namespace(go_id) == NAMESPACES[ont]:
                        w.write(prot_id + ', ' + go_id + ', ' + go.get_namespace(go_id) + ', ' + go.get_term(go_id)['name'] + ', %.3f\n' % score)

            w.write('\n')


    w.close()
    total_time = time.time() - start_time
    print('Total prediction time for %d sequences is %d' % (total_seq, total_time))



def to_onehot(seq, PROT_LETTER_len, PROT_INDEX, start=0):
    onehot = np.zeros((params_local['MAXLEN'], PROT_LETTER_len + 1), dtype=np.int32)
    l = min(params_local['MAXLEN'], len(seq))
    for i in range(start, start + l):
        onehot[i, PROT_INDEX.get(seq[i - start], 0)] = 1
    onehot[0:start, 0] = 1
    onehot[start + l:, 0] = 1
    return onehot



class ProteinGODataset(Dataset):
    def __init__(self, df, terms_dict, PROT_LETTER_len, PROT_INDEX):
        self.df = df
        self.terms_dict = terms_dict
        self.PROT_LETTER_len = PROT_LETTER_len
        self.PROT_INDEX = PROT_INDEX

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        terms_classes = len(self.terms_dict)  #
        row = self.df.iloc[idx]
        seq = row.sequences
        onehot = to_onehot(seq, self.PROT_LETTER_len, self.PROT_INDEX, start=0)
        data_onehot = torch.zeros((params_local['MAXLEN'], self.PROT_LETTER_len + 1), dtype=torch.float32)
        labels = torch.zeros(terms_classes,
                             dtype=torch.float)
        data_onehot[:onehot.shape[0], :] = torch.from_numpy(
            onehot)  #

        for t_id in row.prop_annotations:
            if t_id in self.terms_dict:
                labels[self.terms_dict[t_id]] = 1

        return data_onehot, labels



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
            self.convs.append(conv)  # 排排队，吃果果
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
                net = self.fc_layers[k](net)

        net = self.output_layer(net)
        net = torch.sigmoid(net)

        return net




def get_data(sequences, PROT_LETTER_len, PROT_INDEX):
    pred_seqs = []
    ids = []
    for i, seq in enumerate(sequences):
        if len(seq) > MAXLEN:
            st = 0
            while st < len(seq):
                pred_seqs.append(seq[st: st + MAXLEN])
                ids.append(i)
                st += MAXLEN - 128
        else:
            pred_seqs.append(seq)
            ids.append(i)
    n = len(pred_seqs)
    data = np.zeros((n, MAXLEN, PROT_LETTER_len + 1), dtype=np.float32)

    for i in range(n):
        seq = pred_seqs[i]
        data[i, :, :] = to_onehot(seq, PROT_LETTER_len, PROT_INDEX)
    return ids, data

def read_fasta(filename, chunk_size):
    seqs = list()
    info = list()
    seq = ''
    inf = ''
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq != '':
                    seqs.append(seq)
                    info.append(inf)
                    if len(info) == chunk_size:
                        yield (info, seqs)
                        seqs = list()
                        info = list()
                    seq = ''
                inf = line[1:].split()[0]
            else:
                seq += line
        seqs.append(seq)
        info.append(inf)
    yield (info, seqs)


if __name__ == '__main__':
    main()
