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
from step0_TrainTestSetting_local_aa import params_local as params_local_aa
from step0_TrainTestSetting_local_ss8 import params_local as params_local_ss8

MAXLEN = params_local_aa['MAXLEN']

@ck.command()
@ck.option('--in-file-aa', '-if', default='data/test_data_aa.fa', help='Input FASTA file', required=True)
@ck.option('--in-file-ss8', '-if', default='data/test_data_ss8.fa', help='Input FASTA file', required=True)
@ck.option('--out-file-bp', '-of', default='data/results_bp.csv', help='Output result file')
@ck.option('--out-file-cc', '-of', default='data/results_cc.csv', help='Output result file')
@ck.option('--out-file-mf', '-of', default='data/results_mf.csv', help='Output result file')
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')
@ck.option('--model-file-aa', '-mf', default='data/model_checkpoint_aa.pth', help='Tensorflow model file')
@ck.option('--model-file-ss8', '-mf', default='data/model_checkpoint_ss8.pth', help='Tensorflow model file')
@ck.option('--terms-file-aa', '-tf', default='data/terms_gominre_trxte_aa.pkl', help='List of predicted terms')
@ck.option('--terms-file-ss8', '-tf', default='data/terms_gominre_trxte_ss8.pkl', help='List of predicted terms')
@ck.option('--annotations-file', '-tf', default='data/train_data_aa.pkl', help='Experimental annotations')
@ck.option('--chunk-size', '-cs', default=100000, help='Number of sequences to read at a time')
@ck.option('--diamond-file', '-df', default='data/diamond_aa.res', help='Diamond Mapping file')
@ck.option('--threshold', '-t', default=0.1, help='Prediction threshold')
@ck.option('--batch-size', '-bs', default=32, help='Batch size for prediction model')
@ck.option('--test-data-file-aa', '-tedf', default='data/test_data_aa.pkl', help='XX')
@ck.option('--test-data-file-ss8', '-tedf', default='data/test_data_ss8.pkl', help='XX')
@ck.option('--prot-letter-aa', '-plaa', default=params_local_aa['PROT_LETTER_aa'], type=list, help='XX')
@ck.option('--prot-letter-ss8', '-plss8', default=params_local_aa['PROT_LETTER_ss8'], type=list, help='XX')
@ck.option('--prot-letter-ss3', '-plss3', default=params_local_aa['PROT_LETTER_ss3'], type=list, help='XX')
@ck.option('--kernels-list-aa', '-klaa', default=params_local_aa['kernels'], type=list, help='XX')
@ck.option('--kernels-list-ss8', '-klss8', default=params_local_ss8['kernels'], type=list, help='XX')
@ck.option('--filters-list-aa', '-flaa', default=params_local_aa['filters'], type=list, help='XX')
@ck.option('--filters-list-ss8', '-flss8', default=params_local_ss8['filters'], type=list, help='XX')
@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--alpha', '-a', default='json', help='alpha = json(with quote) or 0-1(without quote, eg 0.3 float)')
@ck.option('--beta', '-b', default='json', help='beta = json(with quote) or 0-1(without quote, eg 0.3 float)')

def main(in_file_aa, in_file_ss8, out_file_bp, out_file_cc, out_file_mf, go_file, model_file_aa, model_file_ss8,
         terms_file_aa, terms_file_ss8,  # out_file
         annotations_file, chunk_size, diamond_file, threshold, batch_size,
         test_data_file_aa, test_data_file_ss8, prot_letter_aa, prot_letter_ss8, prot_letter_ss3,
         kernels_list_aa, kernels_list_ss8, filters_list_aa, filters_list_ss8, ont, alpha, beta):


    ##################### pre_setting  ######################
    PROT_LETTER_aa = prot_letter_aa  # list
    PROT_LETTER_len_aa = 20
    PROT_INDEX_aa = dict()
    for i in range(len(PROT_LETTER_aa)):
        PROT_INDEX_aa[PROT_LETTER_aa[i]] = i + 1

    PROT_LETTER_ss8 = prot_letter_ss8  # list
    PROT_LETTER_len_ss8 = 8
    PROT_INDEX_ss8 = dict()
    for i in range(len(PROT_LETTER_ss8)):
        PROT_INDEX_ss8[PROT_LETTER_ss8[i]] = i + 1

    kernels_tuple_aa = [(value, PROT_LETTER_len_aa + 1) for value in kernels_list_aa]
    kernels_tuple_ss8 = [(value, PROT_LETTER_len_ss8 + 1) for value in kernels_list_ss8]
    go = Ontology(go_file, with_rels=True)

    terms_df_aa = pd.read_pickle(terms_file_aa)
    terms_aa = terms_df_aa['terms'].values.flatten()
    terms_dict_aa = {v: i for i, v in enumerate(terms_aa)}
    terms_classes_aa = len(terms_dict_aa)

    terms_df_ss8 = pd.read_pickle(terms_file_ss8)
    terms_ss8 = terms_df_ss8['terms'].values.flatten()
    terms_dict_ss8 = {v: i for i, v in enumerate(terms_ss8)}
    terms_classes_ss8 = len(terms_dict_ss8)


    ############ step1: Diamond ##################
    # IC
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

    print(type(diamond_preds))




    ############ step2: check_model.pth

    ############ step2.1 aa

    test_df_aa = pd.read_pickle(test_data_file_aa)  # 这里需要test_data.pkl文件，重新进行类似于step3_Test.py的步骤，根据训练好的model再预测
    test_dataset_aa = ProteinGODataset(test_df_aa, terms_dict_aa, PROT_LETTER_len_aa, PROT_INDEX_aa)
    test_dataloader_aa = DataLoader(test_dataset_aa, batch_size=batch_size, shuffle=False)


    device_ids_aa = [0, 1]

    if isinstance(device_ids_aa, str):
        print('--- single GPU = ', device_ids_aa)
        device_aa = torch.device(device_ids_aa if torch.cuda.is_available() else "cpu")
        model_aa = Model1(terms_classes_aa, params_local_aa, kernels_tuple_aa, filters_list_aa)
        model_aa.to(device_aa)

    else:
        print('--- multi GPU = ', device_ids_aa)
        device_aa = torch.device("cuda:{}".format(device_ids_aa[0]) if torch.cuda.is_available() else "cpu")

        model_aa = Model1(terms_classes_aa, params_local_aa, kernels_tuple_aa, filters_list_aa)
        model_aa.to(device_aa)  #
        model_aa = nn.DataParallel(model_aa, device_ids=device_ids_aa)
        model_aa.to(device_aa)

    model_aa.load_state_dict(torch.load(model_file_aa))
    model_aa.to(device_aa)
    model_aa.eval()




    ###################### step2.2 ss8

    test_df_ss8 = pd.read_pickle(test_data_file_ss8)
    test_dataset_ss8 = ProteinGODataset(test_df_ss8, terms_dict_ss8, PROT_LETTER_len_ss8, PROT_INDEX_ss8)
    test_dataloader_ss8 = DataLoader(test_dataset_ss8, batch_size=batch_size, shuffle=False)

    device_ids_ss8 = [2, 3]

    if isinstance(device_ids_ss8, str):
        print('--- single GPU = ', device_ids_ss8)
        device_ss8 = torch.device(device_ids_ss8 if torch.cuda.is_available() else "cpu")
        model_ss8 = Model1(terms_classes_ss8, params_local_ss8, kernels_tuple_ss8, filters_list_ss8)
        model_ss8.to(device_ss8)

    else:
        print('--- multi GPU = ', device_ids_ss8)
        device_ss8 = torch.device("cuda:{}".format(device_ids_ss8[0]) if torch.cuda.is_available() else "cpu")

        model_ss8 = Model1(terms_classes_ss8, params_local_ss8, kernels_tuple_ss8, filters_list_ss8)  #
        model_ss8.to(device_ss8)
        model_ss8 = nn.DataParallel(model_ss8, device_ids=device_ids_ss8)
        model_ss8.to(device_ss8)

    model_ss8.load_state_dict(torch.load(model_file_ss8))
    model_ss8.to(device_ss8)
    model_ss8.eval()


    ########## step3 integration

    # (1 - alpha - beta) * diamond + alpha * preds_aa + beta * preds_ss8
    alphas = {NAMESPACES['bp']: 0, NAMESPACES['cc']: 0, NAMESPACES['mf']: 0}
    betas = {NAMESPACES['bp']: 0, NAMESPACES['cc']: 0, NAMESPACES['mf']: 0}

    if alpha == 'json' and beta == 'json':
        print('alpha/beta is from json, alpha/beta = ', alpha, beta)
        last_release_metadata = 'Alpha_last_release.json'
        with open(last_release_metadata, 'r') as f:
            print('Reading file from json')
            last_release_data = json.load(f)

            alpha = last_release_data["alphas"][ont]
            print('111111111', type(alpha))
            alphas[NAMESPACES[ont]] = alpha

            beta = last_release_data["betas"][ont]
            print('22222222', type(beta))
            betas[NAMESPACES[ont]] = beta

    else:
        print('alpha/beta is from click, alpha/beta = ', alpha, beta)
        print('type_alpha/beta = ', type(alpha), type(beta))
        alphas[NAMESPACES[ont]] = alpha
        betas[NAMESPACES[ont]] = beta

    print('lol')
    print('alphas, beta = ', alphas, betas)

    start_time = time.time()
    total_seq_aa = 0
    total_seq_ss8 = 0

    # output
    if ont == 'bp':
        out_file_ont = out_file_bp
    elif ont == 'cc':
        out_file_ont = out_file_cc
    elif ont == 'mf':
        out_file_ont = out_file_mf

    w = open(out_file_ont, 'wt')


    ################## deep pred_aa

    for prot_ids_aa, sequences_aa in read_fasta(in_file_aa, chunk_size):
        total_seq_aa += len(prot_ids_aa)
        deep_preds_aa = {}
        ids, data = get_data(sequences_aa, PROT_LETTER_len_aa, PROT_INDEX_aa)


        preds_list_aa = []
        with torch.no_grad():
            for inputs_onehot_aa, target_labels_aa in test_dataloader_aa:
                inputs_onehot_aa = inputs_onehot_aa.to(device_aa)

                outputs_aa = model_aa(inputs_onehot_aa)
                preds_list_aa.append(outputs_aa.cpu().numpy())
        preds_aa = np.concatenate(preds_list_aa)

        assert preds_aa.shape[1] == len(terms_aa)
        for i, j in enumerate(ids):
            prot_id_aa = prot_ids_aa[j]
            if prot_id_aa not in deep_preds_aa:
                deep_preds_aa[prot_id_aa] = {}
            for l in range(len(terms_aa)):
                if preds_aa[i, l] >= 0.01:
                    if terms_aa[l] not in deep_preds_aa[prot_id_aa]:
                        deep_preds_aa[prot_id_aa][terms_aa[l]] = preds_aa[i, l]
                    else:
                        deep_preds_aa[prot_id_aa][terms_aa[l]] = max(deep_preds_aa[prot_id_aa][terms_aa[l]], preds_aa[i, l])


        ################## deep pred_ss8 ######################

        for prot_ids_ss8, sequences_ss8 in read_fasta(in_file_ss8, chunk_size):
            total_seq_ss8 += len(prot_ids_ss8)
            deep_preds_ss8 = {}
            ids, data = get_data(sequences_ss8, PROT_LETTER_len_ss8, PROT_INDEX_ss8)

            preds_list_ss8 = []
            # labels_list_ss8 = []
            with torch.no_grad():
                for inputs_onehot_ss8, target_labels_ss8 in test_dataloader_ss8:
                    inputs_onehot_ss8 = inputs_onehot_ss8.to(device_ss8)

                    outputs_ss8 = model_ss8(inputs_onehot_ss8)
                    preds_list_ss8.append(outputs_ss8.cpu().numpy())
            preds_ss8 = np.concatenate(preds_list_ss8)


            assert preds_ss8.shape[1] == len(terms_ss8)
            for i, j in enumerate(ids):
                prot_id_ss8 = prot_ids_ss8[j]
                if prot_id_ss8 not in deep_preds_ss8:
                    deep_preds_ss8[prot_id_ss8] = {}
                for l in range(len(terms_ss8)):
                    if preds_ss8[i, l] >= 0.01:
                        if terms_ss8[l] not in deep_preds_ss8[prot_id_ss8]:
                            deep_preds_ss8[prot_id_ss8][terms_ss8[l]] = preds_ss8[i, l]
                        else:
                            deep_preds_ss8[prot_id_ss8][terms_ss8[l]] = max(deep_preds_ss8[prot_id_ss8][terms_ss8[l]], preds_ss8[i, l])


            for prot_id in prot_ids_aa:
                annots = {}
                if prot_id in diamond_preds:
                    for go_id, score in diamond_preds[prot_id].items():
                        annots[go_id] = score * (1 - alphas[go.get_namespace(go_id)])


                for go_id, score in deep_preds_aa[prot_id].items():
                    if go_id in annots:
                        annots[go_id] += alphas[go.get_namespace(go_id)] * score
                    else:
                        annots[go_id] = alphas[go.get_namespace(go_id)] * score

                for go_id, score in deep_preds_ss8[prot_id].items():
                    if go_id in annots:
                        annots[go_id] += betas[go.get_namespace(go_id)] * score
                    else:
                        annots[go_id] = betas[go.get_namespace(go_id)] * score



                # Propagate scores with ontology structure
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
    print('Total prediction time for %d sequences is %d' % (total_seq_aa, total_time))



######################## class/func

def to_onehot(seq, PROT_LETTER_len, PROT_INDEX, start=0):
    onehot = np.zeros((params_local_aa['MAXLEN'], PROT_LETTER_len + 1), dtype=np.int32)
    l = min(params_local_aa['MAXLEN'], len(seq))
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
        terms_classes = len(self.terms_dict)
        row = self.df.iloc[idx]
        seq = row.sequences
        onehot = to_onehot(seq, self.PROT_LETTER_len, self.PROT_INDEX, start=0)
        data_onehot = torch.zeros((params_local_aa['MAXLEN'], self.PROT_LETTER_len + 1), dtype=torch.float32)
        labels = torch.zeros(terms_classes, dtype=torch.float)
        data_onehot[:onehot.shape[0], :] = torch.from_numpy(onehot)

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

            pool = nn.MaxPool2d(kernel_size=(params_local_aa['MAXLEN'] - self.kernels_tuple[i][0] + 1, 1))

            flat = nn.Flatten()
            self.convs.append(conv)
            self.pools.append(pool)
            self.flats.append(flat)

        if params_local_aa['FC_depth'] > 0:
            self.fc_layers = nn.ModuleList()
            for i in range(params_local_aa['FC_depth']):
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


        if params_local_aa['FC_depth'] > 0:
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
