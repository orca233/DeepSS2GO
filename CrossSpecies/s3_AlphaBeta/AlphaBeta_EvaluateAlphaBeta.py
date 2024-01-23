

#!/usr/bin/env python
import os

import numpy as np
import pandas as pd
import click as ck
# import logging
import math
from utils import FUNC_DICT, Ontology, NAMESPACES
from matplotlib import pyplot as plt
import json
# from step0_TrainTestSetting_global import *
# from step0_TrainTestSetting_local_aa import *
from step0_TrainTestSetting_local_aa import params_local as params_local_aa
from step0_TrainTestSetting_local_ss8 import params_local as params_local_ss8
from step0_TrainTestSetting_global import path_base



os.system('mkdir -p results_alphabeta')

@ck.command()
@ck.option('--train-data-file', '-trdf', default='data/train_data_aa.pkl', help='Data file with training features')
@ck.option('--test-data-file-aa-ss8', '-tsdf', default='data/predictions_aa_ss8.pkl', help='Test data file')
@ck.option('--terms-file-aa', '-tf', default='data/terms_gominre_trxte_aa.pkl', help='Data file with sequences and complete set of annotations')
@ck.option('--terms-file-ss8', '-tf', default='data/terms_gominre_trxte_ss8.pkl', help='Data file with sequences and complete set of annotations')
@ck.option('--diamond-scores-file', '-dsf', default='data/diamond_aa.res', help='Diamond output')
@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--alpha', '-a', default='json', help='alpha = json(with quote) or 0-1(without quote, eg 0.3 float)')
@ck.option('--beta', '-b', default='json', help='beta = json(with quote) or 0-1(without quote, eg 0.3 float)')
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
@ck.option('--run-label-block', '-rbx', default='T', help='judge whether run block X')



def main(train_data_file, test_data_file_aa_ss8, terms_file_aa, terms_file_ss8, diamond_scores_file, ont, alpha, beta,
         go_file, run_label_block):

    go_rels = Ontology(go_file, with_rels=True)
    terms_df_aa = pd.read_pickle(terms_file_aa)
    terms_aa = terms_df_aa['terms'].values.flatten()


    terms_df_ss8 = pd.read_pickle(terms_file_ss8)
    terms_ss8 = terms_df_ss8['terms'].values.flatten()


    train_df = pd.read_pickle(train_data_file)
    test_df_aa_ss8 = pd.read_pickle(test_data_file_aa_ss8)
    print("Length of test set: " + str(len(test_df_aa_ss8)))

    annotations = train_df['prop_annotations'].values
    annotations = list(map(lambda x: set(x), annotations))
    test_annotations = test_df_aa_ss8['prop_annotations'].values
    test_annotations = list(map(lambda x: set(x), test_annotations))
    go_rels.calculate_ic(annotations + test_annotations)


    ics = {}
    for term in terms_aa:
        ics[term] = go_rels.get_ic(term)

    prot_index = {}
    for i, row in enumerate(train_df.itertuples()):
        prot_index[row.proteins] = i

    diamond_scores = {}
    with open(diamond_scores_file) as f:
        for line in f:
            it = line.strip().split()
            if it[0] not in diamond_scores:
                diamond_scores[it[0]] = {}
            diamond_scores[it[0]][it[1]] = float(it[2])  #

    blast_preds = []
    # print('Diamond preds')
    for i, row in enumerate(test_df_aa_ss8.itertuples()):
        annots = {}
        prot_id = row.proteins
        # BlastKNN
        if prot_id in diamond_scores:
            sim_prots = diamond_scores[prot_id]
            allgos = set()
            total_score = 0.0
            for p_id, score in sim_prots.items():
                allgos |= annotations[prot_index[p_id]]
                total_score += score
            allgos = list(sorted(allgos))
            sim = np.zeros(len(allgos), dtype=np.float32)
            for j, go_id in enumerate(allgos):
                s = 0.0
                for p_id, score in sim_prots.items():
                    if go_id in annotations[prot_index[p_id]]:
                        s += score
                sim[j] = s / total_score
            ind = np.argsort(-sim)
            for go_id, score in zip(allgos, sim):
                annots[go_id] = score
        blast_preds.append(annots)



    deep_preds = []
    go_set = go_rels.get_namespace_terms(NAMESPACES[ont])
    go_set.remove(FUNC_DICT[ont])

    labels = test_df_aa_ss8['prop_annotations'].values
    labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))

    if run_label_block == 'T':
        terms_set = set(terms_aa)
        labels = [s.intersection(terms_set) for s in labels]


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

    for i, row in enumerate(test_df_aa_ss8.itertuples()):

        annots_dict = blast_preds[i].copy()
        for go_id in annots_dict:
            annots_dict[go_id] *= 1 - alphas[go_rels.get_namespace(go_id)] - betas[go_rels.get_namespace(go_id)]

        for j, score_aa in enumerate(row.preds_aa):
            go_id = terms_aa[j]
            score_aa *= alphas[go_rels.get_namespace(go_id)]

            if go_id in annots_dict:
                annots_dict[go_id] += score_aa
            else:
                annots_dict[go_id] = score_aa
                
                
        for k, score_ss8 in enumerate(row.preds_ss8):
            go_id = terms_ss8[k]
            score_ss8 *= betas[go_rels.get_namespace(go_id)]
            if go_id in annots_dict:
                annots_dict[go_id] += score_ss8
            else:
                annots_dict[go_id] = score_ss8



        deep_preds.append(annots_dict)


    fmax = 0.0
    tmax = 0.0
    precisions = []
    recalls = []
    smin = 1000000.0
    rus = []
    mis = []
    for t in range(1, 101):
        threshold = t / 100.0
        preds = []
        for i, row in enumerate(test_df_aa_ss8.itertuples()):
            annots = set()
            for go_id, score in deep_preds[i].items():
                if score >= threshold:
                    annots.add(go_id)

            new_annots = set()
            for go_id in annots:
                new_annots |= go_rels.get_ancestors(go_id)
            preds.append(new_annots)

        preds = list(map(lambda x: set(filter(lambda y: y in go_set, x)), preds))

        fscore, prec, rec, s, ru, mi, fps, fns = evaluate_annotations(go_rels, labels, preds)
        avg_fp = sum(map(lambda x: len(x), fps)) / len(fps)
        avg_ic = sum(map(lambda x: sum(map(lambda go_id: go_rels.get_ic(go_id), x)), fps)) / len(fps)
        # print(f'{avg_fp} {avg_ic}')
        precisions.append(prec)
        recalls.append(rec)
        if fmax < fscore:
            fmax = fscore
            tmax = threshold
        if smin > s:
            smin = s

    print(f'threshold: {tmax}')
    print(f'Smin: {smin:0.3f}')
    print(f'Fmax: {fmax:0.3f}')
    precisions = np.array(precisions)
    recalls = np.array(recalls)
    sorted_index = np.argsort(recalls)
    recalls = recalls[sorted_index]
    precisions = precisions[sorted_index]
    aupr = np.trapz(precisions, recalls)
    print(f'AUPR: {aupr:0.3f}')

    with open(f'results_alphabeta/Fmax_AUPR_Smin_{ont}_alpha{alpha:0.2f}_beta{beta:0.2f}.txt', 'w') as f:
        f.write(f'Length of test set: {len(test_df_aa_ss8)} \n')
        f.write(f'threshold: {tmax} \n')
        f.write(f'alpha={alpha:0.2f} \n')
        f.write(f'Smin: {smin:0.3f} \n')
        f.write(f'Fmax: {fmax:0.3f} \n')
        f.write(f'AUPR: {aupr:0.3f} \n')

    df = pd.DataFrame({'precisions': precisions, 'recalls': recalls})
    df.to_pickle(f'results_alphabeta/PR_{ont}_alpha{alpha:0.2f}_beta{beta:0.2f}.pkl')



def evaluate_annotations(go, real_annots, pred_annots):
    total = 0
    p = 0.0
    r = 0.0
    p_total = 0
    ru = 0.0
    mi = 0.0
    fps = []
    fns = []
    for i in range(len(real_annots)):
        if len(real_annots[i]) == 0:
            continue
        tp = set(real_annots[i]).intersection(set(pred_annots[i]))
        fp = pred_annots[i] - tp
        fn = real_annots[i] - tp
        for go_id in fp:
            mi += go.get_ic(go_id)
        for go_id in fn:
            ru += go.get_ic(go_id)
        fps.append(fp)
        fns.append(fn)
        tpn = len(tp)
        fpn = len(fp)
        fnn = len(fn)
        total += 1
        recall = tpn / (1.0 * (tpn + fnn))
        r += recall
        if len(pred_annots[i]) > 0:
            p_total += 1
            precision = tpn / (1.0 * (tpn + fpn))
            p += precision
    ru /= total
    mi /= total
    r /= total
    if p_total > 0:
        p /= p_total
    f = 0.0
    if p + r > 0:
        f = 2 * p * r / (p + r)
    s = math.sqrt(ru * ru + mi * mi)
    return f, p, r, s, ru, mi, fps, fns


if __name__ == '__main__':
    main()