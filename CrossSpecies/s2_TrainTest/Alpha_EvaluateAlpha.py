#!/usr/bin/env python

import numpy as np
import pandas as pd
import click as ck
# import logging
import math
from utils import FUNC_DICT, Ontology, NAMESPACES
from matplotlib import pyplot as plt
import json
# from step0_TrainTestSetting_global import *
from step0_TrainTestSetting_local import *
from step0_TrainTestSetting_global import path_base


@ck.command()
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='Data file with training features')
@ck.option('--test-data-file', '-tsdf', default='data/predictions.pkl', help='Test data file')

@ck.option('--terms-file', '-tf', default='data/terms_gominre_trxte.pkl', help='Data file with sequences and complete set of annotations')  #
@ck.option('--diamond-scores-file', '-dsf', default='data/diamond_aa.res', help='Diamond output')
@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--alpha', '-a', default='json', help='alpha = json(with quote) or 0-1(without quote, eg 0.3 float)')
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')
@ck.option('--run-label-block', '-rbx', default='T', help='judge whether run block X')



def main(train_data_file, test_data_file, terms_file, diamond_scores_file, ont, alpha, go_file, run_label_block):
    go_rels = Ontology(go_file, with_rels=True)
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()

    train_df = pd.read_pickle(train_data_file)
    test_df = pd.read_pickle(test_data_file)
    print("Length of test set: " + str(len(test_df)))

    annotations = train_df['prop_annotations'].values
    annotations = list(map(lambda x: set(x), annotations))
    test_annotations = test_df['prop_annotations'].values
    test_annotations = list(map(lambda x: set(x), test_annotations))
    go_rels.calculate_ic(annotations + test_annotations)




    ########################## Diamond IC ############################
    # Print IC values of terms
    ics = {}
    for term in terms:
        ics[term] = go_rels.get_ic(term)

    prot_index = {}
    for i, row in enumerate(train_df.itertuples()):
        prot_index[row.proteins] = i

    # BLAST Similarity (Diamond)
    diamond_scores = {}
    with open(diamond_scores_file) as f:
        for line in f:
            it = line.strip().split()
            if it[0] not in diamond_scores:
                diamond_scores[it[0]] = {}
            diamond_scores[it[0]][it[1]] = float(it[2])

    blast_preds = []
    for i, row in enumerate(test_df.itertuples()):
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



    # DeepGOPlus
    deep_preds = []
    go_set = go_rels.get_namespace_terms(NAMESPACES[ont])
    go_set.remove(FUNC_DICT[ont])
    labels = test_df['prop_annotations'].values
    labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))

    if run_label_block == 'T':
        terms_set = set(terms)
        labels = [s.intersection(terms_set) for s in labels]

    alphas = {NAMESPACES['bp']: 0, NAMESPACES['cc']: 0, NAMESPACES['mf']: 0}  # 重新初始化？


    # deep_preds
    if alpha == 'json':
        print('alpha is from json, alpha = ', alpha, 'alpha_type = ', type(alpha))

        last_release_metadata = 'Alpha_last_release.json'
        with open(last_release_metadata, 'r') as f:
            print('Reading file from json')
            last_release_data = json.load(f)
            alpha = last_release_data["alphas"][ont]
            print('type_alpha = ', type(alpha))
            alphas[NAMESPACES[ont]] = alpha
    else:
        print('alpha is from click, alpha = ', alpha)
        print('type_alpha = ', type(alpha))
        alpha = float(alpha)
        alphas[NAMESPACES[ont]] = alpha
        print('updated_type_alpha = ', type(alphas[NAMESPACES[ont]]))

    print('lol')
    print('alphas=', alphas)



    # (1 - alpha - beta) * diamond + alpha * preds_aa + beta * preds_ss8
    for i, row in enumerate(test_df.itertuples()):
        annots_dict = blast_preds[i].copy()
        for go_id in annots_dict:
            annots_dict[go_id] *= 1 - alphas[go_rels.get_namespace(go_id)]

        for j, score in enumerate(row.preds):
            go_id = terms[j]
            score *= alphas[go_rels.get_namespace(go_id)]
            if go_id in annots_dict:
                annots_dict[go_id] += score
            else:
                annots_dict[go_id] = score
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
        for i, row in enumerate(test_df.itertuples()):
            annots = set()
            for go_id, score in deep_preds[i].items():
                if score >= threshold:
                    annots.add(go_id)

            new_annots = set()
            for go_id in annots:
                new_annots |= go_rels.get_ancestors(go_id)
            preds.append(new_annots)

        # Filter classes
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

    with open(f'results_alpha/Fmax_AUPR_Smin_{ont}_{alpha:0.2f}.txt', 'w') as f:
        f.write(f'Length of test set: {len(test_df)} \n')
        f.write(f'threshold: {tmax} \n')
        f.write(f'alpha={alpha:0.2f} \n')
        f.write(f'Smin: {smin:0.3f} \n')
        f.write(f'Fmax: {fmax:0.3f} \n')
        f.write(f'AUPR: {aupr:0.3f} \n')

    df = pd.DataFrame({'precisions': precisions, 'recalls': recalls})
    df.to_pickle(f'results_alpha/PR_{ont}_{alpha:0.2f}.pkl')



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