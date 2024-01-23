#!/usr/bin/env python

import numpy as np
import pandas as pd
import click as ck
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
import math
from utils import FUNC_DICT, Ontology, NAMESPACES
from joblib import Parallel, delayed
import json
import os
from step0_TrainTestSetting_global import path_base


@ck.command()
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='Data file with training features')
@ck.option('--predictions-file', '-pf', default='data/predictions.pkl', help='XX')
@ck.option('--terms-file', '-tf', default='data/terms_gominre_trxte.pkl', help='Data file with sequences and complete set of annotations')
@ck.option('--diamond-scores-file', '-dsf', default='data/diamond_aa.res', help='Diamond output')
@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--run-label-block', '-rbx', default='T', help='judge whether run block X')
@ck.option('--alpha-range', '-ar', default='25, 80, 1', type=str, help='xx')
@ck.option('--num-cpu-cores', '-ncc', default=60, type=int, help='xx')
@ck.option('--last-release-metadata', '-lrm', default='Alpha_last_release.json', help='xx')


def main(go_file, train_data_file, predictions_file, terms_file, diamond_scores_file, ont,
         run_label_block, alpha_range, num_cpu_cores, last_release_metadata):


    go_rels = Ontology(go_file, with_rels=True)
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()

    train_df = pd.read_pickle(train_data_file)
    test_df = pd.read_pickle(predictions_file)
    print("Length of test set: " + str(len(test_df)))

    annotations = train_df['prop_annotations'].values
    annotations = list(map(lambda x: set(x), annotations))
    test_annotations = test_df['prop_annotations'].values
    test_annotations = list(map(lambda x: set(x), test_annotations))
    go_rels.calculate_ic(annotations + test_annotations)


    ####################### Diamond IC ##############################
    # Print IC values of terms  # information content
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



    with open(last_release_metadata, 'r') as f:
        last_release_data = json.load(f)
    print('111', last_release_data)


    last_release_data['alphas'][ont] = find_alpha(ont, test_df, blast_preds, go_rels, terms, run_label_block, alpha_range, num_cpu_cores)
    print('2222', last_release_data)
    with open(last_release_metadata, 'w') as f:
        json.dump(last_release_data, f)





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



def eval_alphas(alpha, ont, test_df, blast_preds, go_rels, terms, run_label_block):
    deep_preds = []
    go_set = go_rels.get_namespace_terms(NAMESPACES[ont])
    go_set.remove(FUNC_DICT[ont])

    labels = test_df['prop_annotations'].values
    labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))

    if run_label_block == 'T':
        terms_set = set(terms)
        labels = [s.intersection(terms_set) for s in labels]



    alphas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}
    alphas[NAMESPACES[ont]] = alpha  #
    print('testing alpha =', alpha)




    for i, row in enumerate(test_df.itertuples()):
        # score = deepgo_preds * alpha + blast * (1-alpha)
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
        precisions.append(prec)
        recalls.append(rec)
        if fmax < fscore:
            fmax = fscore
            tmax = threshold
        if smin > s:
            smin = s

    precisions = np.array(precisions)
    recalls = np.array(recalls)
    sorted_index = np.argsort(recalls)
    recalls = recalls[sorted_index]
    precisions = precisions[sorted_index]
    aupr = np.trapz(precisions, recalls)


    if ont == 'mf':
        smin /= 15
    elif ont == 'bp':
        smin /= 40
    elif ont == 'cc':
        smin /= 10

    return alpha, -fmax


def find_alpha(ont, test_df, blast_preds, go_rels, terms, run_label_block, alpha_range, num_cpu_cores):

    extra = [ont, test_df, blast_preds, go_rels, terms, run_label_block]

    str_values = alpha_range.split(',')
    int_values = [int(value.strip()) for value in str_values]
    inputs = range(*int_values)
    print('inputs=', inputs)

    num_cores = num_cpu_cores
    results = Parallel(n_jobs=num_cores)(delayed(eval_alphas)(i / 100, *extra) for i in inputs)

    chosen = min(results, key=lambda x: x[1])
    print(chosen)

    return chosen[0]




if __name__ == '__main__':
    alphas = {'mf': 0, 'bp': 0, 'cc': 0}
    main()

