#!/usr/bin/env python


import numpy as np
import pandas as pd
import click as ck
from sklearn.metrics import classification_report
from sklearn.metrics.pairwise import cosine_similarity
import sys
from collections import deque
import time
import logging
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
from scipy.spatial import distance
from scipy import sparse
import math
from utils import FUNC_DICT, Ontology, NAMESPACES
from Alpha_FindAlpha import evaluate_annotations  # compute_mcc, compute_roc
from matplotlib import pyplot as plt
from joblib import Parallel, delayed
import multiprocessing
import json



logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


@ck.command()
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')
@ck.option('--train-data-file', '-trdf', default='data/train_data_aa.pkl', help='Data file with training features')
@ck.option('--test-data-file-aa-ss8', '-tsdf', default='data/predictions_aa_ss8.pkl', help='Test data file')
@ck.option('--terms-file-aa', '-tf', default='./data/terms_gominre_trxte_aa.pkl', help='set of annotations')
@ck.option('--terms-file-ss8', '-tf', default='./data/terms_gominre_trxte_ss8.pkl', help='set of annotations')
@ck.option('--diamond-scores-file', '-dsf', default='data/diamond_aa.res', help='Diamond output')
@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--run-label-block', '-rbx', default='T', help='judge whether run block X')
@ck.option('--alpha-range', '-ar', default='0, 101, 10', type=str, help='xx')
@ck.option('--beta-range', '-br', default='0, 101, 10', type=str, help='xx')
@ck.option('--num-cpu-cores', '-ncc', default=80, type=int, help='xx')
@ck.option('--last-release-metadata', '-lrm', default='Alpha_last_release.json', help='xx')


def main(go_file, train_data_file, test_data_file_aa_ss8, terms_file_aa, terms_file_ss8, diamond_scores_file, ont,
         run_label_block, alpha_range, beta_range, num_cpu_cores, last_release_metadata):


    go_rels = Ontology(go_file, with_rels=True)
    terms_df_aa = pd.read_pickle(terms_file_aa)
    terms_aa = terms_df_aa['terms'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms_aa)}


    terms_df_ss8 = pd.read_pickle(terms_file_ss8)
    terms_ss8 = terms_df_ss8['terms'].values.flatten()

    train_df = pd.read_pickle(train_data_file)
    test_df_aa_ss8 = pd.read_pickle(test_data_file_aa_ss8)
    print("Length of test set: " + str(len(test_df_aa_ss8)))

    annotations = train_df['prop_annotations'].values
    annotations = list(map(lambda x: set(x), annotations))

    test_annotations_aa_ss8 = test_df_aa_ss8['prop_annotations'].values
    test_annotations_aa_ss8 = list(map(lambda x: set(x), test_annotations_aa_ss8))

    go_rels.calculate_ic(annotations + test_annotations_aa_ss8)


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
            diamond_scores[it[0]][it[1]] = float(it[2])



    blast_preds = []
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



    with open(last_release_metadata, 'r') as f:
        last_release_data = json.load(f)
    print('111', last_release_data)

    last_release_data['alphas'][ont], last_release_data['betas'][ont] = find_alpha_beta(ont, test_df_aa_ss8, blast_preds, go_rels, terms_aa, terms_ss8, run_label_block, alpha_range, beta_range, num_cpu_cores)
    print('2222', last_release_data)
    with open(last_release_metadata, 'w') as f:
        json.dump(last_release_data, f)



def eval_alphas_betas(alpha, beta, ont, test_df_aa_ss8, blast_preds, go_rels, terms_aa, terms_ss8, run_label_block):
    print('alpha, beta = ', alpha, beta)

    deep_preds = []
    go_set = go_rels.get_namespace_terms(NAMESPACES[ont])
    go_set.remove(FUNC_DICT[ont])

    labels = test_df_aa_ss8['prop_annotations'].values
    labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))



    if run_label_block == 'T':
        terms_set = set(terms_aa)
        labels = [s.intersection(terms_set) for s in labels]


    alphas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}
    alphas[NAMESPACES[ont]] = alpha

    betas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}
    betas[NAMESPACES[ont]] = beta




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

    return alpha, beta, -fmax




def find_alpha_beta(ont, test_df_aa_ss8, blast_preds, go_rels, terms_aa, terms_ss8, run_label_block, alpha_range, beta_range, num_cpu_cores):

    extra = [ont, test_df_aa_ss8, blast_preds, go_rels, terms_aa, terms_ss8, run_label_block]
    str_values_alpha = alpha_range.split(',')
    int_values_alpha = [int(value.strip()) for value in str_values_alpha]
    inputs_alpha = range(*int_values_alpha)  # range(25, 80, 1)

    str_values_beta = beta_range.split(',')
    int_values_beta = [int(value.strip()) for value in str_values_beta]
    inputs_beta = range(*int_values_beta)  # range(25, 80, 1)

    print('inputs_alpha/beta = ', inputs_alpha, inputs_beta)


    num_cores = num_cpu_cores

    results = Parallel(n_jobs=num_cores)(
        delayed(eval_alphas_betas)(i / 100, j / 100, *extra)
        for i in inputs_alpha
        for j in inputs_beta
        if i + j <= 100

    )


    chosen = min(results, key=lambda x: x[2])
    print('chosen')
    print(chosen)

    return chosen[0], chosen[1]




if __name__ == '__main__':
    alphas = {'mf': 0, 'bp': 0, 'cc': 0}
    main()
