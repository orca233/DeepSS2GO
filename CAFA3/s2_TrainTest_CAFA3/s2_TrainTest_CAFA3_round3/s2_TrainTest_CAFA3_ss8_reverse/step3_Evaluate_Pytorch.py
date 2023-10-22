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
from matplotlib import pyplot as plt
import json
from step0_TrainTestSetting_local import *



# logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)




##################### func/class ##############################

def compute_roc(labels, preds):
    # Compute ROC curve and ROC area for each class
    fpr, tpr, _ = roc_curve(labels.flatten(), preds.flatten())
    roc_auc = auc(fpr, tpr)
    return roc_auc


def compute_mcc(labels, preds):
    # Compute ROC curve and ROC area for each class
    mcc = matthews_corrcoef(labels.flatten(), preds.flatten())
    return mcc


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
        tp = set(real_annots[i]).intersection(set(pred_annots[i]))  # true positive
        fp = pred_annots[i] - tp  # false positive
        fn = real_annots[i] - tp  # false negative
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
    # fscore, prec, rec, s, ru, mi, fps, fn

################# func/class done ######################################



######## pre_setting #############


train_data_file = 'data/train_data.pkl'  # train_data 只包含 train & valid
# test_data_file = 'data/test_data.pkl'  # test_data 只包含 test
# terms_file = 'data/terms.pkl'
predictions_file = 'data/predictions.pkl'

# terms_all_file = 'data/terms_all.pkl'  # 同时出现在train&test中，重复次数超过GOMinRepeat，
# terms_bp_file = 'data/terms_bp.pkl'  # terms_all.pkl中，属于bp部分
# terms_cc_file = 'data/terms_cc.pkl'
# terms_mf_file = 'data/terms_mf.pkl'
terms_file ='data/terms_%s.pkl' % params_local['onts']  # 給terms_file 指定对应的terms_bp/cc/mf_file文件夹


# ont = 'bp'  # 这个要放在step0_local.py里

go_rels = Ontology(params_local['go_file'], with_rels=True)
terms_df = pd.read_pickle(terms_file)
terms = terms_df['terms'].values.flatten()
# terms_dict = {v: i for i, v in enumerate(terms)}

train_df = pd.read_pickle(train_data_file)
test_df = pd.read_pickle(predictions_file)  # predictions.pkl的信息
print("Length of test set: " + str(len(test_df)))

annotations = train_df['prop_annotations'].values
annotations = list(map(lambda x: set(x), annotations))  # 变成set
test_annotations = test_df['prop_annotations'].values
test_annotations = list(map(lambda x: set(x), test_annotations))
go_rels.calculate_ic(annotations + test_annotations)

# Print IC values of terms
ics = {}
for term in terms:
    ics[term] = go_rels.get_ic(term)

prot_index = {}
for i, row in enumerate(train_df.itertuples()):
    prot_index[row.proteins] = i


# DeepGOPlus
subontology = ['bp', 'cc', 'mf']
for ont in subontology:
    print('\n ------ ont=%s' % ont)
    # 如果训练用的 params_local['onts']也是 bp/cc/mf中的一个，或者 训练所有all，则进行该 subontology的 evaluate。
    if params_local['onts'] == ont or params_local['onts'] == 'all':


        go_set = go_rels.get_namespace_terms(NAMESPACES[ont])  # 所有在当前属于：bp/cc/mf下的GO？
        go_set.remove(FUNC_DICT[ont])
        labels = test_df['prop_annotations'].values
        labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))


        # # block X # 这里添加一下，把labels的每一个元素{GO,GO,GO},{GO,GO},{GO}和terms的2088做交集
        # terms_set = set(terms)
        # labels = [s.intersection(terms_set) for s in labels]
        # # print(labels)
        # print('len(labels & terms_set) = ', len(labels))





        deep_preds = []
        for i, row in enumerate(test_df.itertuples()):  # 这是怎么个循环的。。。循环662遍，对不？
            annots_dict = {}
            for j, score in enumerate(row.preds):
                go_id = terms[j]
                if go_id in annots_dict:  # annots_dict 是从blast/diamond初始化来的
                    annots_dict[go_id] += score  # 这是怎么积分的？两个简单相加？
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
        for t in range(1, 101):  # the range in this loop has influence in the AUPR output
            threshold = t / 100.0
            preds = []
            for i, row in enumerate(test_df.itertuples()):
                annots = set()
                for go_id, score in deep_preds[i].items():  # 一共662行，逐行i看
                    if score >= threshold:
                        annots.add(go_id)
                new_annots = set()

                for go_id in annots:
                    new_annots |= go_rels.get_ancestors(go_id)  # 求并集
                preds.append(new_annots)

            # Filter classes
            preds = list(map(lambda x: set(filter(lambda y: y in go_set, x)), preds))

            fscore, prec, rec, s, ru, mi, fps, fns = evaluate_annotations(go_rels, labels, preds)
            avg_fp = sum(map(lambda x: len(x), fps)) / len(fps)
            avg_ic = sum(map(lambda x: sum(map(lambda go_id: go_rels.get_ic(go_id), x)), fps)) / len(fps)
            # print(f'{avg_fp} {avg_ic}')
            precisions.append(prec)
            recalls.append(rec)
            # print(f'Fscore: {fscore}, Precision: {prec}, Recall: {rec} S: {s}, RU: {ru}, MI: {mi} threshold: {threshold}')
            if fmax < fscore:
                fmax = fscore
                tmax = threshold
            if smin > s:
                smin = s
        print(f'Smin: {smin:0.3f}')
        print(f'Fmax: {fmax:0.3f}')
        precisions = np.array(precisions)
        recalls = np.array(recalls)
        sorted_index = np.argsort(recalls)
        recalls = recalls[sorted_index]
        precisions = precisions[sorted_index]
        aupr = np.trapz(precisions, recalls)
        print(f'AUPR: {aupr:0.3f}')

        ### FS 保存Smin,Fmax,AUPR ###
        with open(f'results/Fmax_AUPR_Smin_{ont}.txt', 'w') as f:
            f.write("Length of test set: %s  " % len(test_df))
            f.write(f'threshold: {tmax} \n')
            f.write(f'Smin: {smin:0.3f} \n')
            f.write(f'Fmax: {fmax:0.3f} \n')
            f.write(f'AUPR: {aupr:0.3f} \n')

        lw = 2
        # plt.figure()
        # plt.plot(recalls, precisions, color='darkorange', lw=lw, label=f'AUPR curve (area = {aupr:0.2f})')
        # plt.xlim([0.0, 1.0])
        # plt.ylim([0.0, 1.05])
        # plt.xlabel('Recall')
        # plt.ylabel('Precision')
        # plt.title('Area Under the Precision-Recall curve')
        # plt.legend(loc="lower right")
        df = pd.DataFrame({'precisions': precisions, 'recalls': recalls})
        df.to_pickle(f'results/PR_{ont}.pkl')  # 保存 pkl
        df.to_csv(f'results/PR_{ont}.csv', index=False)  # 同时 保存一个csv


    else:
        print(f'HOLY SHIT, this subontology ({ont}) is not trained')


'''

1）如何使用多个cpu并行？？？

2）chatGPT 如何计算Smin， Fmax, AUPR，给出Python代码

from sklearn.metrics import precision_recall_curve

# 假设y_true为真实标签，y_score为模型预测的概率值
# 计算精度和召回率
precision, recall, _ = precision_recall_curve(y_true, y_score)

# 计算Smin和Fmax
smin = 1 - min(precision)
fmax = max(2 * precision * recall / (precision + recall))

# 计算AUPR
from sklearn.metrics import auc
aupr = auc(recall, precision)

'''



