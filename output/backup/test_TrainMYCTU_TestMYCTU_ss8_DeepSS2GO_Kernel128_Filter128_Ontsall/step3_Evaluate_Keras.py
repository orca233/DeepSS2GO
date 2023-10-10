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
# from step0_TrainTestSetting import *


import configparser
config_local = configparser.ConfigParser()  # 创建 ConfigParser 对象
config_local.read('config_local.ini')  # 读取配置文件
# read param from config_local.ini
# 循环读取每个配置项并添加到字典中
config_local_dict = {}
for section in config_local.sections():
    section_dict = {}
    for option in config_local.options(section):
        section_dict[option] = config_local.get(section, option)
    config_local_dict[section] = section_dict


logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

@ck.command()
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='Data file with training features')
@ck.option('--test-data-file', '-tsdf', default='data/predictions.pkl', help='Test data file')
@ck.option('--terms-file', '-tf', default='data/terms_all.pkl', help='Data file with sequences and complete set of annotations')
# @ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--subontology', '-sot', default='???', help='GO subontology (bp, mf, cc)')
@ck.option('--go-file', '-gf', default='go.obo', help='Gene Ontology file in OBO Format')  # FS 添加

def main(train_data_file, test_data_file, terms_file, subontology, go_file):  # FS 添加go_file
    go_file = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO_Pytorch/pub_data/go.obo'
    go_rels = Ontology(go_file, with_rels=True)
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()
    # terms_dict = {v: i for i, v in enumerate(terms)}

    train_df = pd.read_pickle(train_data_file)
    test_df = pd.read_pickle(test_data_file)  # predictions.pkl的信息
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

    subontology = ['bp', 'cc', 'mf']
    # DeepGOPlus
    for ont in subontology:
        print('\n ------ ont=%s--------' % ont)
        go_set = go_rels.get_namespace_terms(NAMESPACES[ont])  # 所有在当前'mf'下的GO？
        go_set.remove(FUNC_DICT[ont])
        labels = test_df['prop_annotations'].values  # 这一步有问题，会影响nf，所以后面要加一个  block X
        labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))

        # labels = [{'GO:0042221', 'GO:1902635', 'GO:0009719', 'GO:0006644', ..},{},{}]
        # len=662


        # block X # 这里添加一下，把labels的每一个元素{},{},{}和terms的2088做交集
        terms_set = set(terms)
        labels = [s.intersection(terms_set) for s in labels]
        print('-----new labels -----')
        # print(labels)
        print('---len----')
        print(len(labels))


        deep_preds = []
        for i, row in enumerate(test_df.itertuples()):  # 这是怎么个循环的。。。循环662遍，横向提取每一行元素组成tuple
            '''
            len(row.prop_annotations) = 215
            len(row.preds) = 2088  对每个蛋白质，都进行2088个所有的terms进行预测。
            
            i, row 第0行
            0 Pandas(Index=49935, index=314240, proteins='PI42B_HUMAN', accessions='P78356; Q5U0E8; Q8TBP2;', 
            sequences='MSSNCTSTTAVAVAP...' annotations=['GO:0005776|IMP', 'GO:0005829|TAS', 'GO:0005789|IEA', ...], 
            interpros=['IPR023610', 'IPR002498', 'IPR027484'], orgs='9606', 
            exp_annotations=['GO:0005776', 'GO:0005829', ...], 
            prop_annotations=['GO:0008150', 'GO:0008641', 'GO:0023057', 'GO:0046627', 'GO:0045017',...], ----这个len=215
            cafa_target=True, labels=array([0, 0, 0, ..., 0, 0, 0], dtype=int32), 
            preds=array([0.13251728, 0.04844077, 0.03128296, ..., 0.04268779, 0.00633362, 0.00904174], dtype=float32))

            
            
            i, row:  最后一行
            661 Pandas(Index=16393, index=94312, proteins='DAZP2_HUMAN', 
            accessions='Q15038; A8K254; B4DDT5; B4E1G3; C9JA96; C9JP84; E9PB45; F8VU62;', 
            sequences='MNSKGQYPTQPTYPVQ...', annotations=['GO:0005737|IDA', 'GO:0010494|IDA', 'GO:0016604|IDA', ...], 
            interpros=['IPR022730'], orgs='9606', exp_annotations=['GO:0005737', 'GO:0010494', 'GO:0016604', ...], 
            prop_annotations=['GO:0008150', 'GO:0031648', 'GO:0044389', 'GO:0003674',...], cafa_target=True, 
            labels=array([0, 0, 0, ..., 0, 0, 0], dtype=int32), 
            preds=array([0.15103509, 0.00911674, 0.0355399 , ...], dtype=float32))

            '''

            annots_dict = {}
            for j, score in enumerate(row.preds):
                go_id = terms[j]
                if go_id in annots_dict:  # annots_dict 是从blast/diamond初始化来的
                    annots_dict[go_id] += score  # 这是怎么积分的？两个简单相加？
                else:
                    annots_dict[go_id] = score
            deep_preds.append(annots_dict)

        # deep_preds = [{'GO:0044281': 0.13251728, 'GO:0071824': 0.048440773, ...}, {}, {}...]  len=662

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
                # 每一个i，是一个蛋白质，在该蛋白质内部逐一看每个预测的GO term  # 一共662行，逐行i看
                annots = set()
                for go_id, score in deep_preds[i].items():
                    if score >= threshold:  # 把所有大于threshold的预测值的 GOterm 放到 annots 集合内
                        annots.add(go_id)

                new_annots = set()
                for go_id in annots:
                    new_annots |= go_rels.get_ancestors(go_id)  # 把祖先扒出来，求并集
                preds.append(new_annots)
            # len(preds)=662, preds = [{'GO:0022402', 'GO:0060255', 'GO:0010467',..},{},{}...]


            # Filter classes
            preds = list(map(lambda x: set(filter(lambda y: y in go_set, x)), preds))  # 过滤，判断每个GO是否在go_set，go_set已经局限于三大类之一比如bp
            # len(preds)=662, preds = [{'GO:0022402', 'GO:0060255', 'GO:0010467',..},{},{}...]

            # 这是对于一个threshold的计算
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

        plt.figure()
        lw = 2
        plt.plot(recalls, precisions, color='darkorange', lw=lw, label=f'AUPR curve (area = {aupr:0.2f})')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title('Area Under the Precision-Recall curve')
        plt.legend(loc="lower right")
        df = pd.DataFrame({'precisions': precisions, 'recalls': recalls})
        df.to_pickle(f'results/PR_{ont}_.pkl')


def compute_roc(labels, preds):
    # Compute ROC curve and ROC area for each class
    fpr, tpr, _ = roc_curve(labels.flatten(), preds.flatten())
    roc_auc = auc(fpr, tpr)
    return roc_auc


def compute_mcc(labels, preds):
    # Compute ROC curve and ROC area for each class
    mcc = matthews_corrcoef(labels.flatten(), preds.flatten())
    return mcc


def evaluate_annotations(go, real_annots, pred_annots):  # real_annots=labels, pred_annots=preds
    total = 0
    p = 0.0
    r = 0.0
    p_total = 0
    ru = 0.0
    mi = 0.0
    fps = []
    fns = []
    for i in range(len(real_annots)):  # len=662，      i是蛋白质index
        if len(real_annots[i]) == 0:
            continue

        # 这是对于一个threshold下，一个蛋白质的 labels & preds 求tp/fp/fn
        tp = set(real_annots[i]).intersection(set(pred_annots[i]))  # true positive
        fp = pred_annots[i] - tp  # false positive
        fn = real_annots[i] - tp  # false negative
        for go_id in fp:
            mi += go.get_ic(go_id)  # 这啥 蛋白质功能预测中的不匹配信息（Mismatch Information）
        for go_id in fn:
            ru += go.get_ic(go_id)  # 这啥  蛋白质功能预测中的不确定性（Resnik Uncertainty）

        # 这两个append是累积用的，把所有test的662条蛋白质清算加到一个list[]里
        fps.append(fp)
        fns.append(fn)

        # 这里计算的还是 蛋白质i 的tp/fp/fn 长度（个数）
        tpn = len(tp)
        fpn = len(fp)
        fnn = len(fn)

        total += 1
        recall = tpn / (1.0 * (tpn + fnn))

        r += recall  # r是累积 recall，后面初一个数求平均
        if len(pred_annots[i]) > 0:  # 对于 蛋白质i，有预测出新的prediction（比如针对bp）
            p_total += 1  # positive total个数
            precision = tpn / (1.0 * (tpn + fpn))  # precision是对于 蛋白质i的 计算
            p += precision  # p是累积，后面求平均

    # 注意：r&p不能初一同一个total，r是初一所有的total个数（why???），p是除以有阳的个数
    # 因为如果没预测出来任何GO term，也是有recall的。但没有预测出任何GO term，tp+fp=0不能做分母
    ru /= total
    mi /= total
    mi /= total
    r /= total

    if p_total > 0:
        p /= p_total  # p是累积，求对于所有“有预测结果”的平均

    f = 0.0  # fscore, beta=1
    if p + r > 0:
        f = 2 * p * r / (p + r)

    s = math.sqrt(ru * ru + mi * mi)  # 平方根？

    return f, p, r, s, ru, mi, fps, fns
    # fscore, prec, rec, s, ru, mi, fps, fn


if __name__ == '__main__':
    main()



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



