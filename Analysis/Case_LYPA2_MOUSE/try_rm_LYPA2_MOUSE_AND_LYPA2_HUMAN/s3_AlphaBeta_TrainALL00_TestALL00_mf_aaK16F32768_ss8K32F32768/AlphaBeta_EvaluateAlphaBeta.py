'''
from DeepGOPlus: evaluate_deepgoplus.py

相比于Alpha_FindAlpha.py，（从25-80逐一过一遍alpha，每一个alpha又过 t(threshold)0-100）
这个py文件只针对上述筛选出的某一个alpha进行，2min
'''


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


# logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

os.system('mkdir -p results_alphabeta')

@ck.command()
@ck.option('--train-data-file', '-trdf', default='data/train_data_aa.pkl', help='Data file with training features')
@ck.option('--test-data-file-aa-ss8', '-tsdf', default='data/predictions_aa_ss8.pkl', help='Test data file')  # test_data_file_aa_ss8 只有一个 prop_annotations, 两个 preds_aa & preds_ss8
# @ck.option('--predictions-file', '-pf', default='data/predictions.pkl', help='XX')

@ck.option('--terms-file-aa', '-tf', default='data/terms_gominre_trxte_aa.pkl', help='Data file with sequences and complete set of annotations')  # original: data/terms_all.pkl
@ck.option('--terms-file-ss8', '-tf', default='data/terms_gominre_trxte_ss8.pkl', help='Data file with sequences and complete set of annotations')  # original: data/terms_all.pkl

@ck.option('--diamond-scores-file', '-dsf', default='data/diamond_aa.res', help='Diamond output')
@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--alpha', '-a', default='json', help='alpha = json(with quote) or 0-1(without quote, eg 0.3 float)')  # 如果alpha='json'，则采用json数据，否则alpha=数字，或外来click引入
@ck.option('--beta', '-b', default='json', help='beta = json(with quote) or 0-1(without quote, eg 0.3 float)')  # 如果beta='json'，则采用json数据，否则beta=数字，或外来click引入

@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
# @ck.option('--go-file', '-gf', default=path_base + 'pub_data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
# @ck.option('--go-file', '-gf', default='../../pub_data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
@ck.option('--run-label-block', '-rbx', default='T', help='judge whether run block X')



def main(train_data_file, test_data_file_aa_ss8, terms_file_aa, terms_file_ss8, diamond_scores_file, ont, alpha, beta,
         go_file, run_label_block):  # FS 添加go_file
    # # 从 last_release_metadata 文件中获取alpha ###
    # last_release_metadata = 'Alpha_last_release.json'
    # with open(last_release_metadata, 'r') as f:
    #     last_release_data = json.load(f)
    #     alpha = last_release_data['alphas'][ont]  # # original 这里更新了alpha所以用不到ck中的引入信息
    #
    # print('alpha initial=', alpha)

    go_rels = Ontology(go_file, with_rels=True)
    terms_df_aa = pd.read_pickle(terms_file_aa)
    terms_aa = terms_df_aa['terms'].values.flatten()
    # terms_dict = {v: i for i, v in enumerate(terms)}  # 没用上？

    # FS:
    terms_df_ss8 = pd.read_pickle(terms_file_ss8)  # original: data/terms_all_ss8.pkl
    terms_ss8 = terms_df_ss8['terms'].values.flatten()


    train_df = pd.read_pickle(train_data_file)
    test_df_aa_ss8 = pd.read_pickle(test_data_file_aa_ss8)  # predictions.pkl的信息, original: test_data_file
    print("Length of test set: " + str(len(test_df_aa_ss8)))

    annotations = train_df['prop_annotations'].values
    annotations = list(map(lambda x: set(x), annotations))
    test_annotations = test_df_aa_ss8['prop_annotations'].values
    test_annotations = list(map(lambda x: set(x), test_annotations))
    go_rels.calculate_ic(annotations + test_annotations)

    #####################################################################
    ############################## Diamond ##############################
    #####################################################################
        # Print IC values of terms
    ics = {}
    for term in terms_aa:
        ics[term] = go_rels.get_ic(term)

    prot_index = {}
    for i, row in enumerate(train_df.itertuples()):
        prot_index[row.proteins] = i

    # BLAST Similarity (Diamond)  拿到Diamond-blast的结果
    diamond_scores = {}
    with open(diamond_scores_file) as f:
        for line in f:
            it = line.strip().split()
            if it[0] not in diamond_scores:
                diamond_scores[it[0]] = {}
            diamond_scores[it[0]][it[1]] = float(it[2])  # [2]这个是 diamond.res 第三列，也就是分数

    blast_preds = []
    # print('Diamond preds')
    for i, row in enumerate(test_df_aa_ss8.itertuples()):
        annots = {}
        prot_id = row.proteins
        # BlastKNN
        if prot_id in diamond_scores:
            sim_prots = diamond_scores[prot_id]  # FS: diamond score
            allgos = set()
            total_score = 0.0
            for p_id, score in sim_prots.items():
                allgos |= annotations[prot_index[p_id]]
                total_score += score           # FS: diamond score
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


    ############ 这里和 FindAlpha 是不一样的 ######################
    # 这里和 Alpha_FindAlpha.py 中的 eval_alpha 很像
    # DeepGOPlus
    deep_preds = []  # 都用在哪儿了？
    go_set = go_rels.get_namespace_terms(NAMESPACES[ont])
    go_set.remove(FUNC_DICT[ont])  # 把最根的那三个之一给删了？

    labels = test_df_aa_ss8['prop_annotations'].values
    labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))
    # print(len(go_set))

    # ####  FS 添加 # block X # 这里添加一下，把labels的每一个元素{GO,GO,GO},{GO,GO},{GO}和terms的2088做交集
    if run_label_block == 'T':
        # print('it is running block X !!!!!!!!!!!!!!!!!')
        terms_set = set(terms_aa)
        labels = [s.intersection(terms_set) for s in labels]
        # print(labels)
        # print('len(labels & terms_set) = ', len(labels))


    # alphas = {NAMESPACES['mf']: 0.55, NAMESPACES['bp']: 0.59, NAMESPACES['cc']: 0.46}
    alphas = {NAMESPACES['bp']: 0, NAMESPACES['cc']: 0, NAMESPACES['mf']: 0}  # 重新初始化？
    betas = {NAMESPACES['bp']: 0, NAMESPACES['cc']: 0, NAMESPACES['mf']: 0}  # 初始化，original=0


    # 这里又和 FindAlpha 一样了  （这里比较卡--慢）

    ##################################################
    ############## 引入 json 参数 #####################
    ##################################################

    # 第一步，建库 deep_preds
    # 如果alpha=NA，则采用json数据，否则alpha=数字，或外来click引入
    if alpha == 'json' and beta == 'json':
        # 从last_release_metadata文件中获取alpha ###
        print('alpha/beta is from json, alpha/beta = ', alpha, beta)
        last_release_metadata = 'Alpha_last_release.json'
        with open(last_release_metadata, 'r') as f:
            print('Reading file from json')
            last_release_data = json.load(f)

            alpha = last_release_data["alphas"][ont]  # 从 json中读取数据
            print('111111111', type(alpha))
            alphas[NAMESPACES[ont]] = alpha   # ????????? 方便下面的迭代中使用alpha

            beta = last_release_data["betas"][ont]  # 从 json中读取数据
            print('22222222', type(beta))
            betas[NAMESPACES[ont]] = beta  # ????????? 方便下面的迭代中使用alpha




    else:  # alpha = int，也就是在click中又指定
        print('alpha/beta is from click, alpha/beta = ', alpha, beta)
        print('type_alpha/beta = ', type(alpha), type(beta))
        alphas[NAMESPACES[ont]] = alpha
        betas[NAMESPACES[ont]] = beta





    ### FS 从find_alpha的json文件中获取 alpha
    print('lol')
    print('alphas, beta = ', alphas, betas)  # eg. alphas {'molecular_function': 0.5, 'biological_process': 0, 'cellular_component': 0}

    # 这是核心 !!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i, row in enumerate(test_df_aa_ss8.itertuples()):

        annots_dict = blast_preds[i].copy()
        for go_id in annots_dict:
            # annots_dict[go_id] *= 1 - alphas[go_rels.get_namespace(go_id)]  # blast/diamond * (1-alpha)    original...
            annots_dict[go_id] *= 1 - alphas[go_rels.get_namespace(go_id)] - betas[go_rels.get_namespace(go_id)]  # blast/diamond * (1-alpha)

        for j, score_aa in enumerate(row.preds_aa):  # test_df_aa_ss8 才有 preds_aa
            go_id = terms_aa[j]
            score_aa *= alphas[go_rels.get_namespace(go_id)]  # 深度学习deepgo_preds * (alpha)

            if go_id in annots_dict:
                annots_dict[go_id] += score_aa
            else:
                annots_dict[go_id] = score_aa
                
                
        for k, score_ss8 in enumerate(row.preds_ss8):  # deep深度学习中算的预测值preds
            go_id = terms_ss8[k]   # FS: 这一点修改的很 有趣，简洁。因为ss8是按它自己的terms排序的!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            score_ss8 *= betas[go_rels.get_namespace(go_id)] # 深度学习preds_ss8 * beta
            if go_id in annots_dict:  # if 深度学习的go_id在diamond对比数据中
                annots_dict[go_id] += score_ss8
            else:
                annots_dict[go_id] = score_ss8


                
        deep_preds.append(annots_dict)




    ### 第二步：###
    # 上面是统计估分，下面开始 评估evaluate
    fmax = 0.0
    tmax = 0.0
    precisions = []
    recalls = []
    smin = 1000000.0  # 尽可能大
    rus = []
    mis = []
    for t in range(1, 101):  # the range in this loop has influence in the AUPR output
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

    # 下面和 Alpha_FindAlpha.py就不一样了，要保存结果了
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

    ### FS 保存Smin,Fmax,AUPR ###
    with open(f'results_alphabeta/Fmax_AUPR_Smin_{ont}_alpha{alpha:0.2f}_beta{beta:0.2f}.txt', 'w') as f:
        f.write(f'Length of test set: {len(test_df_aa_ss8)} \n')
        f.write(f'threshold: {tmax} \n')
        f.write(f'alpha={alpha:0.2f} \n')
        f.write(f'Smin: {smin:0.3f} \n')
        f.write(f'Fmax: {fmax:0.3f} \n')
        f.write(f'AUPR: {aupr:0.3f} \n')
    ### FS 保存Smin,Fmax,AUPR  end ####



    # plt.figure()
    # lw = 2
    # plt.plot(recalls, precisions, color='darkorange',
    #          lw=lw, label=f'AUPR curve (area = {aupr:0.2f})')
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    # plt.xlabel('Recall')
    # plt.ylabel('Precision')
    # plt.title('Area Under the Precision-Recall curve')
    # plt.legend(loc="lower right")
    #### plt.savefig(f'results_alphabeta/aupr_{ont}_{alpha:0.2f}.pdf')  # 不保存pdf了，就是个precission-recall图
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