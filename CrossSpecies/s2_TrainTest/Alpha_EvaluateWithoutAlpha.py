'''
from DeepGOPlus: evaluate_deepgoplus.py

相比于Alpha_FindAlpha.py，（从25-80逐一过一遍alpha，每一个alpha又过 t(threshold)0-100）
这个py文件只针对上述筛选出的某一个alpha进行，2min

'''


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

# logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


@ck.command()
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='Data file with training features')
@ck.option('--test-data-file', '-tsdf', default='data/predictions.pkl', help='Test data file')
# @ck.option('--predictions-file', '-pf', default='data/predictions.pkl', help='XX')

@ck.option('--terms-file', '-tf', default='data/terms_gominre_trxte.pkl', help='Data file with sequences and complete set of annotations')  # original: data/terms_all.pkl
# @ck.option('--diamond-scores-file', '-dsf', default='data/diamond_aa.res', help='Diamond output')
@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--alpha', '-a', default='json', help='alpha = json(with quote) or 0-1(without quote, eg 0.3 float)')  # 如果alpha='json'，则采用json数据，否则alpha=数字，或外来click引入
@ck.option('--go-file', '-gf', default=params_local['path_base'] + 'pub_data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
# @ck.option('--go-file', '-gf', default='../../pub_data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
@ck.option('--run-label-block', '-rbx', default='T', help='judge whether run block X')



def main(train_data_file, test_data_file, terms_file, ont, alpha, go_file, run_label_block):  # FS 添加go_file  diamond_scores_file
    # # 从last_release_metadata文件中获取alpha ###
    # last_release_metadata = 'Alpha_last_release.json'
    # with open(last_release_metadata, 'r') as f:
    #     last_release_data = json.load(f)
    #     alpha = last_release_data['alphas'][ont]  # # original 这里更新了alpha所以用不到ck中的引入信息
    #
    # print('alpha initial=', alpha)

    go_rels = Ontology(go_file, with_rels=True)
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()
    # terms_dict = {v: i for i, v in enumerate(terms)}  # 没用上？

    train_df = pd.read_pickle(train_data_file)
    test_df = pd.read_pickle(test_data_file)  # predictions.pkl的信息, original: test_data_file
    print("Length of test set: " + str(len(test_df)))

    annotations = train_df['prop_annotations'].values
    annotations = list(map(lambda x: set(x), annotations))
    test_annotations = test_df['prop_annotations'].values
    test_annotations = list(map(lambda x: set(x), test_annotations))
    go_rels.calculate_ic(annotations + test_annotations)
    # print(go_rels)





    ''' withoutAlpha 注释掉
    #####################################################################
    ########################## Diamond 计算IC ############################
    #####################################################################
    # Print IC values of terms
    ics = {}
    for term in terms:
        ics[term] = go_rels.get_ic(term)

    # print('3333333')
    # print(ics)
    # print(len(ics))


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
    for i, row in enumerate(test_df.itertuples()):
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
    # blast_preds, len = 662, [{'GO:0000045': 0.6616454, 'GO:0000285': 0.19235227, 'GO:0001505': 0.091541134,
    '''


    ############ 这里和 FindAlpha 是不一样的 ######################
    # 这里和 Alpha_FindAlpha.py 中的 eval_alpha 很像
    # DeepGOPlus
    deep_preds = []  # 都用在哪儿了？
    go_set = go_rels.get_namespace_terms(NAMESPACES[ont])  # 这里只留下对应ont的go（不确定哈，猜的）!!!!!!!!!!!!!
    # print('11111111111')
    # print(len(go_set))
    # go_set = {'GO:0061685', 'GO:0034909', 'GO:0016167', 'GO:0004902',  len = 12407

    go_set.remove(FUNC_DICT[ont])  # 把最根的那三个之一给删了？

    labels = test_df['prop_annotations'].values
    labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))
    # print(len(go_set))

    # ####  FS 添加 # block X # 这里添加一下，把labels的每一个元素{GO,GO,GO},{GO,GO},{GO}和terms的2088做交集
    if run_label_block == 'T':
        # print('it is running block X !!!!!!!!!!!!!!!!!')
        terms_set = set(terms)
        labels = [s.intersection(terms_set) for s in labels]
        # print(labels)
        # print('len(labels & terms_set) = ', len(labels))




    # alphas = {NAMESPACES['mf']: 0.55, NAMESPACES['bp']: 0.59, NAMESPACES['cc']: 0.46}
    alphas = {NAMESPACES['bp']: 0, NAMESPACES['cc']: 0, NAMESPACES['mf']: 0}  # 重新初始化？



    # 这里又和 FindAlpha 一样了  （这里比较卡--慢）

    ##################################################
    ############## 引入 json 参数 #####################
    ##################################################


    # 第一步，建库 deep_preds
    # 如果alpha=NA，则采用json数据，否则alpha=数字，或外来click引入
    if alpha == 'json':
        # 从last_release_metadata文件中获取alpha ###
        print('alpha is from json, alpha = ', alpha, 'alpha_type = ', type(alpha))

        last_release_metadata = 'Alpha_last_release.json'
        with open(last_release_metadata, 'r') as f:
            print('Reading file from json')
            last_release_data = json.load(f)
            alpha = last_release_data["alphas"][ont]  # 从 json中读取数据
            print('type_alpha = ', type(alpha))
            alphas[NAMESPACES[ont]] = alpha   # ????????? 方便下面的迭代中使用alpha
    else:  # e.g. alpha = 0.6，也就是在click中又指定。因为这个0.6是str，所以要转格式
        print('alpha is from click, alpha = ', alpha)
        print('type_alpha = ', type(alpha))
        alpha = float(alpha)
        alphas[NAMESPACES[ont]] = alpha
        print('updated_type_alpha = ', type(alphas[NAMESPACES[ont]]))

    ### FS 从find_alpha的json文件中获取 alpha
    print('lol')
    print('alphas=', alphas)  # eg. alphas {'molecular_function': 0.5, 'biological_process': 0, 'cellular_component': 0}


    ###########################################################
    # 这是核心 !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ###########################################################
    # (1 - alpha - beta) * diamond + alpha * preds_aa + beta * preds_ss8
    # alpha=0: 全由 diamond 统计
    # alpha=1: 全由 deepSS2GO_aa 统计
    # beta=1: 全由 deepSS2GO_ss8 统计

    blast_preds = [{} for _ in range(10000)]  # blast_preds = [{},{},{},{},{},...10000个{}]  # 专为Evaluate Without Alpha准备！！！

    for i, row in enumerate(test_df.itertuples()):
        annots_dict = blast_preds[i].copy()
        for go_id in annots_dict:
            annots_dict[go_id] *= 1 - alphas[go_rels.get_namespace(go_id)]  # blast/diamond * (1-alpha)
            # print('go.get_namespace(go_id) = ', go_rels.get_namespace(go_id))  # 其实这里还是混着的,3种GO都有
            # print('alphas[go.get_namespace(go_id)] = ', alphas[go_rels.get_namespace(go_id)])

        for j, score in enumerate(row.preds):
            go_id = terms[j]
            score *= alphas[go_rels.get_namespace(go_id)]  # 深度学习deepgo_preds * (alpha)
            if go_id in annots_dict:
                annots_dict[go_id] += score
            else:
                annots_dict[go_id] = score
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
        # print('threshold =', t)
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
            # print('new_annots = ', new_annots)
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
    with open(f'results_alpha/Fmax_AUPR_Smin_{ont}_{alpha:0.2f}.txt', 'w') as f:
        f.write(f'Length of test set: {len(test_df)} \n')
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
    #### plt.savefig(f'results_alpha/aupr_{ont}_{alpha:0.2f}.pdf')  # 不保存pdf了，就是个precission-recall图
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