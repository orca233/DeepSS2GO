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
# from step0_TrainTestSetting_global import *  # 这一行有问题，一旦取消注释 找不到data/terms_....pkl

# from step0_TrainTestSetting_local_aa import *
# from step0_TrainTestSetting_global import path_base


logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


@ck.command()
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
# @ck.option('--go-file', '-gf', default=path_base + 'pub_data/go.obo', help='Gene Ontology file in OBO Format')  # original: 'data/go.obo'  此行原本没有，FS加
@ck.option('--train-data-file', '-trdf', default='data/train_data_aa.pkl', help='Data file with training features')
# @ck.option('--predictions-file', '-pf', default='data/predictions.pkl', help='XX')
@ck.option('--test-data-file-aa-ss8', '-tsdf', default='data/predictions_aa_ss8.pkl', help='Test data file')  # test_data_file_aa_ss8 只有一个 prop_annotations, 两个 preds_aa & preds_ss8

@ck.option('--terms-file-aa', '-tf', default='./data/terms_gominre_trxte_aa.pkl', help='set of annotations')  # original 'data/terms_all.pkl'
@ck.option('--terms-file-ss8', '-tf', default='./data/terms_gominre_trxte_ss8.pkl', help='set of annotations')  # original 'data/terms_all.pkl'

@ck.option('--diamond-scores-file', '-dsf', default='data/diamond_aa.res', help='Diamond output')  # test_diamond.res / blast_output_diamond.res
@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')

# new
@ck.option('--run-label-block', '-rbx', default='T', help='judge whether run block X')
@ck.option('--alpha-range', '-ar', default='0, 101, 10', type=str, help='xx')  # [25, 80, 1]   params_local['alpha_range']， 初筛：'0, 101, 10'，细筛：'30, 50, 2'
@ck.option('--beta-range', '-br', default='0, 101, 10', type=str, help='xx')  # [25, 80, 1]   params_local['beta_range']，初筛：'0, 101, 10'，细筛：'10, 30, 2'
@ck.option('--num-cpu-cores', '-ncc', default=80, type=int, help='xx')  # 60
@ck.option('--last-release-metadata', '-lrm', default='Alpha_last_release.json', help='xx')


def main(go_file, train_data_file, test_data_file_aa_ss8, terms_file_aa, terms_file_ss8, diamond_scores_file, ont,
         run_label_block, alpha_range, beta_range, num_cpu_cores, last_release_metadata):  # FS加的go_file

    # 预处理
    # test_data_aa_file = 'data/predictions_aa_ss8.pkl'  # 用的是综合的 aa_ss8.pkl
    # test_data_ss8_file = 'data/predictions_ss8.pkl'

    go_rels = Ontology(go_file, with_rels=True)
    terms_df_aa = pd.read_pickle(terms_file_aa)
    terms_aa = terms_df_aa['terms'].values.flatten()  # ['GO:0051254' 'GO:0045637' 'GO:1905368' ...  'GO:0030665']
    terms_dict = {v: i for i, v in enumerate(terms_aa)}  # {'GO:0051254': 0, 'GO:0045637': 1, 'GO:1905368': 2, 'GO:0042129': 3, ... 'GO:0030665': 2087}

    # FS:
    terms_df_ss8 = pd.read_pickle(terms_file_ss8)  # original: data/terms_all_ss8.pkl
    terms_ss8 = terms_df_ss8['terms'].values.flatten()

    train_df = pd.read_pickle(train_data_file)
    test_df_aa_ss8 = pd.read_pickle(test_data_file_aa_ss8)
    # test_df_ss8 = pd.pdread_pickle(test_data_ss8_file)
    print("Length of test set: " + str(len(test_df_aa_ss8)))

    annotations = train_df['prop_annotations'].values  # propagate annotation?????
    annotations = list(map(lambda x: set(x), annotations))

    test_annotations_aa_ss8 = test_df_aa_ss8['prop_annotations'].values
    test_annotations_aa_ss8 = list(map(lambda x: set(x), test_annotations_aa_ss8))

    # test_annotations_ss8 = test_df_ss8['prop_annotations'].values
    # test_annotations_ss8 = list(map(lambda x: set(x), test_annotations_ss8))  # 好像没啥用

    go_rels.calculate_ic(annotations + test_annotations_aa_ss8)


    #####################################################################
    ####################### Diamond 计算 IC ##############################
    #####################################################################
    # Print IC values of terms  # information content
    ics = {}
    for term in terms_aa:
        ics[term] = go_rels.get_ic(term)

    prot_index = {}
    for i, row in enumerate(train_df.itertuples()):
        '''
        row:
        Pandas(Index=60795, index=437551, 
        proteins='S6A11_HUMAN', 
        accessions='P48066; B2R6U6; Q8IYC9;', 
        sequences='MTAEKALPLGNGKAAEEARESEAPGGGCSSGGAAPARHPRVKRDKAVHERGHWNNKVEFVLSVAGEIIGLGNV。。。', 
        annotations=['GO:0005737|IEA', 'GO:0098982|IEA', 'GO:0016020|TAS', 'GO:0043005|IBA', 'GO:0005886|IBA', 。。。], 
        interpros=['IPR000175', 'IPR002982', 'IPR037272'], orgs='9606', 
        exp_annotations=['GO:0016020', 'GO:0005332', 'GO:0008028', 'GO:0051936', 'GO:0015718'], 
        prop_annotations=['GO:1905039', 'GO:0098662', 'GO:0015698', 'GO:0034220', 'GO:0098661',。。。], cafa_target=True)
        '''
        prot_index[row.proteins] = i
    #prot_index  # len = 12576， {'S6A11_HUMAN': 0, 'MBOA5_HUMAN': 1, 'SMPX_HUMAN': 2, 'IFIX_HUMAN': 3,。。。}


    ####################### calc BLAST diamond #########################################
    # BLAST Similarity (Diamond)  这里用aa，不用ss8
    diamond_scores = {}
    with open(diamond_scores_file) as f:
        for line in f:
            it = line.strip().split()
            if it[0] not in diamond_scores:
                diamond_scores[it[0]] = {}
            diamond_scores[it[0]][it[1]] = float(it[2])
    # diamond_scores  # len=503,
    # {'PI42B_HUMAN': {'PI42A_HUMAN': 622.0, 'PI42C_HUMAN': 520.0, 'PI51B_HUMAN': 170.0,。。}
    #  'MACD2_HUMAN': {'MACD1_HUMAN': 238.0, 'GDAP2_HUMAN': 115.0, 。。。} 。。。 {} 。。。}


    blast_preds = []
    # print('Diamond preds')
    for i, row in enumerate(test_df_aa_ss8.itertuples()):
        annots = {}
        prot_id = row.proteins
        # BlastKNN
        if prot_id in diamond_scores:  # 在 diamond_score中，逐一过test 503个
            sim_prots = diamond_scores[prot_id]  # len=6,
            # {'PI42A_HUMAN': 622.0, 'PI42C_HUMAN': 520.0, 'PI51B_HUMAN': 170.0, 'PI51A_HUMAN': 162.0,
            # 'PI51C_HUMAN': 158.0, 'PI5L1_HUMAN': 94.0}   加和=1726
            allgos = set()
            total_score = 0.0
            for p_id, score in sim_prots.items():  # 在 diamond_score每个test所匹配的不定数目的train的结果，逐一过
                allgos |= annotations[prot_index[p_id]]  # 把train.pkl (prot_index)对应的annotations GO term并集到allgos
                total_score += score  # 这是diamond中每个文库(train)的匹配score之和  第一组的6个加和=1726.0
            # print('666666666        ', total_score)
            allgos = list(sorted(allgos))  # ['GO:0000045', 'GO:0000285', 'GO:0001505', 'GO:0001726'...]

            sim = np.zeros(len(allgos), dtype=np.float32)
            # sim.shape = (404,), 对于第一个test蛋白一共有6个train的蛋白质对应，这6个又对应 404 个GO prop_annotation

            # 下一个for循环，下一个test蛋白(1506,)...
            # print('------------------' ,type(sim))
            # print(sim.shape)

            for j, go_id in enumerate(allgos):  # 逐一计算每个GO term对应train文库给出的评分
                s = 0.0
                for p_id, score in sim_prots.items():
                    if go_id in annotations[prot_index[p_id]]:
                        s += score
                sim[j] = s / total_score

            # print('666666666        ', sim)
            ind = np.argsort(-sim)
            for go_id, score in zip(allgos, sim):  # 没看懂
                annots[go_id] = score
        blast_preds.append(annots)
    # blast_preds, len = 662, [{'GO:0000045': 0.6616454, 'GO:0000285': 0.19235227, 'GO:0001505': 0.091541134,


    ##################################################
    ########## 引入 json 参数， find alpha & beta ######
    ##################################################

    # last_release_metadata = 'Alpha_last_release.json'

    with open(last_release_metadata, 'r') as f:
        last_release_data = json.load(f)
    print('111', last_release_data)

    #### find alpha 方程里调用了 eval_alphas
    last_release_data['alphas'][ont], last_release_data['betas'][ont] = find_alpha_beta(ont, test_df_aa_ss8, blast_preds, go_rels, terms_aa, terms_ss8, run_label_block, alpha_range, beta_range, num_cpu_cores)
    print('2222', last_release_data)
    with open(last_release_metadata, 'w') as f:
        json.dump(last_release_data, f)


    # last_release_data['alphas'][ont] = find_alpha_beta(ont, test_df_aa_ss8, blast_preds, go_rels, terms)
    # print('2222', last_release_data)
    # with open(last_release_metadata, 'w') as f:
    #     json.dump(last_release_data, f)


# ################ FS: ALGORITHM TO FIND BEST ALPHAS & BETAS PARAMETER ##################################

def eval_alphas_betas(alpha, beta, ont, test_df_aa_ss8, blast_preds, go_rels, terms_aa, terms_ss8, run_label_block):
    print('alpha, beta = ', alpha, beta)

    deep_preds = []
    go_set = go_rels.get_namespace_terms(NAMESPACES[ont])
    go_set.remove(FUNC_DICT[ont])

    labels = test_df_aa_ss8['prop_annotations'].values
    labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))



    # ####  FS 添加 # block X # 这里添加一下，把labels的每一个元素{GO,GO,GO},{GO,GO},{GO}和terms的2088做交集    这个block仅限于此！！！！！！！！！！！！！！！！！！！！！！！
    if run_label_block == 'T':
        # print('it is running block X !!!!!!!!!!!!!!!!!')
        terms_set = set(terms_aa)
        labels = [s.intersection(terms_set) for s in labels]
        # print(labels)
        # print('len(labels & terms_set) = ', len(labels))


    alphas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}  # 初始化，original=0
    alphas[NAMESPACES[ont]] = alpha  # 比如ont=bp，则alphas[biological_process]=alpha外部引入
    # print('77777777777')
    # print(NAMESPACES[ont])
    # print(alphas)
    # print(alpha)

    betas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}  # 初始化，original=0
    betas[NAMESPACES[ont]] = beta  # 比如ont=bp，则alphas[biological_process]=alpha外部引入



    '''
    row: 这是对于某一个i，也就是test中的一个蛋白质的row：
    Pandas(Index=49935, index=314240, proteins='PI42B_HUMAN', 
    accessions='P78356; Q5U0E8; Q8TBP2;', 
    sequences='MSSNCTSTTAVAVAPLSASKTKTKKKHFVCQKVKLFRASEPILSVLMWGVNHTINELSNVPVPVMLMPDDFKAYSKIKVDNHLFNK...', 
    annotations=['GO:0005776|IMP', 'GO:0005829|TAS', 'GO:0005789|IEA', 'GO:0005654|IDA', ...], 
    interpros=['IPR023610', 'IPR002498', 'IPR027484'], 
    orgs='9606', 
    exp_annotations=['GO:0005776', 'GO:0005829', 'GO:0005654', 'GO:0005634', 'GO:0005886', 'GO:0016308'...], 
    cafa_target=True, labels=array([0., 0., 0., ..., 0., 0., 0.], 
    dtype=float32), 
    !!! preds=array([0.0264879 , 0.00491404, 0.01317153, ..., 0.01190471, 0.01042569, 0.00605111], !!!
    dtype=float32))

    len(row.preds) = 2088
    '''

    # 这是核心 !!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i, row in enumerate(test_df_aa_ss8.itertuples()):

        annots_dict = blast_preds[i].copy()  # test中的第 i 个元素，在blast算的值
        for go_id in annots_dict:
            annots_dict[go_id] *= 1 - alphas[go_rels.get_namespace(go_id)] - betas[go_rels.get_namespace(go_id)]  # blast * (1-alpha- beta)

        for j, score_aa in enumerate(row.preds_aa):  # deep深度学习中算的预测值preds  # test_df_aa_ss8 才有 preds_aa
            go_id = terms_aa[j]
            #original: score_aa *= 1 - alphas[go_rels.get_namespace(go_id)]  # 深度学习preds * (1-alpha)
            score_aa *= alphas[go_rels.get_namespace(go_id)]  # 深度学习preds_aa * alpha

            if go_id in annots_dict:  # if 深度学习的go_id在diamond对比数据中
                annots_dict[go_id] += score_aa
            else:
                annots_dict[go_id] = score_aa

        # 明天接着改，row需要改row的范围，然后把上面j换成beta，这里换成 1-alpha - beta。下面alpha取值范围间隔可设为0.03!!!!
        for k, score_ss8 in enumerate(row.preds_ss8):  # deep深度学习中算的预测值preds
            go_id = terms_ss8[k]   # FS: 这一点修改的很 有趣，简洁。因为ss8是按它自己的terms排序的!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            score_ss8 *= betas[go_rels.get_namespace(go_id)] # 深度学习preds_ss8 * beta
            if go_id in annots_dict:  # if 深度学习的go_id在diamond对比数据中
                annots_dict[go_id] += score_ss8
            else:
                annots_dict[go_id] = score_ss8


        # alpha=0.2:  'GO:0000045' = 0.6616454 * 0.01 + score = 0.042534434720873836。也可能*0.00修改alphas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}
        # annots_dict = {'GO:0000045': 0.042534434720873836, 'GO:0000285': 0.03847045302391053, 'GO:0001505': 0.005691665941849351,...
        deep_preds.append(annots_dict)
        # deep_preds = [{'GO:0000045': 0.042534434720873836, 'GO:0000285': 0.03847045302391053, ...},{}。。。]


    # 上面是统计估分，下面开始评估evaluate
    fmax = 0.0
    tmax = 0.0
    precisions = []
    recalls = []
    smin = 1000000.0
    rus = []
    mis = []
    for t in range(1, 101):  # t=threshold the range in this loop has influence in the AUPR output
        threshold = t / 100.0
        preds = []
        for i, row in enumerate(test_df_aa_ss8.itertuples()):  # 对于test中逐一过每一个蛋白
            annots = set()  # 超过阈值threshold的GO合集
            for go_id, score in deep_preds[i].items():
                if score >= threshold:
                    annots.add(go_id)

            new_annots = set()
            for go_id in annots:
                new_annots |= go_rels.get_ancestors(go_id)  # 把爸爸们找出来
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

    precisions = np.array(precisions)
    recalls = np.array(recalls)
    sorted_index = np.argsort(recalls)
    recalls = recalls[sorted_index]
    precisions = precisions[sorted_index]
    aupr = np.trapz(precisions, recalls)

    # Fmax,AUPR 的 mf & bp ~0.5, cc ~0.6.
    # smin 的：mf~16, cc~13, bp ~50
    # 因为输出为三者加权 np.sum([smin, -fmax, -aupr])，越小越好
    if ont == 'mf':
        smin /= 15
    elif ont == 'bp':
        smin /= 40
    elif ont == 'cc':
        smin /= 10

    return alpha, beta, -fmax  # original: alpha, np.sum([smin, -fmax, -aupr])




def find_alpha_beta(ont, test_df_aa_ss8, blast_preds, go_rels, terms_aa, terms_ss8, run_label_block, alpha_range, beta_range, num_cpu_cores):
    '''
    逐一过所有的alpha （范围：alpha_range）
    并将每一个alpha传到 eval_alphas 函数中，得到最大的 Fmax
    每一个 (alpha, Fmax) 组成 results，最后选出最大的那个Fmax对应的alpha
    '''
    extra = [ont, test_df_aa_ss8, blast_preds, go_rels, terms_aa, terms_ss8, run_label_block]

    '''
    # # original:
    # 初筛
    # inputs_alpha = range(0, 101, 10)  ## ?????   # original range(45, 75, 1)
    # inputs_beta = range(0, 101, 10)
    # 细筛
    # inputs_alpha = range(20, 51, 2)
    # inputs_beta = range(10, 31, 2)
    '''

    # FS change
    # list2range: '25, 80, 1' to range(25, 80, 1)  # 把 alpha_range 这个list转成eange
    str_values_alpha = alpha_range.split(',')
    int_values_alpha = [int(value.strip()) for value in str_values_alpha]
    inputs_alpha = range(*int_values_alpha)  # range(25, 80, 1)

    str_values_beta = beta_range.split(',')
    int_values_beta = [int(value.strip()) for value in str_values_beta]
    inputs_beta = range(*int_values_beta)  # range(25, 80, 1)

    print('inputs_alpha/beta = ', inputs_alpha, inputs_beta)  # FS
    # FS change end

    num_cores = num_cpu_cores  # most cpu can be used, e.g. 60

    # original: results = Parallel(n_jobs=num_cores)(delayed(eval_alphas_betas())(i / 100, *extra) for i in inputs)
    results = Parallel(n_jobs=num_cores)(
        delayed(eval_alphas_betas)(i / 100, j / 100, *extra)
        for i in inputs_alpha
        for j in inputs_beta
        if i + j <= 100  # 这个很重要哈哈哈哈哈，其实是左下角的三角形

        # ???? beta=1 fmax=0.28???
        # bug: inputs_alpha = range(0, 1, 5)
        #     inputs_beta = range(0, 20, 1)

    )

    # 该函数以i / 100作为第一个参数，后跟extra列表中的所有元素
    # i的取值从20到79，步长为1，通过并行化处理，将这些函数调用的结果存储在results列表中

    print('22222222222222222222222222 results:')
    print(results)
    # results = [(0.2, -0.2244738075593926), (0.21, -0.2323958248068838), (0.22, -0.2419488147142912), , (0.23, -0.2503843465951797)...

    chosen = min(results, key=lambda x: x[2])  # original: x[1]    x = [alpha, beta, score(-fmax)] 从eval_alphas_betas中return出来的
    # 在results列表中选择具有最小值的元素。这里使用了lambda表达式作为key函数，将列表中的每个元素x的第二个元素作为比较的依据。
    print('chosen')
    print(chosen)
    # chosen = (0.54, -0.33197674380787534)

    return chosen[0], chosen[1]  # chosen[0]=alpha, chosen[1]=beta
    # 返回chosen列表的第一个元素作为函数的结果。就是alpha。chosen第二个值是最小的 smin-fmax-aupr








# ################ original: ALGORITHM TO FIND BEST ALPHAS PARAMETER ##################################
#
# def eval_alphas(alpha, ont, test_df_aa_ss8, blast_preds, go_rels, terms):
#     deep_preds = []
#     go_set = go_rels.get_namespace_terms(NAMESPACES[ont])
#     go_set.remove(FUNC_DICT[ont])
#
#     labels = test_df_aa_ss8['prop_annotations'].values
#     labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))
#
#     alphas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}  # 初始化，original=0
#     alphas[NAMESPACES[ont]] = alpha  # 比如ont=bp，则alphas[biological_process]=alpha外部引入
#     print('77777777777')
#     print(NAMESPACES[ont])
#     print(alphas)
#     print(alpha)
#
#     for i, row in enumerate(test_df_aa_ss8.itertuples()):
#         annots_dict = blast_preds[i].copy()  # test中的第 i 个元素，在blast算的值
#         for go_id in annots_dict:
#             annots_dict[go_id] *= alphas[go_rels.get_namespace(go_id)]  # blast * alpha
#
#         for j, score in enumerate(row.preds):  # deep深度学习中算的预测值preds
#             go_id = terms[j]
#             score *= 1 - alphas[go_rels.get_namespace(go_id)]  # 深度学习preds * (1-alpha)
#             if go_id in annots_dict:  # if 深度学习的go_id在diamond对比数据中
#                 annots_dict[go_id] += score
#             else:
#                 annots_dict[go_id] = score
#
#
#         # 明天接着改，row需要改row的范围，然后把上面j换成beta，这里换成 1-alpha - beta。下面alpha取值范围间隔可设为0.03!!!!
#         for k, score in enumerate(row.preds):  # deep深度学习中算的预测值preds
#             go_id = terms[k]
#             score *= 1 - alphas[go_rels.get_namespace(go_id)]  # 深度学习preds * (1-alpha)
#             if go_id in annots_dict:  # if 深度学习的go_id在diamond对比数据中
#                 annots_dict[go_id] += score
#             else:
#                 annots_dict[go_id] = score
#
#
#
#
#
#         # alpha=0.2:  'GO:0000045' = 0.6616454 * 0.01 + score = 0.042534434720873836。也可能*0.00修改alphas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}
#         # annots_dict = {'GO:0000045': 0.042534434720873836, 'GO:0000285': 0.03847045302391053, 'GO:0001505': 0.005691665941849351,...
#         deep_preds.append(annots_dict)
#         # deep_preds = [{'GO:0000045': 0.042534434720873836, 'GO:0000285': 0.03847045302391053, ...},{}。。。]
#
#
#     # 上面是统计估分，下面开始评估evaluate
#     fmax = 0.0
#     tmax = 0.0
#     precisions = []
#     recalls = []
#     smin = 1000000.0
#     rus = []
#     mis = []
#     for t in range(1, 101):  # t=threshold the range in this loop has influence in the AUPR output
#         threshold = t / 100.0
#         preds = []
#         for i, row in enumerate(test_df_aa_ss8.itertuples()):  # 对于test中逐一过每一个蛋白
#             annots = set()  # 超过阈值threshold的GO合集
#             for go_id, score in deep_preds[i].items():
#                 if score >= threshold:
#                     annots.add(go_id)
#
#             new_annots = set()
#             for go_id in annots:
#                 new_annots |= go_rels.get_ancestors(go_id)  # 把爸爸们找出来
#             preds.append(new_annots)
#
#         # Filter classes
#         preds = list(map(lambda x: set(filter(lambda y: y in go_set, x)), preds))
#
#         fscore, prec, rec, s, ru, mi, fps, fns = evaluate_annotations(go_rels, labels, preds)
#         avg_fp = sum(map(lambda x: len(x), fps)) / len(fps)
#         avg_ic = sum(map(lambda x: sum(map(lambda go_id: go_rels.get_ic(go_id), x)), fps)) / len(fps)
#         # print(f'{avg_fp} {avg_ic}')
#         precisions.append(prec)
#         recalls.append(rec)
#         # print(f'Fscore: {fscore}, Precision: {prec}, Recall: {rec} S: {s}, RU: {ru}, MI: {mi} threshold: {threshold}')
#         if fmax < fscore:
#             fmax = fscore
#             tmax = threshold
#         if smin > s:
#             smin = s
#
#     precisions = np.array(precisions)
#     recalls = np.array(recalls)
#     sorted_index = np.argsort(recalls)
#     recalls = recalls[sorted_index]
#     precisions = precisions[sorted_index]
#     aupr = np.trapz(precisions, recalls)
#
#     # Fmax,AUPR 的 mf & bp ~0.5, cc ~0.6.
#     # smin 的：mf~16, cc~13, bp ~50
#     # 因为输出为三者加权 np.sum([smin, -fmax, -aupr])，越小越好
#     if ont == 'mf':
#         smin /= 15
#     elif ont == 'bp':
#         smin /= 40
#     elif ont == 'cc':
#         smin /= 10
#
#     return alpha, -fmax  # original: alpha, np.sum([smin, -fmax, -aupr])



###############################################################################

# def find_alpha(ont, test_df_aa_ss8, blast_preds, go_rels, terms):
#     extra = [ont, test_df_aa_ss8, blast_preds, go_rels, terms]
#     inputs = range(45, 75, 1)  ## ????? alpha 只能从这个范围找么？？？？？？？？？  # original range(45, 75, 1)
#     print('inputs=', inputs)  # FS
#
#     num_cores = 30  # most cpu can be used
#
#     results = Parallel(n_jobs=num_cores)(delayed(eval_alphas)(i / 100, *extra) for i in inputs)
#     # 该函数以i / 100作为第一个参数，后跟extra列表中的所有元素
#     # i的取值从20到79，步长为1，通过并行化处理，将这些函数调用的结果存储在results列表中
#
#     print('22222222222222222222222222')
#     print(results)
#     # results = [(0.2, -0.2244738075593926), (0.21, -0.2323958248068838), (0.22, -0.2419488147142912), , (0.23, -0.2503843465951797)...
#
#     chosen = min(results, key=lambda x: x[1])
#     # 在results列表中选择具有最小值的元素。这里使用了lambda表达式作为key函数，将列表中的每个元素x的第二个元素作为比较的依据。
#     print(chosen)
#     # chosen = (0.54, -0.33197674380787534)
#
#     return chosen[0]  # 返回chosen列表的第一个元素作为函数的结果。












#####################################################################################

if __name__ == '__main__':
    alphas = {'mf': 0, 'bp': 0, 'cc': 0}
    main()
