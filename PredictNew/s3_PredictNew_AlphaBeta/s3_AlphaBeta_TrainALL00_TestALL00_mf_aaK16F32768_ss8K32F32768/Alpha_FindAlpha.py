'''
from DeepGOPlus, find_alphas.py
'''

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
from step0_TrainTestSetting_local_aa import params_local as params_local_aa
from step0_TrainTestSetting_local_ss8 import params_local as params_local_ss8
from step0_TrainTestSetting_global import path_base

# logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

# os.system('mkdir -p results')
# os.system('mkdir -p results_alpha')
# os.system('mkdir -p results_alphabeta')

@ck.command()
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
# @ck.option('--go-file', '-gf', default=path_base + 'pub_data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='Data file with training features')
@ck.option('--predictions-file', '-pf', default='data/predictions.pkl', help='XX')
# @ck.option('--test-data-file', '-tsdf', default='data/predictions.pkl', help='Test data file') # original for test_data
@ck.option('--terms-file', '-tf', default='data/terms_gominre_trxte.pkl', help='Data file with sequences and complete set of annotations')  # original: 'data/terms_all.pkl'
@ck.option('--diamond-scores-file', '-dsf', default='data/diamond_aa.res', help='Diamond output')  # test_diamond.res / blast_output_diamond.res
@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
# new
@ck.option('--run-block-x', '-rbx', default='T', help='judge whether run block X')
@ck.option('--alpha-range', '-ar', default='25, 80, 1', type=str, help='xx')  # [25, 80, 1]   params_local['alpha_range']
@ck.option('--num-cpu-cores', '-ncc', default=60, type=int, help='xx')  # 60
@ck.option('--last-release-metadata', '-lrm', default='Alpha_last_release.json', help='xx')


def main(go_file, train_data_file, predictions_file, terms_file, diamond_scores_file, ont, run_block_x, alpha_range,
         num_cpu_cores, last_release_metadata):

    # 预处理
    go_rels = Ontology(go_file, with_rels=True)
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()
    # terms_dict = {v: i for i, v in enumerate(terms)}  # 没用上？

    train_df = pd.read_pickle(train_data_file)
    test_df = pd.read_pickle(predictions_file)  # orignal: test_data_file
    print("Length of test set: " + str(len(test_df)))

    annotations = train_df['prop_annotations'].values  # propagate annotation?????
    annotations = list(map(lambda x: set(x), annotations))
    test_annotations = test_df['prop_annotations'].values
    test_annotations = list(map(lambda x: set(x), test_annotations))
    go_rels.calculate_ic(annotations + test_annotations)


    #####################################################################
    ####################### Diamond 计算 IC ##############################
    #####################################################################
    # Print IC values of terms  # information content
    ics = {}
    for term in terms:
        ics[term] = go_rels.get_ic(term)

    prot_index = {}
    for i, row in enumerate(train_df.itertuples()):
        '''
        row:
        Pandas(Index=60795, index=437551, proteins='S6A11_HUMAN', 
        accessions='P48066; B2R6U6; Q8IYC9;', 
        sequences='MTAEKALPLGNGKAAEEARESEAPGGGCSSGGAAPARHPRVKRDKAVHERGHWNNKVEFVLSVAGEIIGLGNV。。。', 
        annotations=['GO:0005737|IEA', 'GO:0098982|IEA', 'GO:0016020|TAS', 'GO:0043005|IBA', 'GO:0005886|IBA', 。。。], 
        interpros=['IPR000175', 'IPR002982', 'IPR037272'], orgs='9606', 
        exp_annotations=['GO:0016020', 'GO:0005332', 'GO:0008028', 'GO:0051936', 'GO:0015718'], 
        prop_annotations=['GO:1905039', 'GO:0098662', 'GO:0015698', 'GO:0034220', 'GO:0098661',。。。], cafa_target=True)
        '''
        prot_index[row.proteins] = i
    #prot_index  # len = 12576， {'S6A11_HUMAN': 0, 'MBOA5_HUMAN': 1, 'SMPX_HUMAN': 2, 'IFIX_HUMAN': 3,。。。}

    # BLAST Similarity (Diamond)
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
    for i, row in enumerate(test_df.itertuples()):
        annots = {}
        prot_id = row.proteins
        # BlastKNN
        if prot_id in diamond_scores:  # 在 diamond_score中，逐一过test 503个
            sim_prots = diamond_scores[prot_id]  # len=6,
            # {'PI42A_HUMAN': 622.0, 'PI42C_HUMAN': 520.0, 'PI51B_HUMAN': 170.0, 'PI51A_HUMAN': 162.0,
            # 'PI51C_HUMAN': 158.0, 'PI5L1_HUMAN': 94.0}   加和=1726
            allgos = set()
            total_score = 0.0

            # 这下面和 Alpha_Predict.py 148 行一样
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
    ############## 引入 json 参数 #####################
    ##################################################

    # last_release_metadata = 'Alpha_last_release.json'  # 上面定义了
    with open(last_release_metadata, 'r') as f:
        last_release_data = json.load(f)
    print('111', last_release_data)


    # find alpha 执行下面的三个嵌套函数
    last_release_data['alphas'][ont] = find_alpha(ont, test_df, blast_preds, go_rels, terms, run_block_x, alpha_range, num_cpu_cores)
    print('2222', last_release_data)
    with open(last_release_metadata, 'w') as f:
        json.dump(last_release_data, f)



############################################
############## 其他函数 #####################
############################################

################ ALGORITHM TO FIND BEST ALPHAS PARAMETER ##################################

'''
三个主要函数，嵌套
find_alphas(eval_alphas(evaluate_annotations))

find_alphas：对于一定范围的alpha (e.g. range(30, 70, 2))，逐一筛过每一个alpha （0.30, 0.32, 0.34...），并把这个指定的alpha赋值给eval_alphas函数
eval_alphas: 针对alpha传来的确定的 alpha 数值，先将

evaluate_annotations: 根据 真实结果和预测结果 (labels, pred)，计算一系列的评估值f, p, r, s, ru, mi, fps, fns
'''


def evaluate_annotations(go, real_annots, pred_annots):  # (go_rels, labels, preds)
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
    # fscore, precision, recall, s?, ru?, mi?, FalsePositive, FalseNegative



def eval_alphas(alpha, ont, test_df, blast_preds, go_rels, terms, run_block_x):

    '''
    第一步：建库 deep_preds (list格式)，根据 母函数 find_alpha传来的固定alpha数值，结合deepgo & diamond(IC) 计算建库。
    [{'GO:0000045': 0.042534434720873836, 'GO:0000285': 0.03847045302391053, ...},{}。。。]    格式为: GO:score, GO:score, GO:score...
    比如test一共有100个蛋白质，待估算terms有2088个，上述list就是 有100个{}，每个{}中又2088个GO标签，
    每一个GO标签的数值score=deepgo_preds * alpha + diamond * (1-alpha)

    第二步：evaluate
    对每一个threshold （0.01, 0.02, 0.03, ... 1 一共100个），拆解上面得到的deep_preds，逐一比对其中的score
    如果score > threshold，则将该GO统计入 preds，相当于是预测出来的阳性，pred_positive。
    然后把 （labels, preds）转入函数 evaluate_annotations(go, real_annots, pred_annots)，计算统计学参数 fscore, prec, rec, s, ru, mi, fps, fns

    也就是说，这些统计学数字，是在 find_alpha指定的alpha，在eval_alphas指定的threshold下，得到的。
    筛过所有的threshold之后，在所有fscore里选出最大的fmax（不关注在哪个threshold取得的fmax）

    最后返回对应的 alpha, Fmax
    '''

    ### 第一步 建库 deep_preds：###
    deep_preds = []
    go_set = go_rels.get_namespace_terms(NAMESPACES[ont])
    go_set.remove(FUNC_DICT[ont])  # 把最根的那三个之一给删了？

    labels = test_df['prop_annotations'].values
    labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))

    # ####  FS 添加 # block X # 这里添加一下，把labels的每一个元素{GO,GO,GO},{GO,GO},{GO}和terms的2088做交集
    if run_block_x == 'T':
        # print('it is running block X !!!!!!!!!!!!!!!!!')
        terms_set = set(terms)
        labels = [s.intersection(terms_set) for s in labels]
        # print(labels)
        # print('len(labels & terms_set) = ', len(labels))


    alphas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}  # 初始化，original=0
    alphas[NAMESPACES[ont]] = alpha  # 比如ont=bp，则alphas[biological_process]=alpha外部引入
    # print(NAMESPACES[ont])
    # print(alphas)
    print('testing alpha =', alpha)  # 这个alpha 是从下面函数 find_alpha中引入的 (i/100)



    '''
    README: 这是给下面 test_df.itertuples() 中的row做注释：
    
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
    preds=array([0.0264879 , 0.00491404, 0.01317153, ..., 0.01190471, 0.01042569, 0.00605111], 
    dtype=float32))

    len(row.preds) = 2088
    '''
    # 这是核心 !!!!!!!!!!!!!!!!!!!!!!!!!!!
    for i, row in enumerate(test_df.itertuples()):
        # print(row)
        # print(row.preds)
        # print(len(row.preds))

        # 这是核心 !!!!!!!!!!!!!!!!!!!!!!!!!!!
        # score = deepgo_preds * alpha + blast * (1-alpha)
        annots_dict = blast_preds[i].copy()  # test中的第 i 个元素，在blast算的值
        for go_id in annots_dict:
            annots_dict[go_id] *= 1 - alphas[go_rels.get_namespace(go_id)]  # blast * (1-alpha)
        for j, score in enumerate(row.preds):  # deep深度学习中算的预测值preds
            go_id = terms[j]
            score *= alphas[go_rels.get_namespace(go_id)]  # 深度学习deepgo_preds * (alpha)
            if go_id in annots_dict:  # if 深度学习的go_id在diamond对比数据中
                annots_dict[go_id] += score
            else:
                annots_dict[go_id] = score

        # alpha=0.2:  'GO:0000045' = 0.6616454 * 0.01 + score = 0.042534434720873836。也可能*0.00修改alphas = {NAMESPACES['mf']: 0, NAMESPACES['bp']: 0, NAMESPACES['cc']: 0}
        # annots_dict = {'GO:0000045': 0.042534434720873836, 'GO:0000285': 0.03847045302391053, 'GO:0001505': 0.005691665941849351,...
        deep_preds.append(annots_dict)
        # deep_preds = [{'GO:0000045': 0.042534434720873836, 'GO:0000285': 0.03847045302391053, ...},{}。。。]




    ### 第二步：###
    # 上面是统计估分，下面开始 评估evaluate
    fmax = 0.0
    tmax = 0.0
    precisions = []
    recalls = []
    smin = 1000000.0  # 反正很大就是了
    rus = []
    mis = []
    for t in range(1, 101):  # t=threshold the range in this loop has influence in the AUPR output
        threshold = t / 100.0
        preds = []
        for i, row in enumerate(test_df.itertuples()):  # 对于test中逐一过每一个蛋白
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

        fscore, prec, rec, s, ru, mi, fps, fns = evaluate_annotations(go_rels, labels, preds)  # (go, real_annots, pred_annots) ！！！！！！！！！！！！！！！！！！
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
    # 因为输出为三者加权 np.sum([smin, -fmax, -aupr])，越小越好，所以分别除以15,40,10相当于是归一化
    if ont == 'mf':
        smin /= 15
    elif ont == 'bp':
        smin /= 40
    elif ont == 'cc':
        smin /= 10

    return alpha, -fmax  # original: alpha, np.sum([smin, -fmax, -aupr])  这里修改后，只根据fmax来判断选哪个alpha


def find_alpha(ont, test_df, blast_preds, go_rels, terms, run_block_x, alpha_range, num_cpu_cores):
    '''
    逐一过所有的alpha （范围：alpha_range）
    并将每一个alpha传到 eval_alphas 函数中，得到最大的 Fmax
    每一个 (alpha, Fmax) 组成 results，最后选出最大的那个Fmax对应的alpha

    '''
    extra = [ont, test_df, blast_preds, go_rels, terms, run_block_x]
    # inputs = range(25, 80, 1)  ## ????? alpha 只能从？？？？？这个范围找么？？？？  # original range(45, 75, 1)

    # list2range: '25, 80, 1' to range(25, 80, 1)
    str_values = alpha_range.split(',')
    int_values = [int(value.strip()) for value in str_values]
    inputs = range(*int_values)  # range(25, 80, 1)
    print('inputs=', inputs)  # FS
    num_cores = num_cpu_cores  # 60  # most cpu can be used

    results = Parallel(n_jobs=num_cores)(delayed(eval_alphas)(i / 100, *extra) for i in inputs)
    # eval_alphas是个函数，定义在上面。第一个参数 i/100是eval_alphas中的第一个参数 alpha
    # 该函数以i / 100作为第一个参数，后跟extra列表中的所有元素
    # i的取值从20到79，步长为1，通过并行化处理，将这些函数调用的结果存储在results列表中

    print('22222222222222222222222222')
    print(results)
    # results = [(0.2, -0.2244738075593926), (0.21, -0.2323958248068838), (0.22, -0.2419488147142912), , (0.23, -0.2503843465951797)...


    chosen = min(results, key=lambda x: x[1])
    # 在results列表中选择具有最小值的元素。这里使用了lambda表达式作为key函数，将列表中的每个元素x的第二个元素作为比较的依据。
    print(chosen)
    # chosen = (0.54, -0.33197674380787534)

    return chosen[0]  # 返回chosen列表的第一个元素作为函数的结果。






#####################################################################################








if __name__ == '__main__':
    alphas = {'mf': 0, 'bp': 0, 'cc': 0}
    main()

