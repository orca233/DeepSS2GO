#!/usr/bin/env python

# README
#
# 固定threshold, 固定ont=bp
# 摘出 TP FP FN

# 这里的labels & preds 就是这个step4_analysis想要的，各有662个(test中protein的个数)


import numpy as np
import pandas as pd
import click as ck

import logging
from sklearn.metrics import roc_curve, auc, matthews_corrcoef

import math
from utils import FUNC_DICT, Ontology, NAMESPACES

from step0_TrainTestSetting import *
from collections import Counter
from functools import reduce


# import configparser
# config_local = configparser.ConfigParser()  # 创建 ConfigParser 对象
# config_local.read('config_local.ini')  # 读取配置文件
# # read param from config_local.ini
# # 循环读取每个配置项并添加到字典中
# config_local_dict = {}
# for section in config_local.sections():
#     section_dict = {}
#     for option in config_local.options(section):
#         section_dict[option] = config_local.get(section, option)
#     config_local_dict[section] = section_dict



logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

@ck.command()
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='Data file with training features')
@ck.option('--test-data-file', '-tsdf', default='data/predictions.pkl', help='Test data file')
@ck.option('--terms-file', '-tf', default='data/terms.pkl', help='Data file with sequences and complete set of annotations')
# @ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--subontology', '-sot', default='???', help='GO subontology (bp, mf, cc)')
@ck.option('--go-file', '-gf', default=path_pub_data+'go.obo', help='Gene Ontology file in OBO Format')  # FS 添加

def main(train_data_file, test_data_file, terms_file, subontology, go_file):  # FS 添加go_file




    go_rels = Ontology(go_file, with_rels=True)
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()
    # print(terms)  # numpy.ndarray， 2088个， ['GO:0006275' 'GO:0002698' 'GO:0044843' ...  'GO:0031327' 'GO:0048856']

    train_df = pd.read_pickle(train_data_file)
    test_df = pd.read_pickle(test_data_file)  # predictions.pkl的信息
    print("Len(test_set)= ", len(test_df))

    annotations = train_df['prop_annotations'].values  # train_data的GO数据集
    annotations = list(map(lambda x: set(x), annotations))  # 变list
    test_annotations = test_df['prop_annotations'].values  # test_data的GO数据集
    test_annotations = list(map(lambda x: set(x), test_annotations))


    ##### 这个block 计算ic值，生成字典 ######
    go_rels.calculate_ic(annotations + test_annotations)  # 计算information content
    # Print IC values of terms
    ics = {}
    for term in terms:
        ics[term] = go_rels.get_ic(term)
    # print('11111111111111  ics:')
    # print(ics)
    ##### 这个block 计算ic值，生成字典 ---- done ######




    ##### 这个模块生成 train & test 所有注释annotation的统计，并从大到小排成字典 ######
    # annotations + test_annotations = [{GO, GO, GO..}, {GO..}, {}, {}..]  # len=13238
    train_test_annot = []
    for i in (annotations + test_annotations):  # mix all GO into one list,
        train_test_annot.extend(list(i))  # final len = 989741

    def annot_count_to_dict(my_list):
        # 使用Counter计数，并转换为字典格式
        count_dict = dict(Counter(my_list))
        dict_with_list = {key: count_dict[key] for key in count_dict}  # 将字典的值转换为列表形式。
        sorted_dict = dict(
            sorted(dict_with_list.items(), key=lambda x: x[1], reverse=True))  # 按值列表中的第一个元素从大到小排序的字典 x[1][0]
        return sorted_dict

    train_test_annot_dict = annot_count_to_dict(train_test_annot)
    # print(train_test_annot_dict)
    ##### 这个模块生成 train & test 所有注释annotation的统计，并从大到小排成字典 -- done ######



    ###### 这个模块用于找到指定GO term到祖先根节点的最短路径并返回该值 #####
    # annotations + test_annotations = [{GO, GO, GO..}, {GO..}, {}, {}..]  # len=13238
    train_test_annot = []
    for i in (annotations + test_annotations):  # mix all GO into one list,
        train_test_annot.extend(list(i))  # final len = 989741
    train_test_annot_set = set(train_test_annot)  # 去冗余，len = 17265

    short_path_to_root = {}
    for term in train_test_annot_set:
        short_path_to_root[term] = go_rels.calculate_shortest_path_to_root(term)

    short_path_to_root_sorted = dict(
        sorted(short_path_to_root.items(), key=lambda x: x[1], reverse=True))  # 按值列表中的第一个元素从大到小排序的字典 x[1][0]
    # print(short_path_to_root_sorted)  # len = 17265

    ###### 这个模块用于找到指定GO term到祖先根节点的最短路径并返回该值 -- done #####





    ##########################################
    ########## param:  ######################
    ##########################################

    '''
    output/test_TrainHUMAN_TestHUMAN_KernelsXFiltersY_aa/DeepSS2GO_Kernel32_Filter512
    
    ==> results/Fmax_AUPR_Smin_bp.txt <==
    Length of test set: 662  threshold: 0.16
    Smin: 53.027
    Fmax: 0.396
    AUPR: 0.342
    
    ==> results/Fmax_AUPR_Smin_cc.txt <==
    Length of test set: 662  threshold: 0.32
    Smin: 12.649
    Fmax: 0.651
    AUPR: 0.691
    
    ==> results/Fmax_AUPR_Smin_mf.txt <==
    Length of test set: 662  threshold: 0.15
    Smin: 16.353
    Fmax: 0.422
    AUPR: 0.378
    '''

    threshold = 0.05  # 固定
    # subontology = ['bp', 'cc', 'mf']
    subontology = ['cc']
    TopNB = 50

    # test_df = test_df[0:3]  # 测试前n行，正式运行时注释掉-------------------------------


    print('########### param ############')
    print('threshold = ', threshold)
    print('subontology = ', subontology)
    print('TopNB = ', TopNB)


    ########## param: done ##############




    # DeepGOPlus
    for ont in subontology:
        print('\n ----------------- ont = %s -----------------------\n' % ont)
        go_set = go_rels.get_namespace_terms(NAMESPACES[ont])  # 所有在当前'bp'下的GO？
        go_set.remove(FUNC_DICT[ont])  # 除去祖师爷辈的

        ###### ？？？ 下面的labels 为什么要用prop_annotation呢？应该用terms的2088集合啊?因为prop的还包含term出现频率小于50的？
        # 回答：因为test_df用的是'predictions.pkl'，所以这里面的prop_annotation其实已经是筛出的超过50次的terms的了，r u sure?
        ###### 前面这个prop_annotations一共4266条，里面很多因为长度和出现次数小被删掉了啊
        labels = test_df['prop_annotations'].values  # 这个label是真实值，是实际出现在test_df.prop_anno中的GO term。格式：list里的list
        # print(labels, '\n')
        # [list(['GO:0008150', 'GO:0008641', 'GO:0023057', 'GO:0046627',... 'GO:0140657'])
        # list(['GO:0032501', 'GO:0042770', 'GO:0008150'... ])
        # list... list...]  len=662

        labels = list(map(lambda x: set(filter(lambda y: y in go_set, x)), labels))  # go_set是bp/cc/mf中的一个子集。格式：list里的set， go_set go_set go_set go_set go_set ！！！！！
        # print(labels, '\n')
        # 这里每一个{} set的长度都不一样，根据实际情况而定
        # [{'GO:0009893', 'GO:0090217', 'GO:0008286', 'GO:0006661'...,}
        # {'GO:0048856', 'GO:0006950', 'GO:0050789', ..}, ...{}...
        #  set()... ]   为啥是空的，一个set ? 答：因为这个蛋白质没有对应的bp的子集，就空了



        ####################################################################################################
        ############ 如果不加下面这个block，估计会大大增加TN? 因为没有加入label的会score很低，进而将小众的term标为N，成为TN #####
        ####################################################################################################

        print('\n ----- count how many types of GO ------')

        ###### labels包含662个元素（set），每个元素是一个蛋白质的GO集合。把所有出现的terms重新组成一个set，
        ###### 这个模块仅用于统计个数，经过下一步的block x后看看还剩下多少 ######
        labels_1set = reduce(lambda a, b: a.union(b), labels)  # 不断从labels中提取新数据求并集
        # print(labels_1set)
        print('len(labels_1set) before X terms -- contain<50 times = ', len(labels_1set))  # bp=4266,cc=620,mf=1046, total=5932 ---------------------


        # block X # 这里添加一下，把labels的每一个元素{GO,GO,GO},{GO,GO},{GO}和terms的2088做交集
        terms_set = set(terms)
        labels = [s.intersection(terms_set) for s in labels]
        # print(labels)
        print('len(labels & terms_set) -- test.pkl prot nb = ', len(labels))  # 662 是test中蛋白质个数


        ###### labels包含662个元素（set），每个元素是一个蛋白质的GO集合。把所有出现的terms重新组成一个set，
        ###### 这个模块仅用于统计个数，经过下一步的block x后看看还剩下多少 ######
        labels_1set_again = reduce(lambda a, b: a.union(b), labels)  # 不断从labels中提取新数据求并集
        # print(labels_1set_again)
        print('len(labels_1set_again) after X terms = ', len(labels_1set_again))  # bp=1516,cc=268,mf=301, total=2085+3=2088 ---------------------

        ####################################################################################################
        ###################################### done #######################################################
        ####################################################################################################




        print('\n ------ calc preds ------')

        deep_preds = []
        for i, row in enumerate(test_df.itertuples()):  # 这是怎么个循环的。。。循环662遍，对不？
            annots_dict = {}
            for j, score in enumerate(row.preds):
                go_id = terms[j]
                if go_id in annots_dict:  # annots_dict 是从blast/diamond初始化来的?
                    annots_dict[go_id] += score  # 这是怎么积分的？两个简单相加？
                else:
                    annots_dict[go_id] = score
            deep_preds.append(annots_dict)

        # deep_preds是从predictions.pkl （test）生成的list，包含662个{}，这是test中protein个数
        # 每个{annots_dict}包含2088个GO & 对应的分数，这是公共terms的个数
        # [{'GO:0046906': 6.432776e-19, 'GO:0006351': 0.94251436, 'GO:0006790': 5.85886e-09, ...}
        # {'GO:0046906': 2.926681e-20, 'GO:0006351': 0.08532268, 'GO:0006790': 5.5335445e-08, ...}
        # {'GO:0046906': 2.810078e-27, 'GO:0006351': 0.15612392, 'GO:0006790': 1.5042304e-11,...}
        # {}, {}...]


        fmax = 0.0
        tmax = 0.0
        precisions = []
        recalls = []
        smin = 1000000.0
        rus = []
        mis = []

        # for t in range(1, 101):  # the range in this loop has influence in the AUPR output
        #     threshold = t / 100.0




        preds = []
        for i, row in enumerate(test_df.itertuples()):
            # row如下：
            # Pandas(Index=35003,   index=214259,    proteins='KS6A6_HUMAN'
            # accessions='Q9UK32; B2R854; B7ZL90; Q6FHX2; Q8WX28; Q9H4S6;',...
            # sequences='MLPFAPQDEPWDREMEVFS...
            # annotations=['GO:0005737|IBA', 'GO:0005829|IEA',...
            # interpros=['IPR000961', 'IPR011009', 'IPR017892', 'IPR000719', 'IPR017441'...
            # orgs='9606'
            # exp_annotations=['GO:0001650', 'GO:0005739'...
            # prop_annotations=['GO:0032501', 'GO:0042770',...
            # cafa_target=True,
            # labels=array([0, 0, 0, ..., 0, 0, 0], dtype=int32),
            # preds=array([6.4327758e-19, 9.4251436e-01, 5.8588601e-09, ..., 1.5411665e-14,...


            # 这是针对某一个蛋白质的 注释集合term
            annots = set()
            for go_id, score in deep_preds[i].items():  # 一共662行，逐行i看
                if score >= threshold:  # 如果当前评分 > threshold，则判断为positive，收录在annots集合里
                    annots.add(go_id)

            # print('annots = ', annots)
            # print(len(annots))

            # 把父辈也揪出来，估计是因为部分父辈的GO_term因为频率小于50前面被clean了？
            new_annots = set()
            for go_id in annots:
                new_annots |= go_rels.get_ancestors(go_id)  # 求并集，把祖先父类也加进来，相当于把model没有预测出祖先的也给揪出来
            preds.append(new_annots)
            # print(len(new_annots))  比钱面的len(annots)稍微多出0-4个terms

        # preds 不断append新元素，循环结束,每一个set里是“预测”的GO+祖先，因为前面的test_propogation里面也把祖先算进来了
        # preds= [{'GO:0010556', 'GO:0071944', 'GO:0022607',...}, {GO,GO,GO}, {},...{}]
        # preds 应该不需要 block X来进行和2088个terms做交集了，因为预测的池子就是这个terms集合


        ##########################################
        ###### 这个模块仅用于统计个数， ##############
        ##########################################
        preds_set = reduce(lambda a, b: a.union(b), preds)  # 不断从labels中提取新数据求并集
        print('len(preds_set) = ', len(preds_set))  # 恒定=997, 不变，因为这是预测的，包含三大类别

        # Filter classes  过滤种类，筛除只存在于go_set--三大分类其中之一，      。。。 这是源代码内容。。。。
        preds = list(map(lambda x: set(filter(lambda y: y in go_set, x)), preds))  # go_set go_set go_set go_set go_set go_set go_set go_set ！！！！！！

        ###### 这个模块仅用于统计个数， ######
        preds_set_again = reduce(lambda a, b: a.union(b), preds)  # 不断从labels中提取新数据求并集
        print('len(preds_set_again) = ', len(preds_set_again))  # 筛选后bp=649,cc=187, mf=158, total = 994 预测的 -----------

        # preds = [{'GO:2001141', 'GO:0051171', 'GO:0051649',...},{GO,GO,},...}]    len=662
        ######## end ##############################


        fscore, prec, rec, s, ru, mi, fps, fns, tp_lst, fp_lst, fn_lst = evaluate_annotations(go_rels, labels, preds)
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




        ##################################
        ### FS change ###  真正的分析开始啦
        ##################################
        print('\n------------------------- FS CHANGE start ------------------------------')

        # # tp_lst因为有很多重复的GO，所以要count，看哪个重复最多。下面这几个作用不大
        # print('tp_lst= %d' % len(tp_lst))  # 9774 ['GO:0044238', 'GO:0031325', 'GO:0044249'...]
        # print('fp_lst= %d' % len(fp_lst))  # 14741
        # print('fn_lst= %d' % len(fn_lst))  # 18223

        # # 分别把662个test中的 tp，fp， fn 取交集，看看都有哪些错了。求交集，初始化不能为空set()，为lst中的第一个
        # # 。。。还不能单纯去交集，否则所有的一交集，就空了。。。
        # # 应该先把所有出现的GO添加成lst，然后用count，看看哪个出现的多


        # 统计 tp, fp, fn各出现多少次count
        print('\n------ tp fp fn __ {GO:[count, ic], GO, GO...} ----------- \n')


        '''
        方法二，count转字典，并从大到小排序
        对列表中的元素进行计数，将结果转换为字典格式，其中字典的值用列表表示。列表的第一个元素为计数的数字，第二个元素固定为{ics}字典ic值。
        然后，按照字典的值列表中的第一个元素从大到小排序，并展示前top_number个（前100个）也可以选“前无限大个”
        '''

        def count_to_dict_TopNB(my_list, TopNB):
            # 使用Counter计数，并转换为字典格式，key=GOterm, value=[cnt_nb_1, ic, cnt_nb_2], cnt_nb_1=在预测tp/fp/fn中重复次数，cnt_nb_2=在所有train&test数据中重复次数
            count_dict = dict(Counter(my_list))
            # 将字典的值转换为列表形式。 value中，
            # 第0位：预测tp/fp/fn中重复次数，
            # 第1位：ics[key]为该term对应ic值，在上面的字典{ics}可看。
            # 第2位：是在所有train&test数据中重复次数
            # 第3位：距离根节点的最短距离
            dict_with_list = {key: [count_dict[key], ics[key], train_test_annot_dict[key], short_path_to_root_sorted[key]] for key in count_dict}
            sorted_dict = dict(sorted(dict_with_list.items(), key=lambda x: x[1][0], reverse=True))  # 按值列表中的第一个元素从大到小排序的字典
            TopNB_dict = dict(list(sorted_dict.items())[:TopNB])  # 展示前TopNB个结果
            return TopNB_dict



        tp_dict = count_to_dict_TopNB(tp_lst, TopNB)
        print('tp_dict:::::::::: \n', tp_dict, '\n')

        fp_dict = count_to_dict_TopNB(fp_lst, TopNB)
        print('fp_dict:::::::::: \n', fp_dict, '\n')

        fn_dict = count_to_dict_TopNB(fn_lst, TopNB)
        print('fn_dict:::::::::: \n', fn_dict, '\n')



        # 判断交叉数量： tp & fp,   以及  tp & fn
        print('\n------ cross repeat::: tp_fp, tp_fn, fp_fn  ----------- \n')

        tp_x_fp = set(tp_dict.keys()) & set(fp_dict.keys())
        print(tp_x_fp)
        print('tp & fp = %d \n' % len(tp_x_fp))

        tp_x_fn = set(tp_dict.keys()) & set(fn_dict.keys())
        print(tp_x_fn)
        print('tp & fn = %d \n' % len(tp_x_fn))

        fp_x_fn = set(fp_dict.keys()) & set(fn_dict.keys())
        print(fp_x_fn)
        print('fp & fn = %d \n' % len(fp_x_fn))

        print('threshold=', threshold)



        '''
        分析
        'GO:0032502是最大的fp，所以
        手动看seq，分析这个蛋白质seq是保守的还是非保守的，不同蛋白质的这个变化大不大
        这个GO也出现在了tp中，也就是这个threshold的分数其实很敏感
        把前100个fp & fn的GO拿出来，组成新的terms2.pkl，旧模型判断排除掉terms2的案例，用一个新模型进行判断这些新的terms2 的案例
        主要两个思路：1）分析其seq&生物学意义，2）单纯用统计学继续迭代处理。
        生物分析：写个代码，把所有包含fp前100的prot列出来，看看他们是怎么被误判的，

        效果并不好，bp_Fmax = 0.231, AUPR=0.126
        原因可能很多，这些fp和fn其实有很多属于父类GO，删除掉也会影响tp的判断
        解决办法，用全部的训练模型，只对  瘦身后的test.terms进行测试？ -- 那就前后矛盾了，因为前面说模型其实是需要调整的，这里是说模型很好
        tp & fp  threshold = 0.1, 50个有41个重复的。。。
        thre = 0.2，50个有50个重复

        threshold=0.01， top_nb=100,的话，重新训练在step3，最大的fmax的threshold=0.64了，太大了，
        还是不行，GO:0000003也在fp里，说明只要判断错了bp，一定会揪出GO:0000003来。这是一个很高的父类，奇怪的是，GO:0000003没有再tp里？
        尝试把top_nb减少一些，

        有个新想法，因为tp, fp, fn正常排序是最大的在前面，一般会是父类
        那我们倒叙排列，看看哪些 子类 被分错了

        如果不选前100，把所有的tp的term给放到下一轮进行重新判断，937个，看看准确率？估计新的threshold会很大？
        '''

        ### FS change done ###

    print('------------ happy endding -----------------')
























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
    tp_lst = []  # 下面循环会把每个tp_set添加近该list
    fp_lst = []
    fn_lst = []
    for i in range(len(real_annots)):  # 这是针对一个蛋白质的一系列 真实GO & 预测GO的分析，拿到tp, fp, fn
        if len(real_annots[i]) == 0:
            continue
        tp = set(real_annots[i]).intersection(set(pred_annots[i]))  # true positive
        fp = pred_annots[i] - tp  # false positive
        fn = real_annots[i] - tp  # false negative

        ### FS change ###
        tp_lst.extend(list(tp))  # 不能用append，用extend，否则是多个list并列成一个二维的list [[xx,xx],[xx,xxx], ...  ]
        fp_lst.extend(list(fp))
        fn_lst.extend(list(fn))
        ### FS change done ###


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
    return f, p, r, s, ru, mi, fps, fns, tp_lst, fp_lst, fn_lst
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



