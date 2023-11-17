#!/usr/bin/env python

'''
引自deepgoplus中的predict.py
很多函数引自 step2_Test.py

需要model_checkpoint.pth，对test_data.fa/pkl进行预测

如果是全新的未知的测试，一般会给出aa.fa文件，
if 用aa预测，转成pkl，准备好test_data.fa & test_data.pkl即可
elif 用ss8测试，转成ss8.fa，在转成pkl，再准备test_data.fa & test_data.pkl

'''

import click as ck
import numpy as np
import pandas as pd
import torch.nn as nn
import time
from utils import Ontology, NAMESPACES
# from sklearn.metrics import roc_curve, auc, matthews_corrcoef
import torch
from torch.utils.data import Dataset, DataLoader
import json

# 在本轮的step1中提前把 step0_TrainTestSetting_local_aa/ss8 拷贝好
from step0_TrainTestSetting_local_aa import params_local as params_local_aa
from step0_TrainTestSetting_local_ss8 import params_local as params_local_ss8
from step0_TrainTestSetting_global import path_base


MAXLEN = params_local_aa['MAXLEN']

@ck.command()
@ck.option('--in-file-aa', '-if', default='data/test_data_aa.fa', help='Input FASTA file', required=True)  # 提前准备好
@ck.option('--in-file-ss8', '-if', default='data/test_data_ss8.fa', help='Input FASTA file', required=True)  # 提前准备好

# @ck.option('--out-file', '-of', default='data/results.csv', help='Output result file')
@ck.option('--out-file-bp', '-of', default='data/results_bp.csv', help='Output result file')
@ck.option('--out-file-cc', '-of', default='data/results_cc.csv', help='Output result file')
@ck.option('--out-file-mf', '-of', default='data/results_mf.csv', help='Output result file')

@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
# @ck.option('--go-file', '-gf', default=path_base + 'pub_data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
@ck.option('--model-file-aa', '-mf', default='data/model_checkpoint_aa.pth', help='Tensorflow model file')
@ck.option('--model-file-ss8', '-mf', default='data/model_checkpoint_ss8.pth', help='Tensorflow model file')

@ck.option('--terms-file-aa', '-tf', default='data/terms_gominre_trxte_aa.pkl', help='List of predicted terms')  # 这个是从s2_TrainTest/step1中，结合train&test_data交叉的到的
@ck.option('--terms-file-ss8', '-tf', default='data/terms_gominre_trxte_ss8.pkl', help='List of predicted terms')  # 这个是从s2_TrainTest/step1中，结合train&test_data交叉的到的

@ck.option('--annotations-file', '-tf', default='data/train_data_aa.pkl', help='Experimental annotations')
@ck.option('--chunk-size', '-cs', default=100000, help='Number of sequences to read at a time')  # original 1000
@ck.option('--diamond-file', '-df', default='data/diamond_aa.res', help='Diamond Mapping file')
@ck.option('--threshold', '-t', default=0.1, help='Prediction threshold')
@ck.option('--batch-size', '-bs', default=32, help='Batch size for prediction model')
# @ck.option('--maxlen', '-ml', default=params_local['MAXLEN'])
# new
@ck.option('--test-data-file-aa', '-tedf', default='data/test_data_aa.pkl', help='XX')  # 提前准备好
@ck.option('--test-data-file-ss8', '-tedf', default='data/test_data_ss8.pkl', help='XX')  # 提前准备好

# @ck.option('--aa-ss', '-aass', default=params_local['aa_ss'], help='aa/ss8/ss3')
@ck.option('--prot-letter-aa', '-plaa', default=params_local_aa['PROT_LETTER_aa'], type=list, help='XX')
@ck.option('--prot-letter-ss8', '-plss8', default=params_local_aa['PROT_LETTER_ss8'], type=list, help='XX')
@ck.option('--prot-letter-ss3', '-plss3', default=params_local_aa['PROT_LETTER_ss3'], type=list, help='XX')


@ck.option('--kernels-list-aa', '-klaa', default=params_local_aa['kernels'], type=list, help='XX')      # # kernels & filters 需要是list形式，即使只有一个 int，也要转成list
@ck.option('--kernels-list-ss8', '-klss8', default=params_local_ss8['kernels'], type=list, help='XX')      # # kernels & filters 需要是list形式，即使只有一个 int，也要转成list

@ck.option('--filters-list-aa', '-flaa', default=params_local_aa['filters'], type=list, help='XX')
@ck.option('--filters-list-ss8', '-flss8', default=params_local_ss8['filters'], type=list, help='XX')



@ck.option('--ont', '-o', default='mf', help='GO subontology (bp, mf, cc)')
@ck.option('--alpha', '-a', default='json', help='alpha = json(with quote) or 0-1(without quote, eg 0.3 float)')  # 如果alpha='json'，则采用json数据，否则alpha=数字，或外来click引入
@ck.option('--beta', '-b', default='json', help='beta = json(with quote) or 0-1(without quote, eg 0.3 float)')  # 如果beta='json'，则采用json数据，否则beta=数字，或外来click引入


# @ck.option('--maxlen', '-fl', default=params_local['filters'], type=list, help='XX')
# @ck.option('--filters-list', '-fl', default=params_local['filters'], type=list, help='XX')


def main(in_file_aa, in_file_ss8, out_file_bp, out_file_cc, out_file_mf, go_file, model_file_aa, model_file_ss8,
         terms_file_aa, terms_file_ss8,  # out_file
         annotations_file, chunk_size, diamond_file, threshold, batch_size,
         test_data_file_aa, test_data_file_ss8, prot_letter_aa, prot_letter_ss8, prot_letter_ss3,
         kernels_list_aa, kernels_list_ss8, filters_list_aa, filters_list_ss8, ont, alpha, beta):  # maxlen， aa_ss


    #######################################################################
    ##################### pre_setting 数据和参数的预处理 ######################
    ########################################################################

    ##### 根据 params_local的aa/ss8/ss3确定PROT_LETTER，kernels_tuple #####
    # PROT_LETTER = []
    # PROT_LETTER_len = -1
    # if aa_ss == 'aa':
    #     PROT_LETTER = prot_letter_aa
    #     PROT_LETTER_len = len(PROT_LETTER)  # =20
    # elif aa_ss == 'ss8':
    #     PROT_LETTER = prot_letter_ss8
    #     PROT_LETTER_len = len(PROT_LETTER)  # =8
    # elif aa_ss == 'ss3':
    #     PROT_LETTER = prot_letter_ss3
    #     PROT_LETTER_len = len(PROT_LETTER)  # =3
    #
    # PROT_INDEX = dict()
    # for i in range(len(PROT_LETTER)):
    #     PROT_INDEX[PROT_LETTER[i]] = i + 1


    # only for PredictAlphaBeta:
    PROT_LETTER_aa = prot_letter_aa  # list
    PROT_LETTER_len_aa = 20
    PROT_INDEX_aa = dict()
    for i in range(len(PROT_LETTER_aa)):
        PROT_INDEX_aa[PROT_LETTER_aa[i]] = i + 1

    PROT_LETTER_ss8 = prot_letter_ss8  # list
    PROT_LETTER_len_ss8 = 8
    PROT_INDEX_ss8 = dict()
    for i in range(len(PROT_LETTER_ss8)):
        PROT_INDEX_ss8[PROT_LETTER_ss8[i]] = i + 1


    kernels_tuple_aa = [(value, PROT_LETTER_len_aa + 1) for value in kernels_list_aa]
    kernels_tuple_ss8 = [(value, PROT_LETTER_len_ss8 + 1) for value in kernels_list_ss8]




    go = Ontology(go_file, with_rels=True)  # 这里的go即 alpha_EvaluateAlpha.py中的go_rels
    # test_data_file = 'data/test_data.pkl'

    terms_df_aa = pd.read_pickle(terms_file_aa)  # 这个是s2_step1中 train_data X test_data 得到的
    terms_aa = terms_df_aa['terms'].values.flatten()
    terms_dict_aa = {v: i for i, v in enumerate(terms_aa)}
    terms_classes_aa = len(terms_dict_aa)

    terms_df_ss8 = pd.read_pickle(terms_file_ss8)  # 这个是s2_step1中 train_data X test_data 得到的
    terms_ss8 = terms_df_ss8['terms'].values.flatten()
    terms_dict_ss8 = {v: i for i, v in enumerate(terms_ss8)}
    terms_classes_ss8 = len(terms_dict_ss8)


    ############################################################################
    ############ 第一步：计算diamond的IC 还是别的啥，反正和IC相关? ##################
    ############################################################################


    # 计算IC ? # Read known experimental annotations
    annotations = {}
    df = pd.read_pickle(annotations_file)
    for row in df.itertuples():
        annotations[row.proteins] = set(row.prop_annotations)

    go.calculate_ic(annotations.values())
    diamond_preds = {}
    mapping = {}
    with open(diamond_file, 'rt') as f:
        for line in f:
            it = line.strip().split()
            if it[0] not in mapping:
                mapping[it[0]] = {}
            mapping[it[0]][it[1]] = float(it[2])

    for prot_id, sim_prots in mapping.items():
        annots = {}
        allgos = set()
        total_score = 0.0

        # 这下面和Alpha_FindAlpha.py 102行一样
        for p_id, score in sim_prots.items():
            allgos |= annotations[p_id]
            total_score += score

        allgos = list(sorted(allgos))
        sim = np.zeros(len(allgos), dtype=np.float32)

        for j, go_id in enumerate(allgos):
            s = 0.0
            for p_id, score in sim_prots.items():
                if go_id in annotations[p_id]:
                    s += score
            sim[j] = s / total_score
        for go_id, score in zip(allgos, sim):
            annots[go_id] = score
        diamond_preds[prot_id] = annots

    print('22222222')

    print(type(diamond_preds))
    # print(diamond_preds)




    ################################################################################################
    ############ 第二步：通过check_model.pth 计算新的(未知的)test_data.pkl，拿到结果 XXX.XXX ##################
    ################################################################################################

    ###################### 第二步中的 aa #########################

    # 创建ProteinDataset实例，不同于train，这里只有test的实例
    test_df_aa = pd.read_pickle(test_data_file_aa)  # 这里需要test_data.pkl文件，重新进行类似于step3_Test.py的步骤，根据训练好的model再预测
    test_dataset_aa = ProteinGODataset(test_df_aa, terms_dict_aa, PROT_LETTER_len_aa, PROT_INDEX_aa)
    test_dataloader_aa = DataLoader(test_dataset_aa, batch_size=batch_size, shuffle=False)


    # 判断用 single_GPU / multi_GPU 并行运算。      创建model模型实例
    device_ids_aa = [0, 1] # 这里就直接指定了，aa用[0,1], ss8用[2,3]. params_local_aa['device_ids']  # cuda:3"   or    [0, 1, 2, 3]   特注：aa用[0, 1], ss8用[2, 3]

    if isinstance(device_ids_aa, str):  # device_ids是个str="cuda:3"，用单GPU
        print('--- single GPU = ', device_ids_aa)
        device_aa = torch.device(device_ids_aa if torch.cuda.is_available() else "cpu")  # 将模型移动到cuda:2
        model_aa = Model1(terms_classes_aa, params_local_aa, kernels_tuple_aa, filters_list_aa)  # 创建Model1模型实例
        model_aa.to(device_aa)  # 将模型移动到GPU（单个或多个）

    else:  # device_ids是个list=[0, 1, 2, 3]，用多GPU，使用DataParallel进行并行计算
        print('--- multi GPU = ', device_ids_aa)
        # device = torch.device("cuda" if torch.cuda.is_available() else "cpu")  # 对于device_ids = [2, 1]这种有问题，猜的啊：可能是吧model.to(device)默认扔给cuda:0了
        device_aa = torch.device("cuda:{}".format(device_ids_aa[0]) if torch.cuda.is_available() else "cpu")  # 将device_ids[0]的值动态地插入到字符串中，生成一个完整的字符串，例如"cuda:2"

        model_aa = Model1(terms_classes_aa, params_local_aa, kernels_tuple_aa, filters_list_aa)  # 创建Model1模型实例
        model_aa.to(device_aa)  # 将模型移动到GPU（单个或多个）
        model_aa = nn.DataParallel(model_aa, device_ids=device_ids_aa)  # 将模型包装在DataParallel中，实现多GPU的并行计算
        model_aa.to(device_aa)  # 将包装后的模型再次移动到指定设备，以确保模型及其参数和缓冲区都位于正确的设备上

    # 加载模型 ??
    model_aa.load_state_dict(torch.load(model_file_aa))
    model_aa.to(device_aa)
    model_aa.eval()




    ###################### 第二步中的 ss8 #########################

    # 创建ProteinDataset实例，不同于train，这里只有test的实例
    test_df_ss8 = pd.read_pickle(test_data_file_ss8)  # 这里需要test_data.pkl文件，重新进行类似于step3_Test.py的步骤，根据训练好的model再预测
    test_dataset_ss8 = ProteinGODataset(test_df_ss8, terms_dict_ss8, PROT_LETTER_len_ss8, PROT_INDEX_ss8)
    test_dataloader_ss8 = DataLoader(test_dataset_ss8, batch_size=batch_size, shuffle=False)

    # 判断用 single_GPU / multi_GPU 并行运算。      创建model模型实例
    device_ids_ss8 = [2, 3] # 这里就直接指定了，aa用[0,1], ss8用[2,3]. params_local_ss8['device_ids']  # cuda:3"   or    [0, 1, 2, 3]   特注：aa用[0, 1], ss8用[2, 3]

    if isinstance(device_ids_ss8, str):  # device_ids是个str="cuda:3"，用单GPU
        print('--- single GPU = ', device_ids_ss8)
        device_ss8 = torch.device(device_ids_ss8 if torch.cuda.is_available() else "cpu")  # 将模型移动到cuda:2
        model_ss8 = Model1(terms_classes_ss8, params_local_ss8, kernels_tuple_ss8, filters_list_ss8)  # 创建Model1模型实例
        model_ss8.to(device_ss8)  # 将模型移动到GPU（单个或多个）

    else:  # device_ids是个list=[0, 1, 2, 3]，用多GPU，使用DataParallel进行并行计算
        print('--- multi GPU = ', device_ids_ss8)
        # device = torch.device("cuda" if torch.cuda.is_available() else "cpu")  # 对于device_ids = [2, 1]这种有问题，猜的啊：可能是吧model.to(device)默认扔给cuda:0了
        device_ss8 = torch.device("cuda:{}".format(device_ids_ss8[0]) if torch.cuda.is_available() else "cpu")  # 将device_ids[0]的值动态地插入到字符串中，生成一个完整的字符串，例如"cuda:2"

        model_ss8 = Model1(terms_classes_ss8, params_local_ss8, kernels_tuple_ss8, filters_list_ss8)  # 创建Model1模型实例
        model_ss8.to(device_ss8)  # 将模型移动到GPU（单个或多个）
        model_ss8 = nn.DataParallel(model_ss8, device_ids=device_ids_ss8)  # 将模型包装在DataParallel中，实现多GPU的并行计算
        model_ss8.to(device_ss8)  # 将包装后的模型再次移动到指定设备，以确保模型及其参数和缓冲区都位于正确的设备上

    # 加载模型 ??
    model_ss8.load_state_dict(torch.load(model_file_ss8))
    model_ss8.to(device_ss8)
    model_ss8.eval()



    '''
    光加载了，还没预测呢？？？？？？？？？？？？？？？？？？？？？？
    -后面会用这个model进行预测！！！！！！！！！！！！！！
    '''




    #################################################
    ########## 第三部分：整合 #########################
    #################################################

    # (1 - alpha - beta) * diamond + alpha * preds_aa + beta * preds_ss8

    # alphas = {NAMESPACES['mf']: 0.48, NAMESPACES['bp']: 0.4, NAMESPACES['cc']: 0.49}
    alphas = {NAMESPACES['bp']: 0, NAMESPACES['cc']: 0, NAMESPACES['mf']: 0}
    betas = {NAMESPACES['bp']: 0, NAMESPACES['cc']: 0, NAMESPACES['mf']: 0}  # 初始化，original=0


    ##################### 调取alpha beta 统一计算 #######################

    # Combine diamond preds and deepgo
    # 加载alpha & beta
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
            alphas[NAMESPACES[ont]] = alpha  # ????????? 方便下面的迭代中使用alpha

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
    print('alphas, beta = ', alphas, betas)




    start_time = time.time()
    total_seq_aa = 0
    total_seq_ss8 = 0


    # 输出结果，每次只算一个ont，在step8_PredictAlpha.sh中设置，开始啦开始啦开始啦！！！！！
    if ont == 'bp':
        out_file_ont = out_file_bp
    elif ont == 'cc':
        out_file_ont = out_file_cc
    elif ont == 'mf':
        out_file_ont = out_file_mf

    # 这里是加载model，深度学习预测
    w = open(out_file_ont, 'wt')


    ################## 深度学习预测_aa ######################

    for prot_ids_aa, sequences_aa in read_fasta(in_file_aa, chunk_size):  # 这里输入的是*.fa文件，test_data.fa。这个read_fasta是yeld返还，不过这个for循环有啥用？
        # prot_ids = ['MMAA3_MYCTU', 'RL14_MYCTU', 'MDH_MYCTU', 'G6PI_MYCTU'...]
        total_seq_aa += len(prot_ids_aa)
        deep_preds_aa = {}
        ids, data = get_data(sequences_aa, PROT_LETTER_len_aa, PROT_INDEX_aa)


        # original: preds = model.predict(data, batch_size=batch_size)  # 啥格式，啥形式？
        # FS change ...
        preds_list_aa = []
        # labels_list_aa = []
        with torch.no_grad():
            for inputs_onehot_aa, target_labels_aa in test_dataloader_aa:
                inputs_onehot_aa = inputs_onehot_aa.to(device_aa)
                # target_labels_aa = target_labels_aa.to(device_aa)

                outputs_aa = model_aa(inputs_onehot_aa)
                preds_list_aa.append(outputs_aa.cpu().numpy())
                # labels_list_aa.append(target_labels_aa.cpu().numpy())
        preds_aa = np.concatenate(preds_list_aa)
        # [[0.00543423 0.009307   0.02373016 ... 0.09436262 0.18778339 0.00478844] [] [] ... ]
        # labels_aa = np.concatenate(labels_list_aa)  # 后面用了没？
        # FS change end ...

        assert preds_aa.shape[1] == len(terms_aa)
        for i, j in enumerate(ids):
            prot_id_aa = prot_ids_aa[j]
            if prot_id_aa not in deep_preds_aa:
                deep_preds_aa[prot_id_aa] = {}
            for l in range(len(terms_aa)):
                if preds_aa[i, l] >= 0.01:
                    if terms_aa[l] not in deep_preds_aa[prot_id_aa]:
                        deep_preds_aa[prot_id_aa][terms_aa[l]] = preds_aa[i, l]
                    else:
                        deep_preds_aa[prot_id_aa][terms_aa[l]] = max(deep_preds_aa[prot_id_aa][terms_aa[l]], preds_aa[i, l])

        # deep_preds =  {'MMAA3_MYCTU': {'GO:0019222': 0.06975121, 'GO:0046872': 0.07646799, ...},
        #                'RL14_MYCTU': {'GO:0019222': 0.06489174, 'GO:0046872': 0.0781393, ...},
        #                MYCTU: {GO, GO},MYCTU: {GO, GO},...}



        ################## 深度学习预测_ss8 ######################  注意：这个和上面的aa是同列对齐！！！，或者缩进一格呢？？？

        for prot_ids_ss8, sequences_ss8 in read_fasta(in_file_ss8, chunk_size):  # 这里输入的是*.fa文件，test_data.fa
            # prot_ids = ['MMAA3_MYCTU', 'RL14_MYCTU', 'MDH_MYCTU', 'G6PI_MYCTU'...]
            total_seq_ss8 += len(prot_ids_ss8)
            deep_preds_ss8 = {}
            ids, data = get_data(sequences_ss8, PROT_LETTER_len_ss8, PROT_INDEX_ss8)

            # original: preds = model.predict(data, batch_size=batch_size)  # 啥格式，啥形式？
            # FS change ...
            preds_list_ss8 = []
            # labels_list_ss8 = []
            with torch.no_grad():
                for inputs_onehot_ss8, target_labels_ss8 in test_dataloader_ss8:
                    inputs_onehot_ss8 = inputs_onehot_ss8.to(device_ss8)
                    # target_labels_ss8 = target_labels_ss8.to(device_ss8)

                    outputs_ss8 = model_ss8(inputs_onehot_ss8)
                    preds_list_ss8.append(outputs_ss8.cpu().numpy())
                    # labels_list_ss8.append(target_labels_ss8.cpu().numpy())
            preds_ss8 = np.concatenate(preds_list_ss8)
            # [[0.00543423 0.009307   0.02373016 ... 0.09436262 0.18778339 0.00478844] [] [] ... ]
            # labels_ss8 = np.concatenate(labels_list_ss8)  # 后面用了没？
            # FS change end ...

            assert preds_ss8.shape[1] == len(terms_ss8)
            for i, j in enumerate(ids):
                prot_id_ss8 = prot_ids_ss8[j]
                if prot_id_ss8 not in deep_preds_ss8:
                    deep_preds_ss8[prot_id_ss8] = {}
                for l in range(len(terms_ss8)):
                    if preds_ss8[i, l] >= 0.01:
                        if terms_ss8[l] not in deep_preds_ss8[prot_id_ss8]:
                            deep_preds_ss8[prot_id_ss8][terms_ss8[l]] = preds_ss8[i, l]
                        else:
                            deep_preds_ss8[prot_id_ss8][terms_ss8[l]] = max(deep_preds_ss8[prot_id_ss8][terms_ss8[l]], preds_ss8[i, l])



            ##### 注意！！！这个是在上面ss8的缩进中，因为用到了prot_ids_ss8.不过因为这个ids是蛋白质名称，所以部分aa or ss8

            ###########################################################
            # FS   alpha & 1-alpha 对调  这是核心！！！！！！！
            ###########################################################
            for prot_id in prot_ids_aa:  # 逐一过 test_data.fa, prot_ids其实只是一系列的protein名字，不包含aa/ss8
                annots = {}  # 这啥？
                if prot_id in diamond_preds:  # blast/diamond * (1-alpha)
                    for go_id, score in diamond_preds[prot_id].items():
                        # print('go.get_namespace(go_id) = ', go.get_namespace(go_id))  # 其实这里还是混着的,3种GO都有
                        # print('alphas[go.get_namespace(go_id)] = ', alphas[go.get_namespace(go_id)])
                        annots[go_id] = score * (1 - alphas[go.get_namespace(go_id)])

                # # original
                # for go_id, score in deep_preds[prot_id].items():  # 深度学习的预测结果
                #     if go_id in annots:
                #         annots[go_id] += alphas[go.get_namespace(go_id)] * score
                #     else:
                #         annots[go_id] = alphas[go.get_namespace(go_id)] * score

                # FS change, for aa
                for go_id, score in deep_preds_aa[prot_id].items():  # 深度学习的预测结果
                    if go_id in annots:
                        annots[go_id] += alphas[go.get_namespace(go_id)] * score
                    else:
                        annots[go_id] = alphas[go.get_namespace(go_id)] * score

                # FS change, for ss8
                for go_id, score in deep_preds_ss8[prot_id].items():  # 深度学习的预测结果
                    if go_id in annots:
                        annots[go_id] += betas[go.get_namespace(go_id)] * score
                    else:
                        annots[go_id] = betas[go.get_namespace(go_id)] * score




                # print('444444444')
                # print(len(annots))
                # annots = {'GO:0000156': 0.07886640414595605, 'GO:0000160': 0.5202934741973877, ...
                # 每一次prot_id循环中的annots长度都不同， 148, 122, ...


                # Propagate scores with ontology structure
                gos = list(annots.keys())
                for go_id in gos:
                    for g_id in go.get_ancestors(go_id):
                        if g_id in annots:
                            annots[g_id] = max(annots[g_id], annots[go_id])
                        else:
                            annots[g_id] = annots[go_id]

                sannots = sorted(annots.items(), key=lambda x: x[1], reverse=True)
                for go_id, score in sannots:
                    if score >= threshold:
                        # 明天修改，把这里添加一个只针对特定ont的
                        # print('before filtering: go_id, ont = ', go_id, go.get_namespace(go_id))

                        # original:
                        # w.write(prot_id + ', ' + go_id + ', ' + go.get_namespace(go_id) + ', ' + go.get_term(go_id)['name'] + ', %.3f\n' % score)
                        # FS: 想想也没必要只挑出来ont的，索性全部写入results.csv里，反正每一行都有注释是bp/mf/cc，需要什么挑什么即可
                        # FS 又想了想，还是要测的，因为这里只读取了 json文件中的某一个ont-的alpha，其他两个默认为0，所以predict的也不准。
                        if go.get_namespace(go_id) == NAMESPACES[ont]:  # i.e. molecular_function  只挑选属于ont(mf)的写入results.csv
                            # print('after filtering: go_id, ont = ', go_id, go.get_namespace(go_id))
                            w.write(prot_id + ', ' + go_id + ', ' + go.get_namespace(go_id) + ', ' + go.get_term(go_id)['name'] + ', %.3f\n' % score)

                w.write('\n')


    w.close()
    total_time = time.time() - start_time
    print('Total prediction time for %d sequences is %d' % (total_seq_aa, total_time))  # total_seq_aa & total_seq_ss8 一样长



############################################################################
######################## class/func 定义各种 函数 & 类 #######################
############################################################################

def to_onehot(seq, PROT_LETTER_len, PROT_INDEX, start=0):  # change seq to onehot，不在PROT_INDEX里面的，就会在第一列标注1
    onehot = np.zeros((params_local_aa['MAXLEN'], PROT_LETTER_len + 1), dtype=np.int32)  # original:(MAXLEN, 21)
    l = min(params_local_aa['MAXLEN'], len(seq))
    for i in range(start, start + l):
        onehot[i, PROT_INDEX.get(seq[i - start], 0)] = 1
    onehot[0:start, 0] = 1
    onehot[start + l:, 0] = 1
    return onehot



# 自定义数据集类，用于批量想pytorch Model1输入数据
class ProteinGODataset(Dataset):
    def __init__(self, df, terms_dict, PROT_LETTER_len, PROT_INDEX):
        self.df = df
        self.terms_dict = terms_dict
        self.PROT_LETTER_len = PROT_LETTER_len
        self.PROT_INDEX = PROT_INDEX

    def __len__(self):
        return len(self.df)

    def __getitem__(self, idx):
        terms_classes = len(self.terms_dict)  # terms_classes = nb_classes (original)
        row = self.df.iloc[idx]
        seq = row.sequences
        onehot = to_onehot(seq, self.PROT_LETTER_len, self.PROT_INDEX, start=0)  # 将蛋白质序列转换为one-hot编码   # onehot:  'numpy.ndarray, (2000,21)
        data_onehot = torch.zeros((params_local_aa['MAXLEN'], self.PROT_LETTER_len + 1), dtype=torch.float32)  # 初始化一个全零张量，用于存储one-hot编码后的序列，original:((MAXLEN, 21)，？？不确定这样引用对不对self.PROT_LETTER_len
        labels = torch.zeros(terms_classes,
                             dtype=torch.float)  # 初始化一个全零张量，用于存储标签  original:dtype=torch.int32 or torch.long
        data_onehot[:onehot.shape[0], :] = torch.from_numpy(
            onehot)  # 将one-hot编码复制到data_onehot张量中，保证长度一致。 onehot.shape[0]=MAXLEN对不？

        for t_id in row.prop_annotations:
            if t_id in self.terms_dict:
                labels[self.terms_dict[t_id]] = 1  # 将标签设为1，表示该样本对应的GO term存在

        return data_onehot, labels  # 返回数据样本及其标签, original 返回 (data_onehot, labels)



#####################################
######## model 模型框架搭建部分 ########
######################################
class Model1(nn.Module):
    def __init__(self, terms_classes, params, kernels_tuple, filters_list):
        super(Model1, self).__init__()

        self.kernels_tuple = kernels_tuple  # 卷积核的大小，(kernel_height, kernel_width)   [(4, 21)], [(128, 21), (128, 21)]
        # print('\n----------- self.kernels_tuple ------------\n')
        # print(self.kernels_tuple)

        self.filters = filters_list  # 定义卷积核的数量列表  [512, 1024, 32768]
        # print('\n-----------self.filters-------------\n')
        # print(self.filters)

        self.convs = nn.ModuleList()
        self.pools = nn.ModuleList()
        self.flats = nn.ModuleList()

        # 创建卷积层、池化层和展平层的列表
        for i in range(len(self.kernels_tuple)):
            '''
            对于灰度图像来说，通道数为 1，因为它只有一个颜色通道。对于彩色图像来说，通道数通常为 3，因为它包含红、绿、蓝三个通道
            在给定的代码中，卷积层 nn.Conv2d 的输入通道数 in_channels 被设置为 1。
            这是因为 Model1 类的输入数据被表示为一个单通道的输入张量，其中每个输入样本是一个形状为 (2000, 21) 的二维矩阵，表示一个灰度图像。
            '''
            conv = nn.Conv2d(
                in_channels=1,  # 输入通道数为1，可以理解为单层的灰度图。如果是RGP3层，则in_channel=3
                # 接上面：正常3层图像输入是(3, 255, 255)，本例中onehot得到输入模型的初始参数是(1024, 21)。也就是summary里的
                out_channels=self.filters[i],  # 从2维降到1维。(4, 21) 沿着 输入的(32, 2000, 21)中2000的维度扫。32=batch_size
                kernel_size=self.kernels_tuple[i],
                padding=(0,),
                stride=(1, 1)
            )

            # 使用self.kernels_tuple[i][0]获取卷积核大小，用于池化核，(pool_height, pool_width)，这个可能没必要用2D?
            pool = nn.MaxPool2d(kernel_size=(params_local_aa['MAXLEN'] - self.kernels_tuple[i][0] + 1, 1))
            # print(f"MaxPool{i + 1} kernel size:", pool.kernel_size)

            flat = nn.Flatten()  # 展平
            self.convs.append(conv)  # 排排队，吃果果
            self.pools.append(pool)
            self.flats.append(flat)

        # 创建 全连接层（FC）& 输出层（output）
        # FS change, 假定 FC 有多层，
        # 第一层 input = sum(self.filters) (--上述flat展平全长)，output = terms_classes
        # 后面所有层： input = output = terms_classes
        if params_local_aa['FC_depth'] > 0:  # 存在全连接层
            self.fc_layers = nn.ModuleList()  # fc=full connected全连接层
            for i in range(params_local_aa['FC_depth']):
                if i == 0:  # 第一层   # input = sum(self.filters)，output = terms_classes
                    fc = nn.Linear(sum(self.filters), terms_classes)
                else:
                    fc = nn.Linear(terms_classes, terms_classes)  # 后面所有层： input = output = terms_classes 输入输出一样
                self.fc_layers.append(fc)
            self.output_layer = nn.Linear(terms_classes, terms_classes)  # 存在全连接层时候的输出层

        else:  # 不存在全连接层，FC_depth = 0
            self.output_layer = nn.Linear(sum(self.filters), terms_classes)  # 输出层，直接从 ‘flat展平全长’ 蹦到 terms_classes

    def forward(self, x):
        nets = []
        # net = torch.empty(32,x)
        for i in range(len(self.kernels_tuple)):
            # x.unsqueeze(1) 是用来在输入张量 x 的维度上增加一个维度。作用是在索引为 1 的位置上增加一个维度，也就是在通道维度上增加一个维度
            # x.shape=(32, 2000, 21)，这样处理后, shape变成(32, 1, 2000, 21)
            conv = self.convs[i](x.unsqueeze(1))
            pool = self.pools[i](conv)
            flat = self.flats[i](pool)
            nets.append(flat)  # 这只是一个list，如果有多个元素，需要 torch.cat链接成一个完整的 张量 层

        if len(self.kernels_tuple) == 1:  # 如果 kernel只有一个
            net = nets[0]
        elif len(self.kernels_tuple) > 1:  # 如果 kernel的list有多个元素，则需要并排
            net = torch.cat(nets, dim=1)
        else:
            print('error ................')

        # print("net shape:", net.shape)  # 打印net的形状
        # print("nets shape:", [n.shape for n in nets])  # 打印nets中每个元素的形状
        # print(net)
        # print(type(net))

        if params_local_aa['FC_depth'] > 0:  # 如果存在全连接层，则把fc_layers 加在net上，
            for k in range(len(self.fc_layers)):
                net = self.fc_layers[k](net)  #
            # print("net shape:", net.shape)  # 打印net的形状

        net = self.output_layer(net)  # 链接到输出层
        net = torch.sigmoid(net)

        return net


#######

def get_data(sequences, PROT_LETTER_len, PROT_INDEX):
    pred_seqs = []
    ids = []
    for i, seq in enumerate(sequences):
        if len(seq) > MAXLEN:
            st = 0
            while st < len(seq):
                pred_seqs.append(seq[st: st + MAXLEN])
                ids.append(i)
                st += MAXLEN - 128
        else:
            pred_seqs.append(seq)
            ids.append(i)
    n = len(pred_seqs)
    data = np.zeros((n, MAXLEN, PROT_LETTER_len + 1), dtype=np.float32)

    for i in range(n):
        seq = pred_seqs[i]
        data[i, :, :] = to_onehot(seq, PROT_LETTER_len, PROT_INDEX)
    return ids, data


def read_fasta(filename, chunk_size):  # FS
    seqs = list()
    info = list()
    seq = ''
    inf = ''
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq != '':
                    seqs.append(seq)
                    info.append(inf)
                    if len(info) == chunk_size:
                        yield (info, seqs)
                        seqs = list()
                        info = list()
                    seq = ''
                inf = line[1:].split()[0]
            else:
                seq += line
        seqs.append(seq)
        info.append(inf)
    yield (info, seqs)


if __name__ == '__main__':
    main()
