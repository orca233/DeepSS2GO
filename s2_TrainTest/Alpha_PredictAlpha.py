#!/usr/bin/env python

'''
引自deepgoplus中的predict.py
很多函数引自 step2_Test.py

需要model_checkpoint.pth，对test_data.fa/pkl进行预测
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
from step0_TrainTestSetting_local import *

MAXLEN = params_local['MAXLEN']

@ck.command()
@ck.option('--in-file', '-if', default='data/test_data.fa', help='Input FASTA file', required=True)
@ck.option('--out-file', '-of', default='data/results.csv', help='Output result file')
@ck.option('--go-file', '-gf', default=params_local['path_base'] + 'pub_data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
@ck.option('--model-file', '-mf', default='data/model_checkpoint.pth', help='Tensorflow model file')
@ck.option('--terms-file', '-tf', default='data/terms_gominre_trxte.pkl', help='List of predicted terms')
@ck.option('--annotations-file', '-tf', default='data/train_data.pkl', help='Experimental annotations')
@ck.option('--chunk-size', '-cs', default=100000, help='Number of sequences to read at a time')  # original 1000
@ck.option('--diamond-file', '-df', default='data/diamond_aa.res', help='Diamond Mapping file')
@ck.option('--threshold', '-t', default=0.1, help='Prediction threshold')
@ck.option('--batch-size', '-bs', default=32, help='Batch size for prediction model')
# @ck.option('--maxlen', '-ml', default=params_local['MAXLEN'])
# new
@ck.option('--test-data-file', '-tedf', default='data/test_data.pkl', help='XX')
@ck.option('--aa-ss', '-aass', default=params_local['aa_ss'], help='aa/ss8/ss3')
@ck.option('--prot-letter-aa', '-plaa', default=params_local['PROT_LETTER_aa'], type=list, help='XX')
@ck.option('--prot-letter-ss8', '-plss8', default=params_local['PROT_LETTER_ss8'], type=list, help='XX')
@ck.option('--prot-letter-ss3', '-plss3', default=params_local['PROT_LETTER_ss3'], type=list, help='XX')
@ck.option('--kernels-list', '-kl', default=params_local['kernels'], type=list, help='XX')      # # kernels & filters 需要是list形式，即使只有一个 int，也要转成list
@ck.option('--filters-list', '-fl', default=params_local['filters'], type=list, help='XX')

# @ck.option('--maxlen', '-fl', default=params_local['filters'], type=list, help='XX')
# @ck.option('--filters-list', '-fl', default=params_local['filters'], type=list, help='XX')


def main(in_file, out_file, go_file, model_file, terms_file, annotations_file, chunk_size, diamond_file,
         threshold, batch_size,
         test_data_file, aa_ss, prot_letter_aa, prot_letter_ss8, prot_letter_ss3,
         kernels_list, filters_list):  # maxlen


    #######################################################################
    ##################### pre_setting 数据和参数的预处理 ######################
    ########################################################################

    ##### 根据 params_local的aa/ss8/ss3确定PROT_LETTER，kernels_tuple #####
    PROT_LETTER = []
    PROT_LETTER_len = -1
    if aa_ss == 'aa':
        PROT_LETTER = prot_letter_aa
        PROT_LETTER_len = len(PROT_LETTER)  # =20
    elif aa_ss == 'ss8':
        PROT_LETTER = prot_letter_ss8
        PROT_LETTER_len = len(PROT_LETTER)  # =8
    elif aa_ss == 'ss3':
        PROT_LETTER = prot_letter_ss3
        PROT_LETTER_len = len(PROT_LETTER)  # =3

    PROT_INDEX = dict()
    for i in range(len(PROT_LETTER)):
        PROT_INDEX[PROT_LETTER[i]] = i + 1


    kernels_tuple = [(value, PROT_LETTER_len + 1) for value in kernels_list]

    go = Ontology(go_file, with_rels=True)
    # test_data_file = 'data/test_data.pkl'
    terms_df = pd.read_pickle(terms_file)
    terms = terms_df['terms'].values.flatten()
    terms_dict = {v: i for i, v in enumerate(terms)}
    terms_classes = len(terms_dict)


    print('11111111111')

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

    # 创建ProteinDataset实例，不同于train，这里只有test的实例
    test_df = pd.read_pickle(test_data_file)  # 这里需要test_data.pkl文件，重新进行类似于step3_Test.py的步骤，根据训练好的model再预测
    test_dataset = ProteinGODataset(test_df, terms_dict, PROT_LETTER_len, PROT_INDEX)
    test_dataloader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)


    # 判断用 single_GPU / multi_GPU 并行运算。      创建model模型实例
    device_ids = params_local['device_ids']  # cuda:3"   or    [0, 1, 2, 3]

    if isinstance(device_ids, str):  # device_ids是个str="cuda:3"，用单GPU
        print('--- single GPU = ', device_ids)
        device = torch.device(device_ids if torch.cuda.is_available() else "cpu")  # 将模型移动到cuda:2
        model = Model1(terms_classes, params_local)  # 创建Model1模型实例
        model.to(device)  # 将模型移动到GPU（单个或多个）

    else:  # device_ids是个list=[0, 1, 2, 3]，用多GPU，使用DataParallel进行并行计算
        print('--- multi GPU = ', device_ids)
        # device = torch.device("cuda" if torch.cuda.is_available() else "cpu")  # 对于device_ids = [2, 1]这种有问题，猜的啊：可能是吧model.to(device)默认扔给cuda:0了
        device = torch.device("cuda:{}".format(device_ids[0]) if torch.cuda.is_available() else "cpu")  # 将device_ids[0]的值动态地插入到字符串中，生成一个完整的字符串，例如"cuda:2"

        model = Model1(terms_classes, params_local, kernels_tuple, filters_list)  # 创建Model1模型实例
        model.to(device)  # 将模型移动到GPU（单个或多个）
        model = nn.DataParallel(model, device_ids=device_ids)  # 将模型包装在DataParallel中，实现多GPU的并行计算
        model.to(device)  # 将包装后的模型再次移动到指定设备，以确保模型及其参数和缓冲区都位于正确的设备上


    # 加载模型 ????????????????
    model.load_state_dict(torch.load(model_file))
    model.to(device)
    model.eval()


    '''
    光加载了，还没预测呢？？？？？？？？？？？？？？？？？？？？？？
    -后面会用这个model进行预测！！！！！！！！！！！！！！
    '''




    #################################################
    ########## 第三部分：整合 #########################
    #################################################


    # alphas = {NAMESPACES['mf']: 0.48, NAMESPACES['bp']: 0.4, NAMESPACES['cc']: 0.49}
    alphas = {NAMESPACES['mf']: 1, NAMESPACES['bp']: 1, NAMESPACES['cc']: 1}

    start_time = time.time()
    total_seq = 0


    # 输出结果
    w = open(out_file, 'wt')
    for prot_ids, sequences in read_fasta(in_file, chunk_size):  # 这里输入的是*.fa文件，test_data.fa
        total_seq += len(prot_ids)
        deep_preds = {}
        ids, data = get_data(sequences, PROT_LETTER_len, PROT_INDEX)

        # original: preds = model.predict(data, batch_size=batch_size)  # 啥格式，啥形式？
        # FS change ...
        # 下面这一段是吧model带进去用来预测
        preds_list = []
        labels_list = []
        with torch.no_grad():
            for inputs_onehot, target_labels in test_dataloader:
                inputs_onehot = inputs_onehot.to(device)
                target_labels = target_labels.to(device)

                outputs = model(inputs_onehot)
                preds_list.append(outputs.cpu().numpy())
                labels_list.append(target_labels.cpu().numpy())
        preds = np.concatenate(preds_list)
        # [[0.00543423 0.009307   0.02373016 ... 0.09436262 0.18778339 0.00478844] [] [] ... ]
        labels = np.concatenate(labels_list)  # 后面用了没？
        # FS change end ...

        assert preds.shape[1] == len(terms)
        for i, j in enumerate(ids):
            prot_id = prot_ids[j]
            if prot_id not in deep_preds:
                deep_preds[prot_id] = {}
            for l in range(len(terms)):
                if preds[i, l] >= 0.01:
                    if terms[l] not in deep_preds[prot_id]:
                        deep_preds[prot_id][terms[l]] = preds[i, l]
                    else:
                        deep_preds[prot_id][terms[l]] = max(deep_preds[prot_id][terms[l]], preds[i, l])

        # Combine diamond preds and deepgo
        '''
        # original
        for prot_id in prot_ids:
            annots = {}
            if prot_id in diamond_preds:  # blast/diamond * (1-alpha)
                for go_id, score in diamond_preds[prot_id].items():
                    annots[go_id] = score * alphas[go.get_namespace(go_id)]
            for go_id, score in deep_preds[prot_id].items():
                if go_id in annots:
                    annots[go_id] += (1 - alphas[go.get_namespace(go_id)]) * score
                else:
                    annots[go_id] = (1 - alphas[go.get_namespace(go_id)]) * score
        '''

        # FS   alpha & 1-alpha 对调
        for prot_id in prot_ids:
            annots = {}
            if prot_id in diamond_preds:  # blast/diamond * (1-alpha)
                for go_id, score in diamond_preds[prot_id].items():
                    annots[go_id] = score * (1 - alphas[go.get_namespace(go_id)])
            for go_id, score in deep_preds[prot_id].items():
                if go_id in annots:
                    annots[go_id] += alphas[go.get_namespace(go_id)] * score
                else:
                    annots[go_id] = alphas[go.get_namespace(go_id)] * score



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
                    w.write(prot_id + ', ' + go_id + ', ' + go.get_namespace(go_id) + ', ' +
                            go.get_term(go_id)['name'] + ', %.3f\n' % score)

            w.write('\n')
    w.close()
    total_time = time.time() - start_time
    print('Total prediction time for %d sequences is %d' % (total_seq, total_time))







############################################################################
######################## class/func 定义各种 函数 & 类 #######################
############################################################################

def to_onehot(seq, PROT_LETTER_len, PROT_INDEX, start=0):  # change seq to onehot，不在PROT_INDEX里面的，就会在第一列标注1
    onehot = np.zeros((params_local['MAXLEN'], PROT_LETTER_len + 1), dtype=np.int32)  # original:(MAXLEN, 21)
    l = min(params_local['MAXLEN'], len(seq))
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
        data_onehot = torch.zeros((params_local['MAXLEN'], self.PROT_LETTER_len + 1), dtype=torch.float32)  # 初始化一个全零张量，用于存储one-hot编码后的序列，original:((MAXLEN, 21)，？？不确定这样引用对不对self.PROT_LETTER_len
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
            pool = nn.MaxPool2d(kernel_size=(params_local['MAXLEN'] - self.kernels_tuple[i][0] + 1, 1))
            # print(f"MaxPool{i + 1} kernel size:", pool.kernel_size)

            flat = nn.Flatten()  # 展平
            self.convs.append(conv)  # 排排队，吃果果
            self.pools.append(pool)
            self.flats.append(flat)

        # 创建 全连接层（FC）& 输出层（output）
        # FS change, 假定 FC 有多层，
        # 第一层 input = sum(self.filters) (--上述flat展平全长)，output = terms_classes
        # 后面所有层： input = output = terms_classes
        if params_local['FC_depth'] > 0:  # 存在全连接层
            self.fc_layers = nn.ModuleList()  # fc=full connected全连接层
            for i in range(params_local['FC_depth']):
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

        if params_local['FC_depth'] > 0:  # 如果存在全连接层，则把fc_layers 加在net上，
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
