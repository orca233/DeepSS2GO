# README，本篇一共包含4个部分，
# 函数&类定义，包括：普通函数，模型Model1搭建，训练验证函数，测试函数
# pre_setting 数据和参数的预处理
# 训练部分，
# 测试部分


'''
改编自deepgoplus.py
'''


import click as ck
import torch
import torch.nn as nn
import torch.optim as optim
# import torch.nn.functional as F  # 这啥？

import pandas as pd
import numpy as np
import math
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import roc_curve, auc
# from sklearn.model_selection import train_test_split
from utils import EarlyStopping

from step0_TrainTestSetting_local import *
from step0_TrainTestSetting_global import path_base


# from torchsummary import summary
from tqdm import tqdm
import os
from torch.nn.functional import relu  # Import the ReLU function

print('\n\n ----- starting step2 --------------')


@ck.command()
# @ck.option('--go-file', '-gf', default='../../pub_data/go.obo', help='Gene Ontology file in OBO Format')  # ../../pub_data/go.obo
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
# @ck.option('--go-file', '-gf', default=path_base + 'pub_data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='train data file')
@ck.option('--test-data-file', '-tedf', default='data/test_data.pkl', help='XX')
@ck.option('--aa-ss', '-aass', default=params_local['aa_ss'], help='aa/ss8/ss3')
@ck.option('--path-base', '-pb', default=params_local['path_base'], help='..')

# new:
@ck.option('--model-file', '-mf', default='data/model_checkpoint.pth', help='XX')
@ck.option('--predictions-file', '-pf', default='data/predictions.pkl', help='XX')
@ck.option('--logger-file', '-lf', default='data/training.csv', help='XX')  # 每一步epoch保存一个loss
@ck.option('--terms-file', '-tf', default='data/terms_gominre_trxte.pkl', help='XX')  # terms_file = 'data/terms_%s.pkl' % params_local['onts']  # 給terms_file 指定对应的terms_bp/cc/mf_file文件夹
@ck.option('--kernels-list', '-kl', default=params_local['kernels'], type=list, help='XX')      # # kernels & filters 需要是list形式，即使只有一个 int，也要转成list
@ck.option('--filters-list', '-fl', default=params_local['filters'], type=list, help='XX')
@ck.option('--prot-letter-aa', '-plaa', default=params_local['PROT_LETTER_aa'], type=list, help='XX')
@ck.option('--prot-letter-ss8', '-plss8', default=params_local['PROT_LETTER_ss8'], type=list, help='XX')
@ck.option('--prot-letter-ss3', '-plss3', default=params_local['PROT_LETTER_ss3'], type=list, help='XX')

@ck.option('--run-train-switch', '-rtrs', default='F', help='XX')  # step2中默认加载 run_train, 不加载 run_test。 step3 中反过来
@ck.option('--run-test-switch', '-rtes', default='T', help='XX')  # step2中默认加载 run_train, 不加载 run_test。 step3 中反过来


def main(go_file, train_data_file, test_data_file, aa_ss, path_base, model_file, predictions_file,
         kernels_list, filters_list, logger_file, terms_file, prot_letter_aa, prot_letter_ss8, prot_letter_ss3,
         run_train_switch, run_test_switch):
    print('\n################## a long, long time ago ... ##################\n')

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


    # # # kernels & filters 需要是list形式，即使只有一个 int，也要转成list
    # kernels_list = params_local['kernels']
    # filters_list = params_local['filters']

    # 对aa，PROT_LETTER_len=21，[16, 24]-->[(16, 21), (24, 21)], 对ss8, 得到-->[(16, 9), (24, 9)]
    kernels_tuple = [(value, PROT_LETTER_len + 1) for value in kernels_list]

    # go = Ontology(go_file, with_rels=True)  # 训练的时候没用到GO.obo????

    terms_df = pd.read_pickle(terms_file)  # terms_file = 'data/terms_mf.pkl'  或者  'data/terms_all.pkl'
    terms = terms_df['terms'].values.flatten()  # 2088,numpy.ndarray, ['GO:0006796' 'GO:0051017' 'GO:0043065' ... ] # .flatten() 方法将该数组压缩为一维数组
    terms_dict = {v: i for i, v in enumerate(terms)}  # index, value反过来做成dict.
    # terms_dict = {'GO:0006796': 0, 'GO:0051017': 1, 'GO:0043065': 2, 'GO:0002790': 3, ... 'GO:0051707': 2087}
    terms_classes = len(terms_dict)


    #####################################################################################
    ########################### 开始训练 train + valid ####################################
    #####################################################################################

    ################# 数据导入 ##############
    # 根据“TrainValidRatio”分割train,valid数据
    # train : valid : test = 11318 : 1258 : 662 = 13238, (train+valid)/test = 0.95, train/valid = 0.9
    train_df, valid_df = split_train_valid(train_data_file, terms, params_local['TrainValidRatio'])
    test_df = pd.read_pickle(test_data_file)

    # 创建数据集和数据加载器, 创建ProteinDataset的实例
    train_dataset = ProteinGODataset(train_df, terms_dict, PROT_LETTER_len, PROT_INDEX)
    valid_dataset = ProteinGODataset(valid_df, terms_dict, PROT_LETTER_len, PROT_INDEX)
    test_dataset = ProteinGODataset(test_df, terms_dict, PROT_LETTER_len, PROT_INDEX)

    # 创建一个DataLoader来以指定的批次大小加载数据
    train_dataloader = DataLoader(train_dataset, batch_size=params_local['batch_size'], shuffle=True)
    valid_dataloader = DataLoader(valid_dataset, batch_size=params_local['batch_size'], shuffle=False)
    test_dataloader = DataLoader(test_dataset, batch_size=params_local['batch_size'], shuffle=False)  # 注意test_dataloader有用到batch，不知道是 批次算结果？还是逐一算？

    # 判断用 single_GPU / multi_GPU 并行运算。      创建model模型实例
    device_ids = params_local['device_ids']  # cuda:3"   or    [0, 1, 2, 3]

    if isinstance(device_ids, str):  # device_ids是个str="cuda:3"，用单GPU
        print('--- single GPU = ', device_ids)
        device = torch.device(device_ids if torch.cuda.is_available() else "cpu")  # 将模型移动到cuda:2
        model = Model1(terms_classes, params_local, kernels_tuple, filters_list)  # 创建Model1模型实例
        model.to(device)  # 将模型移动到GPU（单个或多个）

    else:  # device_ids是个list=[0, 1, 2, 3]，用多GPU，使用DataParallel进行并行计算
        print('--- multi GPU = ', device_ids)
        # device = torch.device("cuda" if torch.cuda.is_available() else "cpu")  # 对于device_ids = [2, 1]这种有问题，猜的啊：可能是吧model.to(device)默认扔给cuda:0了
        device = torch.device("cuda:{}".format(device_ids[0]) if torch.cuda.is_available() else "cpu")  # 将device_ids[0]的值动态地插入到字符串中，生成一个完整的字符串，例如"cuda:2"

        model = Model1(terms_classes, params_local, kernels_tuple, filters_list)  # 创建Model1模型实例
        model.to(device)  # 将模型移动到GPU（单个或多个）
        model = nn.DataParallel(model, device_ids=device_ids)  # 将模型包装在DataParallel中，实现多GPU的并行计算
        model.to(device)  # 将包装后的模型再次移动到指定设备，以确保模型及其参数和缓冲区都位于正确的设备上

    print('---------------- This is model: ----------------')
    print('model shape by general model ::::')
    print(model)

    '''
    print('----- model shape by summary :::: -----')  # summary在单个GPU只能用于cuda:0。所以放在multiGPU里可通用
    summary(
        model,
        input_size=(params_local['MAXLEN'], PROT_LETTER_len + 1),  # (1024, 21)
        batch_size=params_local['batch_size']  # 32
    )
    '''

    # 使用训练集和验证集训练模型
    optimizer = optim.Adam(model.parameters(), lr=params_local['learning_rate'])  # lr=3e-4
    criterion = nn.BCELoss()

    ######### step3_Test.py 中，把下面两行注释掉 ###########
    if run_train_switch == 'T':  # step2中默认加载 run_train, 不加载 run_test。 step3 中反过来
        run_train(model, train_dataloader, valid_dataloader, device, optimizer, criterion, logger_file)  # 可以控制是否跳过train阶段

    ################ 训练过程结束 ###################


    #####################################################################################
    ########################### 开始测试 test ############################################
    #####################################################################################

    ######### step2_Train.py 中，把下面七八行注释掉 ###########
    if params_local['load_pretrained_model'] == 1:  # 如果需要加载‘已经训练好的模型’，则调用下面的地址
        model_file = params_local['load_pretrained_model_addr']

    model.load_state_dict(torch.load(model_file))
    model.to(device)

    # 测试模型并保存预测结果
    if run_test_switch == 'T':  # step2中默认加载 run_train, 不加载 run_test。 step3 中反过来
        run_test(model, test_dataloader, device, criterion, test_df, predictions_file)



    print('\n################## And they all lived happily ever after! ##################\n')
    ### __main__ 结束





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

# 把 train_data_file 分成 train & valid，
# 注意！！！ 对于TrainHUMAN_TestECOLI，也是train&valid 同源于 HUMAN，只是在step3验证时候才用ECOLI
def split_train_valid(data_file, terms, trainvalidratio):  # trainvalidratio data into train & test/valid session
    df = pd.read_pickle(data_file)
    # split train/valid with trainvalidratio
    n = len(df)
    index = np.arange(n)
    train_n = int(n * trainvalidratio)
    np.random.seed(seed=0)
    np.random.shuffle(index)  # 洗牌
    train_df = df.iloc[index[:train_n]]
    valid_df = df.iloc[index[train_n:]]
    return train_df, valid_df

def compute_roc(labels, preds):  #####计算ROC函数  比较真实值labels与模型预测值predictions的差别
    # Compute ROC curve and ROC area for each class
    fpr, tpr, _ = roc_curve(labels.flatten(), preds.flatten())
    roc_auc = auc(fpr, tpr)
    return roc_auc

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
    def __init__(self, terms_classes, params, kernels_tuple, filters_list):  # terms_classes, params, kernels_tuple, filters_list
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

        # if params_local['FC_depth'] > 0:  # 如果存在全连接层，则把fc_layers 加在net上，
        #     for k in range(len(self.fc_layers)):
        #         net = self.fc_layers[k](net)  #
        #     # print("net shape:", net.shape)  # 打印net的形状
        #
        # net = self.output_layer(net)  # 链接到输出层
        # net = torch.sigmoid(net)

        # FS添加 19/12/2023，flat结束之后，添加relu
        if params_local['FC_depth'] > 0:  # 如果存在全连接层，则把fc_layers 加在net上，
            for k in range(len(self.fc_layers)):
                net = self.fc_layers[k](net)
                net = relu(net)  # Apply the ReLU function after each fully connected layer
            # print("net shape:", net.shape)  # 打印net的形状

        net = self.output_layer(net)  # 链接到输出层
        net = torch.sigmoid(net)

        return net

##########################################################################
####################### train + valid 训练+验证 函数定义 ####################
##########################################################################

def run_train(model, train_dataloader, valid_dataloader, device, optimizer, criterion, logger_file):
    ################## train 训练部分 ###############
    # 定义Early Stopping
    early_stopping = EarlyStopping(patience=params_local['EarlyStopping_patience'], verbose=True)  # original, patience=6

    # 创建一个DataFrame用于保存训练日志
    logger = pd.DataFrame(columns=['epoch', 'train_loss', 'valid_loss'])  # 还想要：'train_accuracy', 'test_accuracy'？？？

    for epoch in range(params_local['epochs']):
        model.train()  # 将模型设置为训练模式
        train_loss = 0

        # 使用tqdm函数包装train_dataloader，并设置了描述字符串为"Epoch {epoch+1}/{epochs}"。
        # 在每个batch的训练过程中，通过set_postfix方法更新进度条上显示的信息，这里选择显示训练损失。在每个epoch结束时，通过close方法关闭进度条。
        train_dataloader_with_progress = tqdm(train_dataloader, desc=f"Epoch {epoch + 1}/{params_local['epochs']}", leave=False)

        for inputs_onehot, target_labels in train_dataloader_with_progress:  # original: for ixx, txx in train_dataloader:
            inputs_onehot = inputs_onehot.to(device)
            target_labels = target_labels.to(device)
            # print('\n1111111111111111111111 train, input_onehot, target_labels')
            # print(inputs_onehot.shape, target_labels.shape)  # torch.Size([32, 1024, 21]) torch.Size([32, 2088])

            optimizer.zero_grad()  # 将优化器中的梯度缓存置零，准备进行反向传播计算新的梯度
            outputs = model(inputs_onehot)  # 通过向模型输入输入数据，获取模型的输出

            loss = criterion(outputs, target_labels)  # 使用指定的损失函数（criterion）计算模型输出与目标标签之间的损失
            loss.backward()  # 执行反向传播，计算损失相对于模型参数的梯度
            optimizer.step()  # 根据计算得到的梯度更新模型参数，执行优化器的参数更新步骤

            train_loss += loss.item()  # 累加每个批次的损失，用于计算平均训练损失

            # 更新tqdm进度条
            train_dataloader_with_progress.set_postfix({'Loss': torch.tensor(train_loss).item()})

        train_loss /= len(train_dataloader)  # # 计算平均训练损失，将总损失除以批次的数量。

        # 在每个epoch结束时，更新tqdm进度条
        train_dataloader_with_progress.close()

        ################## valid 验证部分 ###############
        model.eval()
        valid_loss = 0.0

        # 使用tqdm函数包装valid_dataloader，并设置描述字符串为"Validation"
        valid_dataloader_with_progress = tqdm(valid_dataloader, desc="Validation", leave=False)

        with torch.no_grad():  # 这啥
            for inputs_onehot, target_labels in valid_dataloader_with_progress:  # for ixx, txx in valid_dataloader:
                inputs_onehot = inputs_onehot.to(device)
                target_labels = target_labels.to(device)

                # print('\n2222222222222222 valid, input_onehot, target_labels')
                # print(inputs_onehot.shape, target_labels.shape)  # torch.Size([32, 1024, 21]) torch.Size([32, 2088])


                outputs = model(inputs_onehot)
                valid_loss += criterion(outputs, target_labels).item()  # .item()???

            valid_loss /= len(valid_dataloader)  # 要不要前进一格？

            # 更新tqdm进度条
            valid_dataloader_with_progress.set_postfix({'Loss': torch.tensor(valid_loss).item()})

            # 关闭tqdm进度条
            valid_dataloader_with_progress.close()

            logger = logger.append({
                'epoch': epoch + 1,
                'train_loss': train_loss,
                'valid_loss': valid_loss
            }, ignore_index=True)

            print(f"Epoch {epoch + 1}/{params_local['epochs']} --------------- train_loss: {train_loss:.4f} - Valid_loss: {valid_loss:.4f}")


            # 更新Early Stopping状态并检查是否需要提前停止训练
            if early_stopping(valid_loss, model):
                print("Early Stopping!")
                break

        # 将训练日志保存到文件中
        logger.to_csv(logger_file, index=False)

        if early_stopping.early_stop:
            print("Early Stopping!")
            break


##########################################################
####################### test 函数定义 #####################
###########################################################
# 修改测试部分代码  GPT GPT GPT
def run_test(model, test_dataloader, device, criterion, test_df, predictions_file):
    model.eval()
    test_loss = 0  # 用来和 train_loss, val_loss 对比
    preds_list = []
    labels_list = []

    with torch.no_grad():
        for inputs_onehot, target_labels in test_dataloader:  # 注意test_dataloader有用到batch，不知道是 批次算结果？还是逐一算？
            inputs_onehot = inputs_onehot.to(device)
            target_labels = target_labels.to(device)

            # print('\n333333333333333333 test, input_onehot, target_labels')
            # print(inputs_onehot.shape, target_labels.shape)  # torch.Size([32, 1024, 21]) torch.Size([32, 2088])


            outputs = model(inputs_onehot)
            test_loss += criterion(outputs, target_labels).item()

            preds_list.append(outputs.cpu().numpy())
            labels_list.append(target_labels.cpu().numpy())

    test_loss /= len(test_dataloader)
    print('\n-------------test_loss = ', test_loss)   # 用来和 train_loss, val_loss 对比

    preds = np.concatenate(preds_list)  # numpy.ndarray (662, 2088)
    # [[0.00543423 0.009307   0.02373016 ... 0.09436262 0.18778339 0.00478844] [] [] ... ]
    labels = np.concatenate(labels_list)

    roc_auc = compute_roc(labels, preds)
    print('\n------------  ROC AUC = ', roc_auc)

    # 保存预测结果
    test_df['labels'] = list(labels)   # Length of values (11318) does not match length of index (662) ????????????  Test时候就不要batch了！！！！！！！！！！！！
    test_df['preds'] = list(preds)
    test_df.to_pickle(predictions_file)


# 运行文件
if __name__ == '__main__':
    main()





'''
filters 尺寸的设置来自于，多大的决定蛋白质功能的 氨基酸链的长度。
如果都是100左右的蛋白质决定功能，则将128这个kernel的filter个数设置的大一些。
并不一定是小size的kernel就更能决定GO功能
这里就可做一个交叉实验，尝试不同的 filter尺寸=8, 16, ... 1024，来判断哪个长度的kernel更能决定GO功能

categorical_crossentropy 和 binary_crossentropy 是用于训练神经网络中的损失函数。二者的主要区别在于输出层的形式。
binary_crossentropy 用于二元分类问题，输出层只有一个神经元，其输出是一个介于0和1之间的概率，表示正样本的概率。
categorical_crossentropy 用于多元分类问题，输出层有多个神经元，每个神经元的输出都是一个介于0和1之间的概率，表示属于每个类别的概率。
因此，binary_crossentropy 适用于仅有两个类别的问题，例如正面/负面情感分析，垃圾邮件分类等；
而categorical_crossentropy 适用于有多个类别的问题，例如图像分类、自然语言处理中的文本分类等。

需要注意的是，如果是多元分类问题且标签是 one-hot 编码的，那么可以使用 binary_crossentropy 作为损失函数，!!!!!
但是需要将输出层神经元个数设置为类别数目。!!!!!!!!!!!!
'''


'''
# 检查GPU是否可用
if torch.cuda.is_available():
    device = torch.device("cuda:0")  # 选择第一块GPU, keras中用 gpu:3
    torch.cuda.set_device(device)
    print("已选择GPU:", torch.cuda.get_device_name(device))
else:
    device = torch.device("cpu")
    print("没有可用的GPU，使用CPU.")

torch_device = torch.device(device)  # device = "cuda:0"   选择第一块GPU
torch.cuda.set_device(torch_device)
'''


'''    
keras，'sequence' 是把batch_size带进来，然后划分好所有的数据，idx是每一个batch的代号。还要手动写 next来加载下一个
每一个batch中的data_onehot是三维的，[idx, MAXLEN, 21]，labels也二维的，[idx, 0-1-list].
data_onehot的后两维和labels的第二维度是对应的。

但是！！！ pytorch框架中，'Dataset'负责读取‘每一batch’数据，‘DataLoader’引入batch_size，
然后根据batch_size读取前面'Dataset'中相应行的数据进行train。而后自动加载下一个。

首先定义了一个名为ProteinDataset的自定义数据集类，继承自torch.utils.data.Dataset。在__getitem__方法中，根据给定的索引idx获取相应的数据样本，并将蛋白质序列转换为对应的one-hot编码表示。
然后，根据转换后的数据创建了一个ProteinDataset实例，接着使用torch.utils.data.DataLoader来创建一个数据加载器dataloader，用于以指定的批次大小加载数据。
请注意，需要根据你的需求实现to_onehot函数，以将蛋白质序列转换为one-hot编码表示。此函数应该返回一个numpy数组或torch张量。
'''




'''
使用summary打印model参数时候的一些注意事项：

kernel_size=(16, 21), filter = 128, batchsize=32,时

### 单核模式：
！！！单核模式下，必须用cuda:0，才能使用summary，否则出错。

single GPU =  cuda:0

----------------------------------------------------------------
        Layer (type)               Output Shape         Param #
================================================================
            Conv2d-1         [32, 128, 1009, 1]          43,136
         MaxPool2d-2            [32, 128, 1, 1]               0
           Flatten-3                  [32, 128]               0
            Linear-4                 [32, 1517]         195,693
================================================================
Total params: 238,829
Trainable params: 238,829
Non-trainable params: 0
----------------------------------------------------------------
Input size (MB): 2.62
Forward/backward pass size (MB): 31.96
Params size (MB): 0.91
Estimated Total Size (MB): 35.50
----------------------------------------------------------------



### 多核模式：
！！！不论用哪个核，如果使用summary，则cuda:0都会被调用。而且从始至终都在调用，why？？？？？
！！！不用summary，cuda:0就不会被特别调用了！！！


multi GPU = [2]，多核模式用单核：
----------------------------------------------------------------
        Layer (type)               Output Shape         Param #
================================================================
            Conv2d-1         [32, 128, 1009, 1]          43,136
         MaxPool2d-2            [32, 128, 1, 1]               0
           Flatten-3                  [32, 128]               0
            Linear-4                 [32, 1517]         195,693
            Model1-5                 [32, 1517]               0
================================================================
Total params: 238,829
Trainable params: 238,829
Non-trainable params: 0
----------------------------------------------------------------
Input size (MB): 2.62
Forward/backward pass size (MB): 32.33
Params size (MB): 0.91
Estimated Total Size (MB): 35.87
----------------------------------------------------------------


选用 multi GPU =  [2, 3] 时候： 或者 [1, 2, 3] & [0, 1, 2, 3] 结果一样：
----------------------------------------------------------------
        Layer (type)               Output Shape         Param #
================================================================
            Conv2d-1         [32, 128, 1009, 1]          43,136
         MaxPool2d-2            [32, 128, 1, 1]               0
           Flatten-3                  [32, 128]               0
            Linear-4                 [32, 1517]         195,693
            Model1-5                 [32, 1517]               0
            Conv2d-6         [32, 128, 1009, 1]          43,136
         MaxPool2d-7            [32, 128, 1, 1]               0
           Flatten-8                  [32, 128]               0
            Linear-9                 [32, 1517]         195,693
           Model1-10                 [32, 1517]               0
================================================================
Total params: 477,658
Trainable params: 477,658
Non-trainable params: 0
----------------------------------------------------------------
Input size (MB): 2.62
Forward/backward pass size (MB): 64.67
Params size (MB): 1.82
Estimated Total Size (MB): 69.12
----------------------------------------------------------------


'''








