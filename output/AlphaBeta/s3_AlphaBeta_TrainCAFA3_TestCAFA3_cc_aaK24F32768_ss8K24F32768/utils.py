import os
from collections import deque, Counter
import warnings
import pandas as pd
import numpy as np
from xml.etree import ElementTree as ET
import math
import torch
# from step0_TrainTestSetting_local_aa import *
#


'''
go.calculate_ic
GO:0004693': 0.05062607306996814, 很具体
'GO:0005575': 0.0，很祖先父辈
GO:0004518': 2.027105230162598，也很具体

如果该 GO term 没有祖先节点，那么它的信息内容为 1（即 math.log(1, 2) = 0）；
否则，找到所有祖先节点中在序列中出现次数最少的节点，设其出现次数为 min_n，则该 GO term 的信息内容为 math.log(min_n/n, 2)，其中 n 是该 GO term 在所有序列中出现的总次数。

越大，说明min_n越高，n越少：越子辈
min_n高：很多的父辈，n少自己出现次数少

=0，没祖先，最大爷
'''

BIOLOGICAL_PROCESS = 'GO:0008150'
MOLECULAR_FUNCTION = 'GO:0003674'
CELLULAR_COMPONENT = 'GO:0005575'
FUNC_DICT = {
    'cc': CELLULAR_COMPONENT,
    'mf': MOLECULAR_FUNCTION,
    'bp': BIOLOGICAL_PROCESS}

NAMESPACES = {
    'cc': 'cellular_component',
    'mf': 'molecular_function',
    'bp': 'biological_process'
}

# EXP_CODES = set(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'])
EXP_CODES = {'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}


# CAFA4 Targets
# CAFA_TARGETS = set(['287', '3702', '4577', '6239', '7227', '7955', '9606', '9823', '10090', '10116', '44689',
# '83333', '99287', '226900', '243273', '284812', '559292'])
CAFA_TARGETS = {'287', '3702', '4577', '6239', '7227', '7955', '9606', '9823', '10090', '10116', '44689', '83333',
                '99287', '226900', '243273', '284812', '559292'}


def is_cafa_target(org):
    return org in CAFA_TARGETS


def is_exp_code(code):
    return code in EXP_CODES


class Ontology(object):
    # def __init__(self, filename=params_local['path_base']+'pub_data/go.obo', with_rels=False):  # original: 'data/go.obo'   # '../../pub_data/go.obo'
    def __init__(self, filename='data/go.obo', with_rels=False):  # original: 'data/go.obo'   # '../../pub_data/go.obo'

        # print(os.getcwd())
        os.system('ls %s' % filename)
        self.ont = self.load(filename, with_rels)
        self.ic = None

    def has_term(self, term_id):
        return term_id in self.ont

    def get_term(self, term_id):
        if self.has_term(term_id):
            return self.ont[term_id]
        return None

    def calculate_ic(self, annots):
        cnt = Counter()
        for x in annots:
            cnt.update(x)
        print('---------- cnt -------------')
        # print(cnt)
        # Counter({'GO:0005575': 12037, 'GO:0110165': 11958, 'GO:0008150': 10151,
        # 'GO:0016020': 10086, 'GO:0005622': 10079, 'GO:0003674': 9125,

        self.ic = {}
        for go_id, n in cnt.items():
            # n: 在train & test 中合集出现多少次
            parents = self.get_parents(go_id)  # 找爸爸
            # print('parents = ', parents)
            # parents =  {'GO:0044238', 'GO:1901564', 'GO:0043170'} 分别出现5703,3519,4996次
            # go_id, min_n, n =  GO:0019538 3519 2855

            if len(parents) == 0:
                min_n = n
            else:
                min_n = min([cnt[x] for x in parents])  # min_n 父辈中出现最少的次数

            # print('go_id, min_n, n = ', go_id, min_n, n)  # add by FU
            self.ic[go_id] = math.log(min_n / n, 2)

    def get_ic(self, go_id):
        if self.ic is None:
            raise Exception('Not yet calculated')
        if go_id not in self.ic:
            return 0.0
        return self.ic[go_id]


    ### FS add ### 计算到达根节点最短的距离
    '''
    添加一个距离参数来记录每个节点到指定GO term的距离。
    将节点和距离作为一个元组添加到队列中。
    在每次迭代中，检查当前节点是否等于根节点，如果是，则返回当前距离值。
    如果队列为空时仍未找到根节点，表示指定GO term无法到达根节点，返回None。
    '''

    # 这是对单个的term_id求到根节点最短距离，calculate_ic是对list的所有term求ic值
    def calculate_shortest_path_to_root(self, term_id):
        if term_id not in self.ont:
            return None

        q = deque()
        q.append((term_id, 0))

        while q:
            t_id, distance = q.popleft()

            if t_id == "GO:0008150" or t_id == 'GO:0003674' or t_id =='GO:0005575':
                return distance

            for parent_id in self.ont[t_id]["is_a"]:
                if parent_id in self.ont:
                    q.append((parent_id, distance + 1))

        return None

    def get_shortest_path_to_root(self, term_id):  # 没啥用
        shortest_path = self.calculate_shortest_path_to_root(term_id)
        if shortest_path is None:
            return 0
        return shortest_path

    ### FS add -- done ###

    def load(self, filename, with_rels):
        ont = dict()
        obj = None
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line == '[Term]':
                    if obj is not None:
                        ont[obj['id']] = obj
                    obj = dict()
                    obj['is_a'] = list()
                    obj['part_of'] = list()
                    obj['regulates'] = list()
                    obj['alt_ids'] = list()
                    obj['is_obsolete'] = False
                    continue
                elif line == '[Typedef]':
                    if obj is not None:
                        ont[obj['id']] = obj
                    obj = None
                else:
                    if obj is None:
                        continue
                    l = line.split(": ")
                    if l[0] == 'id':
                        obj['id'] = l[1]
                    elif l[0] == 'alt_id':
                        obj['alt_ids'].append(l[1])
                    elif l[0] == 'namespace':
                        obj['namespace'] = l[1]
                    elif l[0] == 'is_a':
                        obj['is_a'].append(l[1].split(' ! ')[0])
                    elif with_rels and l[0] == 'relationship':
                        it = l[1].split()
                        # add all types of relationships
                        obj['is_a'].append(it[1])
                    elif l[0] == 'name':
                        obj['name'] = l[1]
                    elif l[0] == 'is_obsolete' and l[1] == 'true':
                        obj['is_obsolete'] = True
            if obj is not None:
                ont[obj['id']] = obj
        for term_id in list(ont.keys()):
            for t_id in ont[term_id]['alt_ids']:
                ont[t_id] = ont[term_id]
            if ont[term_id]['is_obsolete']:
                del ont[term_id]
        for term_id, val in ont.items():
            if 'children' not in val:
                val['children'] = set()
            for p_id in val['is_a']:
                if p_id in ont:
                    if 'children' not in ont[p_id]:
                        ont[p_id]['children'] = set()
                    ont[p_id]['children'].add(term_id)
        return ont

    def get_ancestors(self, term_id):
        if term_id not in self.ont:
            return set()
        term_set = set()
        q = deque()
        q.append(term_id)
        while (len(q) > 0):
            t_id = q.popleft()
            if t_id not in term_set:
                term_set.add(t_id)
                for parent_id in self.ont[t_id]['is_a']:
                    if parent_id in self.ont:
                        q.append(parent_id)
        return term_set

    def get_parents(self, term_id):
        if term_id not in self.ont:
            return set()
        term_set = set()
        for parent_id in self.ont[term_id]['is_a']:
            if parent_id in self.ont:
                term_set.add(parent_id)
        return term_set


    '''
    定义一个函数，把一个GO set中所有的父类都去掉，只留下最根部的
    注意：可以get_children，但是不能get_all_children，会无穷多个
    如果GO_set里有一个爷爷辈，一个孙子辈，没有父亲辈分，也要把爷爷挑出来（GO:0009581 & GO:0050896）
    此时就要用get_ancestors了
    
    GPT：
    根据下面两个gene ontology的function，写一个python，把给定GO_set集合中所有的祖先节点都剔除，只保留底层子节点term。
    通过'is_a'判断父辈子辈关联。
    比如给定GO_set={GO:0048609, GO:0060378, GO:0000003, GO:0001964, GO:0009605}
    剔除所有祖先节点剩下集合rest_set={GO:0060378, GO:0001964}    
    '''

    ### FU added ###  删除所有祖先，只保留根节点。
    # 这个要对stage1中的pkl进行选择了，要在swissprot.pkl中添加一列，类似于prop_annot，这个新的一列是rm_ancestors_annot
    # 也可以对swissprot_clean_ALL00_aa.pkl进行添加列，而后进行step7 & 8

    def rm_ancestors(self, go_set):
        root_terms = set()
        for term_id in go_set:
            ancestors = self.get_ancestors(term_id) - {term_id}  # func get_ancestor 是包含底层子节点，所以要减去改term_id
            root_terms.update(ancestors)

        print('root_terms=', root_terms)
        return go_set - root_terms
    ### FU added -- done ###




    def get_namespace_terms(self, namespace):
        terms = set()
        for go_id, obj in self.ont.items():
            if obj['namespace'] == namespace:
                terms.add(go_id)
        return terms

    def get_namespace(self, term_id):
        return self.ont[term_id]['namespace']

    def get_term_set(self, term_id):
        if term_id not in self.ont:
            return set()
        term_set = set()
        q = deque()
        q.append(term_id)
        while len(q) > 0:
            t_id = q.popleft()
            if t_id not in term_set:
                term_set.add(t_id)
                for ch_id in self.ont[t_id]['children']:
                    q.append(ch_id)
        return term_set


def read_fasta(filename):
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
                    seq = ''
                inf = line[1:]
            else:
                seq += line
        seqs.append(seq)
        info.append(inf)
    return info, seqs


class DataGenerator(object):

    def __init__(self, batch_size, is_sparse=False):
        self.batch_size = batch_size
        self.is_sparse = is_sparse

    def fit(self, inputs, targets=None):
        self.start = 0
        self.inputs = inputs
        self.targets = targets
        if isinstance(self.inputs, tuple) or isinstance(self.inputs, list):
            self.size = self.inputs[0].shape[0]
        else:
            self.size = self.inputs.shape[0]
        self.has_targets = targets is not None

    def __next__(self):
        return self.next()

    def reset(self):
        self.start = 0

    def next(self):
        if self.start < self.size:
            batch_index = np.arange(
                self.start, min(self.size, self.start + self.batch_size))
            if isinstance(self.inputs, tuple) or isinstance(self.inputs, list):
                res_inputs = []
                for inp in self.inputs:
                    if self.is_sparse:
                        res_inputs.append(
                            inp[batch_index, :].toarray())
                    else:
                        res_inputs.append(inp[batch_index, :])
            else:
                if self.is_sparse:
                    res_inputs = self.inputs[batch_index, :].toarray()
                else:
                    res_inputs = self.inputs[batch_index, :]
            self.start += self.batch_size
            if self.has_targets:
                if self.is_sparse:
                    labels = self.targets[batch_index, :].toarray()
                else:
                    labels = self.targets[batch_index, :]
                return (res_inputs, labels)
            return res_inputs
        else:
            self.reset()
            return self.next()

# FS add  file_path为go.obo，把所有的go_term提取出来变成一个set
def extract_go_terms(file_path):
    go_terms = set()

    with open(file_path, 'r') as file:
        term_id = None
        for line in file:
            line = line.strip()
            if line == "[Term]":
                if term_id is not None:
                    go_terms.add(term_id)
                term_id = None
            elif line.startswith("id:"):
                term_id = line.split("id:")[1].strip()

    if term_id is not None:
        go_terms.add(term_id)

    return go_terms





class EarlyStopping:
    """Early stops the training if validation loss doesn't improve after a given patience."""

    def __init__(self, patience=7, verbose=False, delta=0):
        """
        Args:
            patience (int): How long to wait after last time validation loss improved.
                            上次验证集损失值改善后等待几个epoch
                            Default: 7
            verbose (bool): If True, prints a message for each validation loss improvement.
                            如果是True，为每个验证集损失值改善打印一条信息
                            Default: False
            delta (float): Minimum change in the monitored quantity to qualify as an improvement.
                            监测数量的最小变化，以符合改进的要求
                            Default: 0
        """
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.delta = delta

    def __call__(self, val_loss, model):

        score = -val_loss

        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
        elif score < self.best_score + self.delta:
            self.counter += 1
            # print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
        else:
            self.best_score = score
            self.save_checkpoint(val_loss, model)
            self.counter = 0

    def save_checkpoint(self, val_loss, model):
        '''
        Saves model when validation loss decrease.
        验证损失减少时保存模型。
        '''
        if self.verbose:
            # print(f'Valid_loss drops ({self.val_loss_min:.6f} --> {val_loss:.6f}).  Saving model ...')
            print(f'Valid_loss drops.  Saving model ...')
        torch.save(model.state_dict(),
                   'data/model_checkpoint.pth')  # 存迄今最优模型参数， model_file='data/model_checkpoint.pth'
        # torch.save(model, 'finish_model.pkl') # 这里会存储迄今最优的模型???
        self.val_loss_min = val_loss




def pkl2fa(input_file, output_file):  # DeepGOPlus 原文件为 diamond_data.py
    # Load interpro data
    df = pd.read_pickle(input_file)
    print('pkl length', len(df))
    with open(output_file, 'w') as f:
        for row in df.itertuples():
            f.write('>' + row.proteins + '\n')
            f.write(row.sequences + '\n')





