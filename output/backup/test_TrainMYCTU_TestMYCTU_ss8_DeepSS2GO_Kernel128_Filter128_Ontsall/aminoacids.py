import numpy as np
from step0_TrainTestSetting import *


# aminoacids读取的应该是local的
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

# MAXLEN = int(config_local_dict['database']['MAXLEN'])
MAXLEN = 2000

### FS ###
aa_letter = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
ss8_letter = ['C', 'S', 'T', 'H', 'G', 'I', 'E', 'B']  # no E in aa_letter
ss3_letter = ['C', 'E', 'H']  # C = C+S+T    E = E+B     H = H+G+I

PROT_LETTER = []
if config_local_dict['database']['aa_ss'] == 'aa':
    PROT_LETTER = aa_letter
    PROT_LETTER_len = 20
elif config_local_dict['database']['aa_ss'] == 'ss8':
    PROT_LETTER = ss8_letter
    PROT_LETTER_len = 8
elif config_local_dict['database']['aa_ss'] == 'ss3':
    PROT_LETTER = ss3_letter
    PROT_LETTER_len = 3


# PROT_LETTER_len = len(PROT_LETTER)

PROT_INDEX = dict()
for i in range(len(PROT_LETTER)):
    PROT_INDEX[PROT_LETTER[i]] = i + 1
INVALID_ACIDS = {'U', 'O', 'B', 'Z', 'J', 'X', '*'}  # 这个是为了n-gram。是一种将蛋白质序列划分为连续的n个氨基酸残基的片段的表示方法


def to_onehot(seq, start=0):  # change seq to onehot，不在PROT_INDEX里面的，就会在第一列标注1
    onehot = np.zeros((MAXLEN, PROT_LETTER_len+1), dtype=np.int32)  # original:(MAXLEN, 21)
    l = min(MAXLEN, len(seq))
    for i in range(start, start + l):
        onehot[i, PROT_INDEX.get(seq[i - start], 0)] = 1
    onehot[0:start, 0] = 1
    onehot[start + l:, 0] = 1
    return onehot



# if config_local_dict['database']['aa_ss'] == 'aa':
#     PROT_LETTER = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
#     PROT_LETTER_len = len(PROT_LETTER)
#
#     PROT_INDEX = dict()
#     for i in range(len(PROT_LETTER)):
#         PROT_INDEX[PROT_LETTER[i]] = i + 1
#     INVALID_ACIDS = {'U', 'O', 'B', 'Z', 'J', 'X', '*'}
#
# else:  # aa_ss == 'ss8' 或者 ‘ss3’
#     PROT_LETTER = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',  'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'B']
#         # FS: 'V' 2 'B'
#     PROT_LETTER_len = len(PROT_LETTER)
#
#     PROT_INDEX = dict()
#     for i in range(len(PROT_LETTER)):
#         PROT_INDEX[PROT_LETTER[i]] = i + 1
#     INVALID_ACIDS = {'U', 'O', '-', 'Z', 'J', 'X', '*'}  # FS: 'B' 2 '-'





#  稍后再慢慢研究
#
# elif aa_ss == 'ss8':
#     PROT_LETTER = ['C', 'S', 'T', 'H', 'G', 'I', 'E', 'B']
#     PROT_LETTER_len = len(PROT_LETTER)
#
#     PROT_INDEX = dict()
#     for i in range(len(PROT_LETTER)):
#         PROT_INDEX[PROT_LETTER[i]] = i + 1
#     INVALID_ACIDS = set(['U', 'O', '-', 'Z', 'J', 'X', '*'])
#
# elif aa_ss == 'ss3':
#     PROT_LETTER = ['C', 'E', 'H']
#     PROT_LETTER_len = len(PROT_LETTER)
#
#     PROT_INDEX = dict()
#     for i in range(len(PROT_LETTER)):
#         PROT_INDEX[PROT_LETTER[i]] = i + 1
#     INVALID_ACIDS = set(['U', 'O', '-', 'Z', 'J', 'X', '*'])

### FS change end ###



# MAXLEN = 2000
NGRAMS = {}
for i in range(20):
    for j in range(20):
        for k in range(20):
            ngram = PROT_LETTER[i] + PROT_LETTER[j] + PROT_LETTER[k]
            index = 400 * i + 20 * j + k + 1
            NGRAMS[ngram] = index

def is_ok(seq):
    for c in seq:
        if c in INVALID_ACIDS:
            return False
    return True

def to_ngrams(seq):
    l = min(MAXLEN, len(seq) - 3)
    ngrams = np.zeros((l,), dtype=np.int32)
    for i in range(l):
        ngrams[i] = NGRAMS.get(seq[i: i + 3], 0)
    return ngrams



