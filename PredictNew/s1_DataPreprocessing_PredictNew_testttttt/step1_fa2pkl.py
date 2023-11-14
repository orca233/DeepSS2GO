# new_aa.fa 是个新文件，未知的，待测的，要测GO。

import os
from step0_DataPreprocessingSetting import *

input_file = path_base + 'pub_data/data_new/new_aa.fa'
output_file = path_base + 'pub_data/data_new/new_aa.pkl'

os.system('python Alpha_fa2pkl.py -inf %s -ouf %s' % (input_file, output_file))





