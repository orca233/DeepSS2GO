#!/usr/bin/env python

import os
from step0_TrainTestSetting_local import *

aa_ss = params_local['aa_ss']


#
if aa_ss == 'aa':
    print('predicting aa bp...')
    os.system("python Alpha_PredictAlpha.py -t 0.02 --in-file 'data/new_clean_aa.fa' --test-data-file 'data/new_clean_aa.pkl' -o bp --alpha 'json'")

    print('predicting aa cc...')
    os.system("python Alpha_PredictAlpha.py -t 0.02 --in-file 'data/new_clean_aa.fa' --test-data-file 'data/new_clean_aa.pkl' -o cc --alpha 'json'")

    print('predicting aa mf...')
    os.system("python Alpha_PredictAlpha.py -t 0.02 --in-file 'data/new_clean_aa.fa' --test-data-file 'data/new_clean_aa.pkl' -o mf --alpha 'json'")



elif aa_ss == 'ss8':
    print('predicting ss8 bp...')
    os.system("python Alpha_PredictAlpha.py -t 0.02 --in-file 'data/new_clean_ss8.fa' --test-data-file 'data/new_clean_ss8.pkl' -o bp --alpha 'json'")

    print('predicting ss8 cc...')
    os.system("python Alpha_PredictAlpha.py -t 0.02 --in-file 'data/new_clean_ss8.fa' --test-data-file 'data/new_clean_ss8.pkl' -o cc --alpha 'json'")

    print('predicting ss8 mf...')
    os.system("python Alpha_PredictAlpha.py -t 0.02 --in-file 'data/new_clean_ss8.fa' --test-data-file 'data/new_clean_ss8.pkl' -o mf --alpha 'json'")





