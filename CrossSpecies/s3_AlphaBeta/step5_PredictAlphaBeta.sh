




echo predicting bp
python AlphaBeta_PredictAlphaBeta.py -t 0.1 --in-file-aa 'data/test_data_aa.fa' --in-file-ss8 'data/test_data_ss8.fa' --test-data-file-aa 'data/test_data_aa.pkl' --test-data-file-ss8 'data/test_data_ss8.pkl' --alpha 'json' --beta 'json' -o bp  # 0.6

#echo predicting cc
#python AlphaBeta_PredictAlphaBeta.py -t 0.1 --in-file-aa 'data/test_data_aa.fa' --in-file-ss8 'data/test_data_ss8.fa' --test-data-file-aa 'data/test_data_aa.pkl' --test-data-file-ss8 'data/test_data_ss8.pkl' --alpha 'json' --beta 'json' -o cc  # 0.6
#
#echo predicting mf
#python AlphaBeta_PredictAlphaBeta.py -t 0.1 --in-file-aa 'data/test_data_aa.fa' --in-file-ss8 'data/test_data_ss8.fa' --test-data-file-aa 'data/test_data_aa.pkl' --test-data-file-ss8 'data/test_data_ss8.pkl' --alpha 'json' --beta 'json' -o mf  # 0.6








