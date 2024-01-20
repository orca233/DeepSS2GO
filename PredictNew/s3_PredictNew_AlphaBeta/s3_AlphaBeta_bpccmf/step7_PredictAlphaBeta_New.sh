
# Predict
echo predicting bp
python AlphaBeta_PredictAlphaBeta.py -t 0.1 --in-file-aa 'data/new_clean_aa.fa' --in-file-ss8 'data/new_clean_ss8.fa' --test-data-file-aa 'data/new_clean_aa.pkl' --test-data-file-ss8 'data/new_clean_ss8.pkl' --alpha 'json' --beta 'json' -o bp  # 0.6

echo predicting cc
python AlphaBeta_PredictAlphaBeta.py -t 0.1 --in-file-aa 'data/new_clean_aa.fa' --in-file-ss8 'data/new_clean_ss8.fa' --test-data-file-aa 'data/new_clean_aa.pkl' --test-data-file-ss8 'data/new_clean_ss8.pkl' --alpha 'json' --beta 'json' -o cc  # 0.6

echo predicting mf
python AlphaBeta_PredictAlphaBeta.py -t 0.02 --in-file-aa 'data/new_clean_aa.fa' --in-file-ss8 'data/new_clean_ss8.fa' --test-data-file-aa 'data/new_clean_aa.pkl' --test-data-file-ss8 'data/new_clean_ss8.pkl' --alpha 'json' --beta 'json' -o mf  # 0.6


# Attention, same as: s3_AlphaBeta_TrainALL00_TestALL00_mf_aaK16F32768_ss8K32F3276
# step7，activate all three echo




