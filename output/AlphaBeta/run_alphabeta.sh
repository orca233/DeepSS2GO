
cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/AlphaBeta/s3_AlphaBeta_TrainCAFA3_TestCAFA3_bp_aaK8F256_ss8K96F65536 || exit
echo step3.1

bash step3.1_FindAlphaBeta_FineScreening.sh


cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/AlphaBeta/s3_AlphaBeta_TrainCAFA3_TestCAFA3_cc_aaK24F32768_ss8K24F32768 || exit
echo 11111111111_s3_AlphaBeta_TrainCAFA3_TestCAFA3_cc_aaK24F32768_ss8K24F32768

bash step3.0_FindAlphaBeta_CoarseScreening.sh



cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/AlphaBeta/s3_AlphaBeta_TrainCAFA3_TestCAFA3_mf_aaK16F32768_ss8K32F65536 || exit
echo 11111111111_s3_AlphaBeta_TrainCAFA3_TestCAFA3_mf_aaK16F32768_ss8K32F65536

bash step3.0_FindAlphaBeta_CoarseScreening.sh



echo happy_ending
