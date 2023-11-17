
# single
cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/test_TrainALL00/test_TrainALL00_TestALL00_ss8/DeepSS2GO_Kernel32_Filter32768_Ontsall || exit
echo 1111_DeepSS2GO_Kernel32_Filter32768_Ontsall
bash step5_Diamond4CrossSpecies.sh
bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh


cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/test_TrainALL00/test_TrainALL00_TestALL00_ss8/DeepSS2GO_Kernel48_Filter16384_Ontsall || exit
echo 1111_DeepSS2GO_Kernel48_Filter16384_Ontsall
bash step5_Diamond4CrossSpecies.sh
bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh


cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/AlphaBeta/s3_AlphaBeta_TrainALL00_TestALL00_bp_aaK16F32768_ss8K32F32768 || exit
echo 22222_s3_AlphaBeta_TrainALL00_TestALL00_bp_aaK16F32768_ss8K32F32768
bash step1_cpData_Combine_Predictions_aass8.sh
bash step2_Diamond4CrossSpecies.sh
bash step3.0_FindAlphaBeta_CoarseScreening.sh  # 粗筛


cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/AlphaBeta/s3_AlphaBeta_TrainALL00_TestALL00_cc_aaK16F32768_ss8K48F16384 || exit
echo 22222_s3_AlphaBeta_TrainALL00_TestALL00_cc_aaK16F32768_ss8K48F16384
bash step1_cpData_Combine_Predictions_aass8.sh
bash step2_Diamond4CrossSpecies.sh
bash step3.0_FindAlphaBeta_CoarseScreening.sh  # 粗筛



cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/AlphaBeta/s3_AlphaBeta_TrainALL00_TestALL00_mf_aaK16F32768_ss8K32F32768 || exit
echo 22222_s3_AlphaBeta_TrainALL00_TestALL00_mf_aaK16F32768_ss8K32F32768
bash step1_cpData_Combine_Predictions_aass8.sh
bash step2_Diamond4CrossSpecies.sh
bash step3.0_FindAlphaBeta_CoarseScreening.sh  # 粗筛



echo happy_ending
echo happy_ending
