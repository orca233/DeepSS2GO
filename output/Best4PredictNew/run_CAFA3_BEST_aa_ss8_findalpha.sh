echo starting_run_CAFA3_BEST


echo 111111111111111_test_TrainCAFA3_TestCAFA3_aa_DeepSS2GO_Kernel8_Filter32768_Ontsall
cd test_TrainCAFA3_TestCAFA3_aa_DeepSS2GO_Kernel8_Filter32768_Ontsall || exit





bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..



echo 111111111111111_test_TrainCAFA3_TestCAFA3_aa_DeepSS2GO_Kernel16_Filter32768_Ontsall
cd test_TrainCAFA3_TestCAFA3_aa_DeepSS2GO_Kernel16_Filter32768_Ontsall || exit





bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..


echo 111111111111111_test_TrainCAFA3_TestCAFA3_aa_DeepSS2GO_Kernel24_Filter32768_Ontsall
cd test_TrainCAFA3_TestCAFA3_aa_DeepSS2GO_Kernel24_Filter32768_Ontsall || exit





bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..






echo starting_run_CAFA3_BEST


echo 111111111111111_test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel24_Filter32768_Ontsall
cd test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel24_Filter32768_Ontsall || exit





bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..



echo 111111111111111_test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel64_Filter4096_Ontsall
cd test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel64_Filter4096_Ontsall || exit





bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..


echo 111111111111111_test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel64_Filter32768_Ontsall
cd test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel64_Filter32768_Ontsall || exit





bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..

