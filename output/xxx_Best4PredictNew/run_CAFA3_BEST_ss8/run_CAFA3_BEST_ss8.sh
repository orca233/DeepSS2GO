echo starting_run_CAFA3_BEST

cd ..

echo 111111111111111_test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel8_Filter512_Ontsall
cd test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel8_Filter512_Ontsall || exit

python step1_SplitTrainTest_Terms.py
python step2_Train.py
python step3_Test.py
python step4_pkl2fa.py
bash step5_Diamond4CrossSpecies.sh
bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..



echo 111111111111111_test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel32_Filter8192_Ontsall
cd test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel32_Filter8192_Ontsall || exit

python step1_SplitTrainTest_Terms.py
python step2_Train.py
python step3_Test.py
python step4_pkl2fa.py
bash step5_Diamond4CrossSpecies.sh
bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..


echo 111111111111111_test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel24_Filter16384_Ontsall
cd test_TrainCAFA3_TestCAFA3_ss8_DeepSS2GO_Kernel24_Filter16384_Ontsall || exit

python step1_SplitTrainTest_Terms.py
python step2_Train.py
python step3_Test.py
python step4_pkl2fa.py
bash step5_Diamond4CrossSpecies.sh
bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..

