echo starting_run_CAFA3_BEST

cd ..

echo 111111111111111_test_TrainCAFA3_TestCAFA3_aa_DeepSS2GO_Kernel16_Filter64_Ontsall
cd test_TrainCAFA3_TestCAFA3_aa_DeepSS2GO_Kernel16_Filter64_Ontsall || exit

python step1_SplitTrainTest_Terms.py
python step2_Train.py
python step3_Test.py
python step4_pkl2fa.py
bash step5_Diamond4CrossSpecies.sh
bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..



echo 111111111111111_test_TrainCAFA3_TestCAFA3_aa_DeepSS2GO_Kernel16_Filter65536_Ontsall
cd test_TrainCAFA3_TestCAFA3_aa_DeepSS2GO_Kernel16_Filter65536_Ontsall || exit

python step1_SplitTrainTest_Terms.py
python step2_Train.py
python step3_Test.py
python step4_pkl2fa.py
bash step5_Diamond4CrossSpecies.sh
bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
bash step8_PredictAlpha.sh
cd ..




