echo TrainARATH NOF65536, orcus, GPU2-3


echo ------s2_TrainTest_TrainARATH_TestMYCTU_ss8_NOF65536_reverse01------
cd ../s2_TrainTest_TrainARATH_TestMYCTU_ss8_NOF65536_reverse01 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainARATH_TestMYCTU_aa_NOF65536_reverse01------
cd ../s2_TrainTest_TrainARATH_TestMYCTU_aa_NOF65536_reverse01 || exit
python run_KernelX_FilterY.py




#
echo ------s2_TrainTest_TrainARATH_TestMOUSE_ss8_NOF65536_reverse01------
cd ../s2_TrainTest_TrainARATH_TestMOUSE_ss8_NOF65536_reverse01 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainARATH_TestMOUSE_aa_NOF65536_reverse01------
cd ../s2_TrainTest_TrainARATH_TestMOUSE_aa_NOF65536_reverse01 || exit
python run_KernelX_FilterY.py




#




#





#








echo happy_ending

