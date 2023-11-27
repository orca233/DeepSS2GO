echo HUMAN_GPU2-3, workstation X

#echo ------s2_TrainTest_TrainHUMAN_TestARATH_aa------
#cd ../s2_TrainTest_TrainHUMAN_TestARATH_aa || exit
#python run_KernelX_FilterY.py
#
#echo ------s2_TrainTest_TrainHUMAN_TestARATH_ss8------
#cd ../s2_TrainTest_TrainHUMAN_TestARATH_ss8 || exit
#python run_KernelX_FilterY.py
#
#echo ------s2_TrainTest_TrainHUMAN_TestECOLI_aa------
#cd ../s2_TrainTest_TrainHUMAN_TestECOLI_aa || exit
#python run_KernelX_FilterY.py
#
#echo ------s2_TrainTest_TrainHUMAN_TestECOLI_ss8------
#cd ../s2_TrainTest_TrainHUMAN_TestECOLI_ss8 || exit
#python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainHUMAN_TestMOUSE_aa_NOF65536------
cd ../s2_TrainTest_TrainHUMAN_TestMOUSE_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestMOUSE_ss8_NOF65536------
cd ../s2_TrainTest_TrainHUMAN_TestMOUSE_ss8_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestMYCTU_aa_NOF65536------
cd ../s2_TrainTest_TrainHUMAN_TestMYCTU_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestMYCTU_ss8_NOF65536------
cd ../s2_TrainTest_TrainHUMAN_TestMYCTU_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainHUMAN_TestYEAST_aa_NOF65536------
cd ../s2_TrainTest_TrainHUMAN_TestYEAST_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestYEAST_ss8_NOF65536------
cd ../s2_TrainTest_TrainHUMAN_TestYEAST_ss8_NOF65536 || exit
python run_KernelX_FilterY.py


echo happy_ending

