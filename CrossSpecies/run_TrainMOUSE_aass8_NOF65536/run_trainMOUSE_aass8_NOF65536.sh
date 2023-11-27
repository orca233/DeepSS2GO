echo MOUSE_GPU0-1

#echo ------s2_TrainTest_TrainMOUSE_TestARATH_aa_NOF65536------
#cd ../s2_TrainTest_TrainMOUSE_TestARATH_aa_NOF65536 || exit
#python run_KernelX_FilterY.py
#
#echo ------s2_TrainTest_TrainMOUSE_TestARATH_ss8_NOF65536------
#cd ../s2_TrainTest_TrainMOUSE_TestARATH_ss8_NOF65536 || exit
#python run_KernelX_FilterY.py
#
#echo ------s2_TrainTest_TrainMOUSE_TestECOLI_aa_NOF65536------
#cd ../s2_TrainTest_TrainMOUSE_TestECOLI_aa_NOF65536 || exit
#python run_KernelX_FilterY.py
#
#echo ------s2_TrainTest_TrainMOUSE_TestECOLI_ss8_NOF65536------
#cd ../s2_TrainTest_TrainMOUSE_TestECOLI_ss8_NOF65536 || exit
#python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainMOUSE_TestHUMAN_aa_NOF65536------
cd ../s2_TrainTest_TrainMOUSE_TestMOUSE_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainMOUSE_TestHUMAN_ss8_NOF65536------
cd ../s2_TrainTest_TrainMOUSE_TestMOUSE_ss8_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainMOUSE_TestMYCTU_aa_NOF65536------
cd ../s2_TrainTest_TrainMOUSE_TestMYCTU_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainMOUSE_TestMYCTU_ss8_NOF65536------
cd ../s2_TrainTest_TrainMOUSE_TestMYCTU_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainMOUSE_TestYEAST_aa_NOF65536------
cd ../s2_TrainTest_TrainMOUSE_TestYEAST_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainMOUSE_TestYEAST_ss8_NOF65536------
cd ../s2_TrainTest_TrainMOUSE_TestYEAST_ss8_NOF65536 || exit
python run_KernelX_FilterY.py


echo happy_ending

