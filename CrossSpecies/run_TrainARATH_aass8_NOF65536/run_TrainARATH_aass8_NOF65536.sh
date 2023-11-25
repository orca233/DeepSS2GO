echo TrainECOLI NOF65536, pluto, GPU0-1,,

echo ------s2_TrainTest_TrainECOLI_TestARATH_aa_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestARATH_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainECOLI_TestARATH_ss8_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestARATH_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainECOLI_TestECOLI_aa_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestECOLI_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainECOLI_TestECOLI_ss8_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestECOLI_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainECOLI_TestHUMAN_aa_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestHUMAN_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainECOLI_TestHUMAN_ss8_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestHUMAN_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainECOLI_TestMOUSE_aa_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestMOUSE_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainECOLI_TestMOUSE_ss8_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestMOUSE_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainECOLI_TestMYCTU_aa_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestMYCTU_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainECOLI_TestMYCTU_ss8_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestMYCTU_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainECOLI_TestYEAST_aa_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestYEAST_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainECOLI_TestYEAST_ss8_NOF65536------
cd ../s2_TrainTest_TrainECOLI_TestYEAST_ss8_NOF65536 || exit
python run_KernelX_FilterY.py


echo happy_ending

