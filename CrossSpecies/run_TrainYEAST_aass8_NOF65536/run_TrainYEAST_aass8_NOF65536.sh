echo TrainMYCTU NOF65536, haumea, GPU0-1

echo ------s2_TrainTest_TrainMYCTU_TestARATH_aa_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestARATH_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainMYCTU_TestARATH_ss8_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestARATH_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainMYCTU_TestECOLI_aa_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestECOLI_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainMYCTU_TestECOLI_ss8_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestECOLI_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainMYCTU_TestHUMAN_aa_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestHUMAN_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainMYCTU_TestHUMAN_ss8_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestHUMAN_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainMYCTU_TestMOUSE_aa_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestMOUSE_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainMYCTU_TestMOUSE_ss8_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestMOUSE_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainMYCTU_TestMYCTU_aa_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestMYCTU_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainMYCTU_TestMYCTU_ss8_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestMYCTU_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainMYCTU_TestYEAST_aa_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestYEAST_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainMYCTU_TestYEAST_ss8_NOF65536------
cd ../s2_TrainTest_TrainMYCTU_TestYEAST_ss8_NOF65536 || exit
python run_KernelX_FilterY.py


echo happy_ending

