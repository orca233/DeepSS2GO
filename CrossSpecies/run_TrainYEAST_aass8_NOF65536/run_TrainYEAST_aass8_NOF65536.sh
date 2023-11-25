echo TrainYEAST NOF65536, haumea, GPU2-3

echo ------s2_TrainTest_TrainYEAST_TestARATH_aa_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestARATH_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainYEAST_TestARATH_ss8_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestARATH_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainYEAST_TestECOLI_aa_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestECOLI_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainYEAST_TestECOLI_ss8_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestECOLI_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainYEAST_TestHUMAN_aa_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestHUMAN_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainYEAST_TestHUMAN_ss8_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestHUMAN_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainYEAST_TestMOUSE_aa_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestMOUSE_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainYEAST_TestMOUSE_ss8_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestMOUSE_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainYEAST_TestMYCTU_aa_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestMYCTU_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainYEAST_TestMYCTU_ss8_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestMYCTU_ss8_NOF65536 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainYEAST_TestYEAST_aa_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestYEAST_aa_NOF65536 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainYEAST_TestYEAST_ss8_NOF65536------
cd ../s2_TrainTest_TrainYEAST_TestYEAST_ss8_NOF65536 || exit
python run_KernelX_FilterY.py


echo happy_ending

