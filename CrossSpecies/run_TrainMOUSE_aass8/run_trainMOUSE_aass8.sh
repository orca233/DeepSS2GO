echo HUMAN_GPU2-3

echo ------s2_TrainTest_TrainHUMAN_TestARATH_aa------
cd ../s2_TrainTest_TrainHUMAN_TestARATH_aa || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestARATH_ss8------
cd ../s2_TrainTest_TrainHUMAN_TestARATH_ss8 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestECOLI_aa------
cd ../s2_TrainTest_TrainHUMAN_TestECOLI_aa || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestECOLI_ss8------
cd ../s2_TrainTest_TrainHUMAN_TestECOLI_ss8 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainHUMAN_TestMOUSE_aa------
cd ../s2_TrainTest_TrainHUMAN_TestMOUSE_aa || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestMOUSE_ss8------
cd ../s2_TrainTest_TrainHUMAN_TestMOUSE_ss8 || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestMYCTU_aa------
cd ../s2_TrainTest_TrainHUMAN_TestMYCTU_aa || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestMYCTU_ss8------
cd ../s2_TrainTest_TrainHUMAN_TestMYCTU_ss8 || exit
python run_KernelX_FilterY.py



echo ------s2_TrainTest_TrainHUMAN_TestYEAST_aa------
cd ../s2_TrainTest_TrainHUMAN_TestYEAST_aa || exit
python run_KernelX_FilterY.py

echo ------s2_TrainTest_TrainHUMAN_TestYEAST_ss8------
cd ../s2_TrainTest_TrainHUMAN_TestYEAST_ss8 || exit
python run_KernelX_FilterY.py


echo happy_ending

