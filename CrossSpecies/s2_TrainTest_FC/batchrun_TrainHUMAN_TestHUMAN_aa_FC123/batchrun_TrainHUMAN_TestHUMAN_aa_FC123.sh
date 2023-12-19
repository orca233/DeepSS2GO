echo TrainHUMAN_aa_FC, apus, GPU0-1

echo ------s2_TrainTest_TrainHUMAN_TestHUMAN_aa_FC1------
cd ../s2_TrainTest_TrainHUMAN_TestHUMAN_aa_FC1 || exit
python run_KernelX_FilterY.py


echo ------s2_TrainTest_TrainHUMAN_TestHUMAN_aa_FC2------
cd ../s2_TrainTest_TrainHUMAN_TestHUMAN_aa_FC2 || exit
python run_KernelX_FilterY.py


echo ------s2_TrainTest_TrainHUMAN_TestHUMAN_aa_FC3------
cd ../s2_TrainTest_TrainHUMAN_TestHUMAN_aa_FC3 || exit
python run_KernelX_FilterY.py


echo happy_ending
