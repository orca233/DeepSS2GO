

cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/test_CAFA3_round5_ss8/DeepSS2GO_Kernel96_Filter65536_Ontsall || exit
echo 11111111111_DeepSS2GO_Kernel96_Filter65536_Ontsall

# bash step5_Diamond4CrossSpecies.sh
bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
#bash step8_PredictAlpha.sh


cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/test_CAFA3_round5_ss8/DeepSS2GO_Kernel24_Filter32768_Ontsall || exit
echo 11111111_DeepSS2GO_Kernel24_Filter32768_Ontsall

# bash step5_Diamond4CrossSpecies.sh
bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
#bash step8_PredictAlpha.sh


cd /home/fsong/work/py_proj/prot_algo/DeepSS2GO/output/test_CAFA3_round5_ss8/DeepSS2GO_Kernel32_Filter65536_Ontsall || exit
echo 111111_DeepSS2GO_Kernel32_Filter65536_Ontsal

# bash step5_Diamond4CrossSpecies.sh
bash step6_FindAlpha.sh
bash step7_EvaluateAlpha.sh
#bash step8_PredictAlpha.sh


echo happy_ending
