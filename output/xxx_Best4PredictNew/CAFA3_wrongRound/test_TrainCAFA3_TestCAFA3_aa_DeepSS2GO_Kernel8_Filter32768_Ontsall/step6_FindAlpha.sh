
mkdir -p results
mkdir -p results_alpha
mkdir -p results_alphabeta


#####################################################
# find_alpha  只是找到最佳alpha
# 根据实际情况，不一定要计算所有的 ont

echo Find Alpha bp
python Alpha_FindAlpha.py -o bp -ncc 50 --alpha-range 25,80,1  # -ar 三个数字之间不要有空格
echo Find Alpha cc
python Alpha_FindAlpha.py -o cc -ncc 50 --alpha-range 25,80,1
echo Find Alpha mf
python Alpha_FindAlpha.py -o mf -ncc 50 --alpha-range 25,80,1



