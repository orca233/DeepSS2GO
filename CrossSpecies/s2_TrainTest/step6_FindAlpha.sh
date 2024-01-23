
mkdir -p results
mkdir -p results_alpha
mkdir -p results_alphabeta




echo Find Alpha bp
python Alpha_FindAlpha.py -o bp -ncc 60 --alpha-range 25,80,1  # -ar
echo Find Alpha cc
python Alpha_FindAlpha.py -o cc -ncc 60 --alpha-range 25,80,1
echo Find Alpha mf
python Alpha_FindAlpha.py -o mf -ncc 60 --alpha-range 25,80,1



