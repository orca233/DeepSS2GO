mkdir -p results
mkdir -p results_alpha
mkdir -p results_alphabeta

echo evaluate bp
python Alpha_EvaluateWithoutAlpha.py -o bp --alpha 1 # 0.6
echo evaluate cc
python Alpha_EvaluateWithoutAlpha.py -o cc --alpha 1
echo evaluate mf
python Alpha_EvaluateWithoutAlpha.py -o mf --alpha 1





