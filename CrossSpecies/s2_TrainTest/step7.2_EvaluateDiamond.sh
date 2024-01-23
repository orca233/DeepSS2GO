mkdir -p results
mkdir -p results_alpha
mkdir -p results_alphabeta

echo evaluate bp
python Alpha_EvaluateAlpha.py -o bp --alpha 0 # 0.6
echo evaluate cc
python Alpha_EvaluateAlpha.py -o cc --alpha 0
echo evaluate mf
python Alpha_EvaluateAlpha.py -o mf --alpha 0


