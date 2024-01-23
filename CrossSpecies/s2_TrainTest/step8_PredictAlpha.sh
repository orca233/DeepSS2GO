
mkdir -p results
mkdir -p results_alpha
mkdir -p results_alphabeta


echo predicting bp
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/test_data.fa' --test-data-file 'data/test_data.pkl' -o bp --alpha 'json' # 0.6

echo predicting cc
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/test_data.fa' --test-data-file 'data/test_data.pkl' -o cc --alpha 'json' # 0.6

echo predicting mf
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/test_data.fa' --test-data-file 'data/test_data.pkl' -o mf --alpha 'json' # 0.6






