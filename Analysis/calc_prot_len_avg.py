# 用于统计每条蛋白质长度并计算平均值、中位数和方差的Python程序


from Bio import SeqIO
import numpy as np

def calculate_statistics(sequence_lengths):
    average_length = np.mean(sequence_lengths)
    median_length = np.median(sequence_lengths)
    variance_length = np.var(sequence_lengths)
    return average_length, median_length, variance_length

def main(fasta_file):
    sequence_lengths = []

    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            sequence_length = len(record.seq)
            sequence_lengths.append(sequence_length)
            print(f"Protein {record.id}: Length {sequence_length}")

    average_length, median_length, variance_length = calculate_statistics(sequence_lengths)

    print("\nSummary:")
    print(f"Average Length: {average_length}")
    print(f"Median Length: {median_length}")
    print(f"Variance Length: {variance_length}")

if __name__ == "__main__":
    fasta_file = "/home/fsong/work/py_proj/prot_algo/DeepSS2GO/pub_data/data_cafa3/CAFA3_train_data_clean_aa.fa"  # Replace with the actual path to your fasta file
    main(fasta_file)

