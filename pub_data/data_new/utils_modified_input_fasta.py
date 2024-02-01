'''
README

1. change symbol '.' 2 '_' in each seq name (fasta from NCBI)
2. multi-line 2 one-line

python utils_modified_input_fasta.py input.fasta output.fasta

Contributed by Dawei Huang @ SUSTech
'''


import argparse

parser = argparse.ArgumentParser(description='Process a FASTA file.')
parser.add_argument('input_file', help='Input FASTA file')
parser.add_argument('output_file', help='Output FASTA file')
args = parser.parse_args()

with open(args.input_file, 'r') as file:
    lines = file.readlines()

with open(args.output_file, 'w') as file:
    for line in lines:
        if line.startswith('>'):
            # remove space
            line = line.split(' ')[0]
            # "." was replaced by "_"
            line = '\n' + line.replace('.', '_').strip() + '\n'
            file.write(line)
        else:
            # removes any leading and trailing whitespace
            file.write(line.strip())

with open(args.output_file, 'r+') as file:
    content = file.read()
    file.seek(0)
    file.write(content.lstrip('\n'))
