# Remove gaps from a fasta file

import argparse

args = argparse.ArgumentParser(description='Remove gaps from a fasta file')
args.add_argument('input_file', help='Path to the input file')
args.add_argument('output_file', help='Path to the output file')
args = args.parse_args()

with open(args.input_file, 'r') as fi:
    with open(args.output_file, 'w') as fo:
        for line in fi:
            if line.startswith('>'):
                fo.write(line)
            else:
                line = line.replace('-', '')
                line = line.replace('.', '')
                fo.write(line)
