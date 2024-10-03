import argparse
import matplotlib.pyplot as plt
from Bio import SeqIO
from receptor_utils import simple_bio_seq as simple

# Plot the distribution of sequence lengths in the nominated fasta file

args = argparse.ArgumentParser(description='Plot the distribution of sequence lengths in the nominated file')
args.add_argument('infile', help='Path to the input file')
args.add_argument('output_file', help='Path to the plot file')

args = args.parse_args()

lengths = []

# find the infile file extension
ext = args.infile.split('.')[-1].lower()

if ext in ['fasta', 'fastq']:
    for rec in SeqIO.parse(args.infile, ext):
        lengths.append(len(rec.seq))
elif ext in ['csv']:
    data = simple.read_csv(args.infile)
    for row in data:
        if 'sequence' not in row:
            print('Sorry, I dont know how to count sequences in that file')
            exit(1)
        lengths.append(len(row['sequence']))
elif ext in ['tsv']:
    data = simple.read_csv(args.infile, delimiter='\t')
    for row in data:
        if 'sequence' not in row:
            print('Sorry, I dont know how to count sequences in that file')
            exit(1)
        lengths.append(len(row['sequence']))
else:
    print('Sorry, I dont know how to count sequences in that file')
    exit(1)

plt.hist(lengths, bins=50)

plt.xlabel('Sequence length')
plt.ylabel('Frequency')
plt.title('Sequence length distribution')

plt.savefig(args.output_file)

print(f'Min: {min(lengths)}, Max: {max(lengths)}, Mean: {sum(lengths) / len(lengths)}, Median: {sorted(lengths)[len(lengths) // 2]}')


