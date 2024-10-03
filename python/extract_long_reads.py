import argparse
from receptor_utils import simple_bio_seq as simple
from Bio import SeqIO

# Extract reads over the specified length into a fasta file

args = argparse.ArgumentParser(description='Plot the distribution of sequence lengths in the nominated fasta file')
args.add_argument('min_length', help='Minimum length of reads to extract')
args.add_argument('infile', help='Path to the input file')
args.add_argument('output_file', help='Path to the plot file')

max_reads = 50

args = args.parse_args()

long_reads = {}

# find the infile file extension
ext = args.infile.split('.')[-1].lower()

for rec in SeqIO.parse(args.infile, ext):
    if len(rec.seq) >= int(args.min_length):
        long_reads[rec.id] = rec.seq
        if len(long_reads) >= max_reads:
            break

simple.write_fasta(args.output_file, long_reads)

print(f'Extracted {len(long_reads)} reads over {args.min_length} bases')


