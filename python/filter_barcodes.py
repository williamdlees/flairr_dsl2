# Enumerate and plot barcode usage - how many barcodes, how many times each barcode was used

import matplotlib.pyplot as plt
from Bio import SeqIO
import random
import argparse


# Parse the command line arguments

parser = argparse.ArgumentParser(description='Enumerate and plot barcode usage')
parser.add_argument('infile', help='FASTQ file to process')
parser.add_argument('outfile', help='Output file for reduced FASTQ')
parser.add_argument('plotfile', help='Output file for plot')
parser.add_argument('--barcode_field', help='Field name for barcode in FASTQ description', default='BARCODE')
parser.add_argument('--cutoff', help='Cutoff for number of records per barcode', type=int, default=100)
parser.add_argument('--limit', help='Limit on number of records to process', type=int, default=0)
args = parser.parse_args()


# Function to read FASTQ records from a file
def read_fastq(file_path):
    records = []
    count = 0
    for record in SeqIO.parse(file_path, "fastq"):
        records.append(record)
        count += 1
        if args.limit > 0 and count >= args.limit:
            break
    return records


def find_barcode(rec):
    for field in rec.description.split('|'):
        if field.startswith(args.barcode_field + '='):
            return field.split('=')[1]
    return None


recs = read_fastq(args.infile)

# cluster the reads by barcode

barcode_recs = {}
rec_count = 0

for rec in recs:
    rec_count += 1
    barcode = find_barcode(rec)
    if barcode not in barcode_recs:
        barcode_recs[barcode] = []
    barcode_recs[barcode].append(rec)

# Extract the counts for plotting
counts = [len(x) for x in barcode_recs.values() if len(x) > args.cutoff]

# Plot the histogram
plt.hist(counts, bins=10, edgecolor='black')
plt.xlabel(f'Number of records per barcode, for barcodes with > {args.cutoff} records')
plt.ylabel('Frequency')
plt.title('Histogram of Barcode Usage')

# plot to a pdf
plt.savefig(args.plotfile)

# print the total number of records

print(f'Total number of records: {rec_count}')

# print the number of records that would be retained after cutoff

retained = 0
for bcs in barcode_recs.values():
    if len(bcs) <= args.cutoff:
        retained += len(bcs)
    else:
        retained += args.cutoff

print(f'Number of records retained after cutoff of {args.cutoff}: {retained}')

# reduce the number of reads in each cluster to the cutoff by drawing reads at random with no replacement

for barcode, recs in barcode_recs.items():
    if len(recs) > args.cutoff:
        barcode_recs[barcode] = random.sample(recs, args.cutoff)

# write the reduced records to a new file

with open(args.outfile, 'w') as fo:
    for barcode, recs in barcode_recs.items():
        for rec in recs:
            SeqIO.write(rec, fo, "fastq")
