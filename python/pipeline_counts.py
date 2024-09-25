# Create a table summarising read counts at each pipeline stage for each sample

import toml
import glob
import os
import argparse
import csv
import re
import fnmatch


# find number of reads in a fastq file
def count_fastq_reads(file):
    with open(file) as f:
        return sum(1 for line in f) // 4


# find number of reads in a fasta file
def count_fasta_reads(file):
    with open(file) as f:
        return sum(1 for line in f if line.startswith('>'))


# find number of reads listed in a tab file
def count_tab_reads(file):
    with open(file) as f:
        return sum(1 for line in f) - 1


# Write csv file given a list of dicts. Fieldnames are taken from the first row
def write_csv(file: str, rows: list, delimiter: str = None, scan_all: bool = False):
    if not rows:
        return

    if scan_all:
        fieldnames = list()
        for row in rows:
            for key in row.keys():
                if key not in fieldnames:
                    fieldnames.append(key)       # use lists rather set so we can retain ordering 
    else:
        fieldnames = rows[0].keys()

    with open(file, 'w', newline='') as fo:
        if delimiter:
            writer = csv.DictWriter(fo, fieldnames=list(fieldnames), delimiter=delimiter)
        else:
            writer = csv.DictWriter(fo, fieldnames=fieldnames)

        writer.writeheader()
        
        for row in rows:
            writer.writerow(row)


def find_samples(dirs, config):
    samples = []
    for _, file_template in config['sample_files']:
        lh, rh = file_template.split('{sample}')
        file_template = file_template.replace('{sample}', '*')
        for file in glob.glob(os.path.join(config['sample_dir'], file_template)):
            # find the value of {sample} in the matched template     
            file = file.replace('\\', '/')  
            sample = file.split(lh)[1].split(rh)[0]
            if sample not in samples:
                samples.append(sample)
 
    return sorted(samples)


args = argparse.ArgumentParser(description='Summarise read counts at each pipeline stage for each sample')
args.add_argument('config_file', help='Path to the config file')
args.add_argument('results_file', help='Path to the results file')
args = args.parse_args()

config = toml.load(args.config_file)

results = []
results_perc = []


print('Processing samples from directory ', config['sample_dir'])
sample_dir = config['sample_dir']

samples = find_samples(sample_dir, config)

print(f'Processing {sample_dir}...')

for sample_name in samples:
    rec = {'sample': sample_name}
    rec_perc = {'sample': sample_name}
    first_count = -1

    for stage, stage_filespec in config['sample_files']:
        stage_filename_glob = glob.glob(os.path.join(sample_dir, stage_filespec).replace('{sample}', sample_name))

        if len(stage_filename_glob) == 0:
            rec[stage] = 0
            rec_perc[stage] = 0
            continue
        elif len(stage_filename_glob) > 1:
            print(f'Error: multiple files match {sample_name} {stage}')
            exit(1)

        stage_filename = stage_filename_glob[0]
        file_type = os.path.splitext(stage_filename)[1]

        if file_type == '.fastq':
            num_reads = count_fastq_reads(os.path.join(stage_filename))
        elif file_type == '.fasta':
            num_reads = count_fasta_reads(os.path.join(stage_filename))
        elif file_type == '.tab' or file_type == '.tsv' or file_type == '.csv':
            num_reads = count_tab_reads(os.path.join(stage_filename))
        else:
            print(f'Error: unsupported file type {file_type} in config file')
            exit(1)

        rec[stage] = num_reads
        if first_count < 0:
            first_count = num_reads

        if first_count > 0:
            rec_perc[stage] = round(100 * num_reads / first_count, 1)
        else:
            rec_perc[stage] = 0

    results.append(rec)
    results_perc.append(rec_perc)

write_csv(args.results_file, results)
fn = os.path.splitext(args.results_file)
write_csv(fn[0] + '_perc' + fn[1], results_perc)
