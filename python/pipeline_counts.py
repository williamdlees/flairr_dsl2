# Create a table summarising read counts at each pipeline stage for each sample

import toml
import glob
import os
import argparse
import csv
import re
from receptor_utils import simple_bio_seq as simple


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


loci = ['IGH', 'IGK', 'IGL', 'IGM', 'IGG', 'IGA', 'IGD', 'IGE', 'TRA', 'TRB', 'TRD', 'TRG']


def find_samples(template):
    samples = []
    
    # Walk through directory structure
    for locus in loci:
        pattern = template
        pattern = pattern.replace("{sample}", "(?P<sample>[^/]+)").replace("{locus}", locus)
        pattern = re.compile(pattern)
        file_template = sample_dir.replace('\\', '/').replace('{locus}', locus).replace('{sample}', '*')
        for p in glob.glob(file_template):
            p = p.replace('\\', '/')
            match = pattern.match(p)
            if match:
                sample = match.group("sample")
                samples.append(sample)
                # print(sample)

    return sorted(list(set(samples)))


args = argparse.ArgumentParser(description='Summarise read counts at each pipeline stage for each sample')
args.add_argument('config_file', help='Path to the config file')
args.add_argument('results_file', help='Path to the results file')
args = args.parse_args()

config = toml.load(args.config_file)

results = []
results_perc = []

sample_dir = config['sample_dir']
print('Processing samples from directory ', sample_dir)
samples = find_samples(sample_dir)

if not samples:
    print("No samples found.")
    exit(0)

for locus in loci:
    for sample_name in samples:
        root = sample_dir.replace('{locus}', locus).replace('{sample}', sample_name)
    
        rec = {'locus': locus, 'sample': sample_name}
        rec_perc = {'locus': locus, 'sample': sample_name}
        first_count = -1
        found_data_for_sample = False

        for stage, stage_filespec in config['sample_files']:
            stage_filename_glob = glob.glob(os.path.join(root, stage_filespec).replace('{sample}', sample_name).replace('{locus}', locus))

            if len(stage_filename_glob) == 0:
                rec[stage] = 'NA'
                rec_perc[stage] = 'NA'
                continue
            elif len(stage_filename_glob) > 1:
                print(f'Error: multiple files match {sample_name} {stage}')
                exit(1)

            found_data_for_sample = True
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

        if found_data_for_sample:
            results.append(rec)
            results_perc.append(rec_perc)

simple.write_csv(args.results_file, results)
fn = os.path.splitext(args.results_file)
simple.write_csv(fn[0] + '_perc' + fn[1], results_perc)
