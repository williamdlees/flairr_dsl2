# Create a table summarising read counts at each pipeline stage for each sample

import toml
import glob
import os
import argparse
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


def find_sample_roots(template):
    sample_roots = []
    template = template.replace('\\', '/')
    locus_pattern = '|'.join(re.escape(locus) for locus in loci)
    pattern = template.replace('{sample}', r'(?P<sample>[^/]+)')
    pattern = pattern.replace(
        '{locus}',
        rf'(?P<locus_dir>(?P<locus>{locus_pattern})(?:_(?P<run>[^/]+))?)'
    )
    pattern = re.compile(rf'^{pattern}$')
    file_template = template.replace('{sample}', '*').replace('{locus}', '*')

    for path_name in glob.glob(file_template):
        normalized_path = path_name.replace('\\', '/')
        match = pattern.match(normalized_path)
        if not match:
            continue

        sample_roots.append({
            'sample': match.group('sample'),
            'locus': match.group('locus'),
            'locus_dir': match.group('locus_dir'),
            'run': match.group('run') or '',
            'root': path_name,
        })

    unique_roots = {}
    for sample_root in sample_roots:
        key = (sample_root['sample'], sample_root['locus_dir'], sample_root['root'])
        unique_roots[key] = sample_root

    return [
        unique_roots[key]
        for key in sorted(unique_roots, key=lambda item: (item[0], item[1], item[2]))
    ]


args = argparse.ArgumentParser(description='Summarise read counts at each pipeline stage for each sample')
args.add_argument('config_file', help='Path to the config file')
args.add_argument('results_file', help='Path to the results file')
args = args.parse_args()

config = toml.load(args.config_file)

results = []
results_perc = []

sample_dir = config['sample_dir']
print('Processing samples from directory ', sample_dir)
sample_roots = find_sample_roots(sample_dir)

if not sample_roots:
    print("No samples found.")
    exit(0)

for sample_root in sample_roots:
    sample_name = sample_root['sample']

    locus = sample_root['locus']
    locus_dir = sample_root['locus_dir']
    run = sample_root['run']
    root = sample_root['root']

    print(f"{sample_name} {locus}")

    rec = {'locus': locus, 'run': run, 'sample': sample_name}
    rec_perc = {'locus': locus, 'run': run, 'sample': sample_name}
    first_count = -1
    found_data_for_sample = False

    for stage, stage_filespec in config['sample_files']:
        stage_pattern_template = os.path.join(root, stage_filespec).replace('{sample}', sample_name)
        stage_patterns = [stage_pattern_template.replace('{locus}', locus_dir)]

        if locus_dir != locus:
            stage_patterns.append(stage_pattern_template.replace('{locus}', locus))

        stage_filename_glob = []
        for stage_pattern in stage_patterns:
            for matched_file in glob.glob(stage_pattern):
                if matched_file not in stage_filename_glob:
                    stage_filename_glob.append(matched_file)

        if len(stage_filename_glob) == 0:
            rec[stage] = ''
            rec_perc[stage] = ''
            continue
        elif len(stage_filename_glob) > 1:
            print(f'Error: multiple files match {sample_name} {locus_dir} {stage}')
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

# sort results by sample, locus, run
results.sort(key=lambda x: (x['sample'], x['locus'], x['run']))
results_perc.sort(key=lambda x: (x['sample'], x['locus'], x['run']))


def _finalise_sample_recs(sample_recs, perc_recs, compacted_results, perc_compacted_results):
    if len(sample_recs) == 2 and sample_recs[0]['run'] == '' and sample_recs[1]['run'] != '':
        print('merging')
        for k in sample_recs[1].keys():
            if sample_recs[1][k] != '' and sample_recs[0][k] == '':
                sample_recs[0][k] = sample_recs[1][k]
        for k in perc_recs[1].keys():
            if perc_recs[1][k] != '' and perc_recs[0][k] == '':
                perc_recs[0][k] = perc_recs[1][k]

        if sample_recs[1]['filterSeq_quality']:
            try:
                total_reads = int(sample_recs[0]['filterSeq_quality'])
            except ValueError:
                total_reads = 0
                
        for k in sample_recs[0].keys():
            if k not in ['sample', 'locus', 'run'] and sample_recs[0][k]:
                try:
                    perc_recs[0][k] = round(100 * int(sample_recs[0][k]) / total_reads, 1)
                except ValueError:
                    perc_recs[0][k] = ''
        compacted_results.append(sample_recs[0])
        perc_compacted_results.append(perc_recs[0])
    else:
        print('extending')
        if sample_recs[0]['run'] == '':
            sample_recs[0]['run'] = 'combined'

        total_reads = 0
        for i in range(1, len(sample_recs)):
            if sample_recs[i]['filterSeq_quality']:
                try:
                    total_reads += int(sample_recs[i]['filterSeq_quality'])
                except ValueError:
                    pass

            for k in sample_recs[0].keys():
                if k not in ['sample', 'locus', 'run'] and sample_recs[0][k]:
                    try:
                        perc_recs[0][k] = round(100 * int(sample_recs[0][k]) / total_reads, 1)
                    except ValueError:
                        perc_recs[0][k] = ''

        compacted_results.extend(sample_recs)
        perc_compacted_results.extend(perc_recs)

    return compacted_results, perc_compacted_results


def compact_sample_results(records, perc_records):
    # if there are just 2 records for a sample and one is an 'annotation' record which has no run, and the other is a 'run' record which has a run,
    # then merge the two records into one
    compacted_results = []
    perc_compacted_results = []
    locus_sample = ''
    sample_recs = []
    perc_recs = []

    for rec, perc_rec in zip(records, perc_records):
        rec_locus_sample = rec['locus'] + '_' + rec['sample']
        if rec_locus_sample != locus_sample:
            if sample_recs:
                print(f"Processing {locus_sample} with {len(sample_recs)} records")
                compacted_results, perc_compacted_results = _finalise_sample_recs(sample_recs, perc_recs, compacted_results, perc_compacted_results)
            locus_sample = rec_locus_sample
            sample_recs = [rec]
            perc_recs = [perc_rec]
        else:
            sample_recs.append(rec)
            perc_recs.append(perc_rec)

    if sample_recs:
        print(f"Processing {locus_sample} with {len(sample_recs)} records")
        compacted_results, perc_compacted_results = _finalise_sample_recs(sample_recs, perc_recs, compacted_results, perc_compacted_results)

    return compacted_results, perc_compacted_results


compacted_results, perc_compacted_results = compact_sample_results(results, results_perc)

simple.write_csv(args.results_file, compacted_results)
fn = os.path.splitext(args.results_file)
simple.write_csv(fn[0] + '_perc' + fn[1], perc_compacted_results)
