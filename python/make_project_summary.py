from math import log
import os
import csv
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description='Concatenate alignment files from multiple samples and loci.')
    parser.add_argument('--input_dir', type=str, default=os.getcwd(),
                        help='Input directory containing sample directories (default: current directory)')
    parser.add_argument('--output_file_path', type=str, default='.',
                        help='Output file path (default: .)')
    return parser.parse_args()


args = parse_arguments()

if not os.path.isfile(os.path.join(args.output_file_path, 'pipeline_counts.csv')):
    print("Error: pipeline_counts.csv not found")
    exit(1)

if not os.path.isfile(os.path.join(args.output_file_path, 'project_alignments.csv')):
    print("Error: project_alignments.csv not found")
    exit(1)

samples = {}
project = None


with open(os.path.join(args.output_file_path, 'pipeline_counts.csv'), newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        sample = row['sample']
        locus = row['locus']
        if sample not in samples:
            samples[sample] = {}
        if locus not in samples[sample]:
            samples[sample][locus] = {
                'sample': sample,
                'locus': locus,
                'reads': row['filterSeq_quality'],
                'sequences': row['unique_atleast2'],
                'clones': row['clones'],
            }

for sample in samples.keys():
    for locus in samples[sample].keys():
        samples[sample][locus]['shannon'] = samples[sample][locus]['simpson'] = ''
        cd_path = f"{args.input_dir}/{sample}/{locus}/clones/{sample}_{locus}_clone_diversity.csv"
        if not os.path.exists(cd_path):
            print(f"Error: file {cd_path} not found")
            continue

        with open(cd_path, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row['q'] == '1':
                    samples[sample][locus]['shannon'] = round(log(float(row['d']), 10), 2)
                elif row['q'] == '2':
                    samples[sample][locus]['simpson'] = round(log(float(row['d']), 10), 2)

annots = {}
with open(os.path.join(args.output_file_path, 'project_alignments.csv'), newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        sample = row['sample']
        locus = row['locus']
        if sample not in annots:
            annots[sample] = {}
        if locus not in annots[sample]:
            annots[sample][locus] = []
        annots[sample][locus].append(row)

for sample in samples.keys():
    for locus in samples[sample].keys():
        if sample not in annots or locus not in annots[sample]:
            continue
        cdr3s = set()
        iso = {'IGG1': 0, 'IGG2': 0, 'IGG3': 0, 'IGG4': 0, 'IGA': 0, 'IGM': 0, 'IGE': 0, 'IGD': 0}
        for row in annots[sample][locus]:
            cdr3s.add(row['cdr3_aa'])
            if locus == 'IGH':
                if 'IGHG' in row['c_call']:
                    isotype = row['c_call'].split(',')[0][:5]
                elif len(row['c_call']) > 4:
                    isotype = row['c_call'].split(',')[0][:4]
                isotype = isotype.replace('IGH', 'IG')
                iso[isotype] += 1
        samples[sample][locus]['cdr3s'] = len(cdr3s)
        for isotype, count in iso.items():
            samples[sample][locus][isotype] = count

header = ['sample', 'locus', 'reads', 'sequences', 'clones', 'shannon', 'simpson', 'cdr3s', 'IGG1', 'IGG2', 'IGG3', 'IGG4', 'IGA', 'IGM', 'IGE', 'IGD']
with open(os.path.join(args.output_file_path, 'flairr_project_summary.csv'), 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=header)
    writer.writeheader()
    for sample in samples.keys():
        for locus in samples[sample].keys():
            writer.writerow(samples[sample][locus])
