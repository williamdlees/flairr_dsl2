# Plot the distribution of barcode alignments for each sample

import toml
import glob
import os
import argparse
import pandas as pd
import re
import plotly.express as px
from receptor_utils import simple_bio_seq as simple


loci = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRD', 'TRG']


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

    return sorted(list(set(samples)))


args = argparse.ArgumentParser(description='Plot barcode statistics')
args.add_argument('config_file', help='Path to the config file')
args.add_argument('output_file', help='Path to the plot file')
args = args.parse_args()

config = toml.load(args.config_file)

as_log = ''
for fd, fn in config['sample_files']:
    if fd == 'alignSets':
        as_log = fn

if as_log == '':
    print('Error: alignSets not found in config file')
    exit(1)

sample_dir = config['sample_dir']
samples = find_samples(sample_dir)

for locus in loci:
    results = []
    root = sample_dir.replace('{locus}', locus).replace('{sample}', '*')
    data = None
    for log_file in glob.glob(os.path.join(root, as_log.replace('{sample}', '*'))):
        sample_name = os.path.basename(log_file).split(as_log.split('{sample}')[1])[0]

        data = simple.read_csv(os.path.join(log_file), delimiter='\t')

        counts = {}
        for row in data:
            count = int(row['SEQCOUNT'])
            if count not in counts:
                counts[count] = 0
            counts[count] += 1

        max_count = max(counts.keys())

        for count in range(1, max_count + 1):
            if count in counts:
                results.append({'sample': sample_name, 'count': count, 'alignments': counts[count]})

    if data:
        df = pd.DataFrame(results)
        fig = px.scatter(df, x='count', y='alignments', color='sample', title='Barcode alignments in each sample',
                        log_x=True, log_y=True, labels={'alignments': 'Number of alignments of this size', 'count': 'Number of reads in alignment'})
        fn = os.path.splitext(args.output_file)
        fig.write_html(fn[0] + '_' + locus + fn[1])
