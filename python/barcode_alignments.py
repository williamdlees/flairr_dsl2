# Plot the distribution of barcode alignments for each sample

import toml
import glob
import os
import argparse
import pandas as pd
import csv
import plotly.express as px


# Read csv file into a list of dicts
def read_csv(file: str, delimiter: str = None):
    """Read a delimited file into a list of dicts (as produced by DictReader)
    """
    ret = []
    with open(file, 'r') as fi:
        if delimiter:
            reader = csv.DictReader(fi, delimiter=delimiter)
        else:
            reader = csv.DictReader(fi)
        for row in reader:
            ret.append(row)

    return ret


args = argparse.ArgumentParser(description='Summarise read counts at each pipeline stage for each sample')
args.add_argument('config_file', help='Path to the config file')
args.add_argument('output_file', help='Path to the plot file')
args = args.parse_args()

config = toml.load(args.config_file)

results = []

as_log = ''
for fd, fn in config['sample_files']:
    if fd == 'alignSets':
        as_log = fn

if as_log == '':
    print('Error: alignSets not found in config file')
    exit(1)

for log_file in glob.glob(os.path.join(config['sample_dir'], as_log.replace('{sample}', '*'))):
    sample_name = os.path.basename(log_file).split(as_log.split('{sample}')[1])[0]

    data = read_csv(os.path.join(log_file), delimiter='\t')

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

df = pd.DataFrame(results)
fig = px.scatter(df, x='count', y='alignments', color='sample', title='Barcode alignments in each sample',
                 log_x=True, log_y=True, labels={'alignments': 'Number of alignments of this size', 'count': 'Number of reads in alignment'})
fig.write_html(args.output_file)
print("")