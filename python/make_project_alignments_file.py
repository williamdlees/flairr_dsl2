import os
import csv
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description='Concatenate alignment files from multiple samples and loci.')
    parser.add_argument('--input_dir', type=str, default=os.getcwd(),
                        help='Input directory containing sample directories (default: current directory)')
    parser.add_argument('--output_file', type=str, default='project_alignments.csv',
                        help='Output file path (default: project_alignments.csv)')
    parser.add_argument('--personalized', action='store_true',
                        help='Whether to use personalized alignment files (default: False)')
    return parser.parse_args()


def scan_directories(input_dir):
    samples = {}

    for sample_name in os.listdir(input_dir):
        sample_path = os.path.join(input_dir, sample_name)
        if os.path.isdir(sample_path):
            loci = []
            for locus_name in os.listdir(sample_path):
                locus_path = os.path.join(sample_path, locus_name)
                if os.path.isdir(locus_path):
                    loci.append(locus_name)
            samples[sample_name] = loci

    return samples


def find_alignment_file(input_dir, sample, locus, personalized):
    alignment_dir = os.path.join(input_dir, sample, locus, 'alignment')
    target_suffix = 'align_personalized.tsv' if personalized else 'align_non-personalized.tsv'
    if os.path.isdir(alignment_dir):
        for file_name in os.listdir(alignment_dir):
            if file_name.endswith(target_suffix):
                return os.path.join(alignment_dir, file_name)
    return None


def find_clone_ids(input_dir, sample, locus):
    clone_file = None
    clone_ids = {}
    clone_dir = os.path.join(input_dir, sample, locus, 'clones')
    if os.path.isdir(clone_dir):
        for file_name in os.listdir(clone_dir):
            if file_name.endswith('_clone_pass.tsv'):
                clone_file = file_name
                break
    if clone_file:
        clone_file_path = os.path.join(clone_dir, clone_file)
        print(f"Found clone file for sample '{sample}', locus '{locus}': {clone_file_path}")
        with open(clone_file_path, newline='') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            for row in reader:
                clone_id = row.get('clone_id')
                if clone_id:
                    clone_ids[row['sequence_id']] = clone_id
    return clone_ids


def concatenate_files(samples, input_dir, output_file, personalized):
    all_data = []
    header = None

    for sample, loci in samples.items():
        print(f"Processing sample '{sample}' with loci: {', '.join(loci)}")
        for locus in loci:
            alignment_file = find_alignment_file(input_dir, sample, locus, personalized)
            if alignment_file:
                clone_ids = {}
                if personalized and 'IG' in locus:
                    clone_ids = find_clone_ids(input_dir, sample, locus)
                with open(alignment_file, newline='') as csvfile:
                    reader = csv.DictReader(csvfile, delimiter='\t')
                    for row in reader:
                        if header is None:
                            # Preserve the original column order from the first file
                            header = ['sample', 'locus']
                            header.extend(list(row.keys()))
                            header.extend(['clone_id', 'isotype_call'])
                        row['sample'] = sample
                        row['locus'] = locus
                        row['clone_id'] = clone_ids.get(row['sequence_id'], '')
                        row['isotype_call'] = ''

                        if 'c_call' in row and 'IGH' in row['c_call']:
                            row['isotype_call'] = row['c_call'].split(',')[0].split('*')[0]

                        all_data.append(row)
            else:
                if '_' not in locus:
                    print(f"Warning: No alignment file found for sample '{sample}', locus '{locus}'")
                continue

    if all_data:
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=header, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(all_data)
        print(f"Output written to '{output_file}'")
    else:
        print("No alignment files found to concatenate.")


if __name__ == "__main__":
    args = parse_arguments()
    samples = scan_directories(args.input_dir)
    k = list(samples.keys())[0]
    samples = {k: samples[k]}  # first sample only for testing
    concatenate_files(samples, args.input_dir, args.output_file, args.personalized)
