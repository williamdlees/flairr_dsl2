import os
import csv
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description='Concatenate alignment files from multiple samples and loci.')
    parser.add_argument('--input_dir', type=str, default=os.getcwd(),
                        help='Input directory containing sample directories (default: current directory)')
    parser.add_argument('--output_file', type=str, default='project_alignments.csv',
                        help='Output file path (default: project_alignments.csv)')
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


def find_alignment_file(input_dir, sample, locus):
    alignment_dir = os.path.join(input_dir, sample, locus, 'alignment')
    if os.path.isdir(alignment_dir):
        for file_name in os.listdir(alignment_dir):
            if file_name.endswith('align_non-personalized.tsv'):
                return os.path.join(alignment_dir, file_name)
    return None


def concatenate_files(samples, input_dir, output_file):
    all_data = []
    header = None

    for sample, loci in samples.items():
        print(f"Processing sample '{sample}' with loci: {', '.join(loci)}")
        for locus in loci:
            alignment_file = find_alignment_file(input_dir, sample, locus)
            if alignment_file:
                with open(alignment_file, newline='') as csvfile:
                    reader = csv.DictReader(csvfile, delimiter='\t')
                    for row in reader:
                        if header is None:
                            # Preserve the original column order from the first file
                            header = ['sample', 'locus']
                            header.extend(list(row.keys()))
                        row['sample'] = sample
                        row['locus'] = locus
                        all_data.append(row)
            else:
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
    concatenate_files(samples, args.input_dir, args.output_file)
