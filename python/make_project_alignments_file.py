import os
import csv

def scan_directories():
    current_dir = os.getcwd()
    samples = {}

    for sample_name in os.listdir(current_dir):
        sample_path = os.path.join(current_dir, sample_name)
        if os.path.isdir(sample_path):
            loci = []
            for locus_name in os.listdir(sample_path):
                locus_path = os.path.join(sample_path, locus_name)
                if os.path.isdir(locus_path):
                    loci.append(locus_name)
            samples[sample_name] = loci

    return samples

def find_alignment_file(sample, locus):
    alignment_dir = os.path.join(os.getcwd(), sample, locus, 'alignment')
    if os.path.isdir(alignment_dir):
        for file_name in os.listdir(alignment_dir):
            if file_name.endswith('align_non-personalized.tsv'):
                return os.path.join(alignment_dir, file_name)
    return None

def concatenate_files(samples):
    all_data = []
    columns = set()

    for sample, loci in samples.items():
        print(f"Processing sample '{sample}' with loci: {', '.join(loci)}")
        for locus in loci:
            alignment_file = find_alignment_file(sample, locus)
            if alignment_file:
                with open(alignment_file, newline='') as csvfile:
                    reader = csv.DictReader(csvfile, delimiter='\t')
                    for row in reader:
                        row['Sample'] = sample
                        row['Locus'] = locus
                        all_data.append(row)
                        columns.update(row.keys())
            else:
                print(f"Warning: No alignment file found for sample '{sample}', locus '{locus}'")

    if all_data:
        with open('project_alignments.csv', 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=columns, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(all_data)
    else:
        print("No alignment files found to concatenate.")

if __name__ == "__main__":
    samples = scan_directories()
    concatenate_files(samples)