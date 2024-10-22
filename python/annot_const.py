# Put constant region sequence and identity into an annotation file
# The constant region sequence is extracted from the input sequence on each row

import argparse
import csv
from Bio import Align, SeqIO
from receptor_utils import simple_bio_seq as simple


def aligned_diff(novel_seq: str, ref_seq: str):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.target_end_gap_score = 0

    return aligner.align(novel_seq, ref_seq)[0]


def get_cigar_from_alignment(alignment):
    target_indices, query_indices = alignment.aligned
    cigar = []

    last_t_end, last_q_end = 0, 0

    for (t_start, t_end), (q_start, q_end) in zip(target_indices, query_indices):
        if t_start > last_t_end:
            cigar.append(f'{t_start - last_t_end}D')
        if q_start > last_q_end:
            cigar.append(f'{q_start - last_q_end}I')

        length = t_end - t_start
        cigar.append(f'{length}M')

        last_t_end, last_q_end = t_end, q_end

    return ''.join(cigar)


def main(infile, outfile, germline_ref):
    germlines = simple.read_fasta(germline_ref)
    fi = open(infile, 'r')
    reader = csv.DictReader(fi, delimiter='\t')
    fo = open(outfile, 'w', newline='')
    writer = csv.DictWriter(fo, reader.fieldnames + ['const', 'const_in_frame', 'const_ref', 'c_cigar', 'c_score', 'c_identity'], delimiter='\t')
    writer.writeheader()

    for rec in reader:
        if not rec['c_call'] or not rec['j_sequence_end']:
            continue

        rec['c_call'] = rec['c_call'].split(',')[0]
        if rec['c_call'] not in germlines:
            print(f"Germline {rec['c_call']} not found in the germline file")
            continue

        try:
            j_sequence_end = int(rec['j_sequence_end'])
        except ValueError:
            continue

        # allow a window around the J end to find the constant region

        rough_seq = rec['sequence'][j_sequence_end - 4:]

        if len(rough_seq) < 10:
            continue

        best_alignment = aligned_diff(germlines[rec['c_call']], rough_seq)
        c_sequence_start = j_sequence_end - 3 + best_alignment.aligned[1][0][0]
        best_alignment = aligned_diff(germlines[rec['c_call']], rec['sequence'][c_sequence_start - 1:])

        rec['const'] = best_alignment._get_row(1)
        rec['const_in_frame'] = 'F' if '*' in simple.translate(best_alignment._get_row(1).replace('-', '')) else 'T'
        rec['const_ref'] = best_alignment._get_row(0)
        rec['c_cigar'] = get_cigar_from_alignment(best_alignment)
        rec['c_score'] = best_alignment.counts().identities
        rec['c_identity'] = round(100 * best_alignment.counts().identities/len(best_alignment[0]), 2)

        writer.writerow(rec)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Input annotation file')
    parser.add_argument('germline_ref', help='Germline reference file')
    parser.add_argument('outfile', help='Output annotation file')
    args = parser.parse_args()
    main(args.infile, args.outfile, args.germline_ref)