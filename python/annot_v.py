# Put v region IMGT sequence alignment into an annotation file

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
    fieldnames = reader.fieldnames
    fieldnames.extend(['cdr1_aligned_end', 'cdr1_aligned_start', 'cdr2_aligned_start', 'cdr2_aligned_end'])
    writer = csv.DictWriter(fo, fieldnames, delimiter='\t')
    writer.writeheader()

    for rec in reader:
        if not rec['v_call']:
            continue

        # In the processing below, we rely on the fact that the alignment spacer is -, while the IMT gap is .
        # Once processing is complete, we change all spacers to .

        if rec['sequence_id'] == '1TTTTCTTATGTTAGATCTTGGG|PRCONS=IGG_CH3_PCR-primer|PRFREQ=1.0|CONSCOUNT=2|DUPCOUNT=1':
            print('debug')

        v_call = rec['v_call'].split(',')[0]
        if v_call not in germlines:
            print(f"{rec['sequence_id']}: {v_call} not found in germline reference: please provide in the reference all sequences that were used by IgBlast")
            exit(1)

        if rec['v_sequence_alignment'] not in rec['sequence_alignment']:
            print(f"{rec['sequence_id']}: v_sequence_alignment not found in sequence_alignment: skipping")
            continue

        if len(rec['v_sequence_alignment']) != len(rec['v_germline_alignment']):
            print(f"{rec['sequence_id']}: v_sequence_alignment and v_germline_alignment have different lengths: skipping")
            continue

        ref_gapped = germlines[v_call]
        alignment = aligned_diff(ref_gapped.replace('.', ''), rec['v_sequence_alignment'].replace('.', '').replace('-', ''))
        rec['v_cigar'] = get_cigar_from_alignment(alignment)

        # thread any insertions in the aligned_ref into the ref_gapped

        aligned_ref_ungapped = alignment._get_row(0)
        aligned_ref_ungapped_i = iter(aligned_ref_ungapped)
        aligned_ref_gapped = ''
        for b in ref_gapped:
            if b == '.':
                aligned_ref_gapped += '.'
            else:
                try:
                    aligned_ref_gapped += next(aligned_ref_ungapped_i)
                except StopIteration:
                    break

        while True:
            try:
                aligned_ref_gapped += next(aligned_ref_ungapped_i)
            except StopIteration:
                break

        aligned_v_ungapped = alignment._get_row(1)
        aligned_v_ungapped_i = iter(aligned_v_ungapped)
        aligned_v_gapped = ''
        for b in aligned_ref_gapped:
            if b == '.':
                aligned_v_gapped += '.'
            else:
                aligned_v_gapped += next(aligned_v_ungapped_i)

        while True:
            try:
                aligned_v_gapped += next(aligned_v_ungapped_i)
            except StopIteration:
                break

        original_sequence_alignment = rec['sequence_alignment']
        original_germline_alignment = rec['germline_alignment']
        original_v_sequence_alignment = rec['v_sequence_alignment']
        rec['sequence_alignment'] = aligned_v_gapped + rec['sequence_alignment'][len(rec['v_sequence_alignment']):]
        rec['germline_alignment'] = aligned_ref_gapped + rec['germline_alignment'][len(rec['v_sequence_alignment']):]
        rec['v_sequence_alignment'] = aligned_v_gapped

        # make a table to translate zero-based coords in the original sequence to zero-based coords in the aligned sequence

        seq_coords = []
        i = 0
        for b in original_sequence_alignment:
            while rec['sequence_alignment'][i] == '.':
                i = i + 1
            seq_coords.append(i)
            i = i + 1

        germline_coords = []
        i = 0
        for b in original_germline_alignment:
            while rec['germline_alignment'][i] == '.':
                i = i + 1
            germline_coords.append(i)
            i = i + 1

        def find_coord(i, coords):
            if i < len(coords):
                return coords[i]
            else:
                return coords[-1] + len(coords) - i + 1

        seq_fields = [
            'v_alignment_start', 'v_alignment_end',
            'd_alignment_start', 'd_alignment_end',
            'j_alignment_start', 'j_alignment_end',
        ]

        for field in seq_fields:
            if rec[field] and rec[field].isdigit():
                rec[field] = find_coord(int(rec[field])-1, seq_coords) + 1

        for field in ('cdr1', 'cdr2'):
            unaligned_pos = original_sequence_alignment.index(rec[field])
            if unaligned_pos < 0:
                print(f"{rec['sequence_id']} {field} not found in sequence_alignment")
                continue

            rec[field + '_aligned_start'] = find_coord(unaligned_pos, seq_coords) + 1
            rec[field + '_aligned_end'] = find_coord(unaligned_pos + len(rec[field]), seq_coords)
            aligned_seq = rec['sequence_alignment'][rec[field + '_aligned_start']-1:rec[field + '_aligned_end']]
            if aligned_seq.replace('.', '') != rec[field]:
                print(f"Error: {rec['sequence_id']} {field} changed from {rec[field]} to non-matching {aligned_seq}")
            rec[field] = aligned_seq

        rec['sequence_alignment'] = rec['sequence_alignment'].replace('-', '.')
        rec['germline_alignment'] = rec['germline_alignment'].replace('-', '.')
        rec['v_sequence_alignment'] = rec['v_sequence_alignment'].replace('-', '.')

        # sanity checks - non gapped characters should match

        if rec['sequence_alignment'].replace('.', '') != original_sequence_alignment.replace('.', '').replace('-', ''):
            print(f"Error: {rec['sequence_id']} revised sequence_alignment does not match original")
        
        if rec['germline_alignment'].replace('.', '') != original_germline_alignment.replace('.', '').replace('-', ''):
            print(f"Error: {rec['sequence_id']} revised germline_alignment does not match original")
        
        if rec['v_sequence_alignment'].replace('.', '') != original_v_sequence_alignment.replace('.', '').replace('-', ''):
            print(f"Error: {rec['sequence_id']} revised v_sequence_alignment does not match original")

        writer.writerow(rec)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Input annotation file')
    parser.add_argument('germline_ref', help='Germline reference file')
    parser.add_argument('outfile', help='Output annotation file')
    args = parser.parse_args()
    main(args.infile, args.outfile, args.germline_ref)