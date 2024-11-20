# Put v region IMGT sequence alignment into an annotation file

import argparse
import csv
from Bio import Align, SeqIO
from receptor_utils import simple_bio_seq as simple


def aligned_diff(novel_seq: str, ref_seq: str):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.query_right_open_gap_score = 1
    aligner.target_right_open_gap_score = 1

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

        # In the processing below, we leave the alignment spacer as -, while the IMT gap in the reference set is .
        # Once processing is complete, we change all spacers to .

        v_call = rec['v_call'].split(',')[0]
        if v_call not in germlines:
            print(f"{rec['sequence_id']}: {v_call} not found in germline reference: please provide in the reference all sequences that were used by IgBlast")
            exit(1)

        if not rec['v_sequence_start'] or not rec['v_sequence_end']:
            print(f"{rec['sequence_id']}: no v-sequence")
            rec['sequence_alignment'] = ''
            rec['germline_alignment'] = ''
            rec['v_sequence_alignment'] = ''
            continue

        v_sequence = rec['sequence'][int(rec['v_sequence_start'])-1:int(rec['v_sequence_end'])]

        if rec['d_sequence_start'] and rec['d_sequence_end']:
            d_sequence = rec['sequence'][int(rec['d_sequence_start'])-1:int(rec['d_sequence_end'])]
        else:
            d_sequence = ''

        if rec['j_sequence_start'] and rec['j_sequence_end']:
            j_sequence = rec['sequence'][int(rec['j_sequence_start'])-1:int(rec['j_sequence_end'])]
        else:
            j_sequence = ''

        if rec['c_sequence_start'] and rec['c_sequence_end']:
            c_sequence = rec['sequence'][int(rec['c_sequence_start'])-1:int(rec['c_sequence_end'])]
        else:
            c_sequence = ''
            rec['c_sequence_alignment'] = ''    # if part of the vdj sequence is missing, igblast may produce an alignment but no coordinates
            rec['c_germline_alignment'] = ''

        vdj_sequence = v_sequence + d_sequence + j_sequence + c_sequence

        ref_gapped = germlines[v_call]
        alignment = aligned_diff(ref_gapped.replace('.', ''), v_sequence)
        rec['v_cigar'] = get_cigar_from_alignment(alignment)

        # thread any insertions in the aligned_ref into the ref_gapped

        aligned_ref_ungapped = alignment._get_row(0)
        aligned_v_ungapped = alignment._get_row(1)

        # trim reference alignment to the length of the v_sequence_alignment

        while aligned_v_ungapped[-1] == '-':
            aligned_v_ungapped = aligned_v_ungapped[:-1]
            aligned_ref_ungapped = aligned_ref_ungapped[:-1]

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

        rec['sequence_alignment'] = aligned_v_gapped + rec['d_sequence_alignment'] + rec['j_sequence_alignment'] + rec['c_sequence_alignment']
        rec['germline_alignment'] = aligned_ref_gapped + rec['d_germline_alignment'] + rec['j_germline_alignment'] + rec['c_germline_alignment']
        rec['v_sequence_alignment'] = aligned_v_gapped

        # make a table to translate zero-based coords in the original sequence to zero-based coords in the aligned sequence

        seq_coords = []
        for i in range(int(rec['v_sequence_start'])):
            if i >= 0:
                seq_coords.append(-1)

        i = 0
        for b in vdj_sequence:
            while rec['sequence_alignment'][i] == '.' or rec['sequence_alignment'][i] == '-':
                i = i + 1
            seq_coords.append(i)
            i = i + 1

        def find_coord(i, coords):
            if i < len(coords):
                return coords[i]
            else:
                return coords[-1] + len(coords) - i + 1

        rec['v_alignment_start'] = 1
        rec['v_alignment_end'] = int(rec['v_alignment_start']) - 1 + len(rec['v_sequence_alignment'])
        if rec['sequence_alignment'][int(rec['v_alignment_start'])-1:rec['v_alignment_end']] != rec['v_sequence_alignment']:
            print(f"Error: {rec['sequence_id']} range of v_sequence_alignment does not match sequence_alignment")
        end = rec['v_alignment_end']
        if d_sequence:
            rec['d_alignment_start'] = rec['v_alignment_end'] + 1
            rec['d_alignment_end'] = rec['v_alignment_end'] + len(rec['d_sequence_alignment'])
            if rec['sequence_alignment'][rec['d_alignment_start']-1:rec['d_alignment_end']] != rec['d_sequence_alignment']:
                print(f"Error: {rec['sequence_id']} range of d_sequence_alignment does not match sequence_alignment")
            end = rec['d_alignment_end']
        else:
            rec['d_alignment_start'] = 0
            rec['d_alignment_end'] = 0
        if j_sequence:
            rec['j_alignment_start'] = end + 1
            rec['j_alignment_end'] = end + len(rec['j_sequence_alignment'])
            if rec['sequence_alignment'][rec['j_alignment_start']-1:rec['j_alignment_end']] != rec['j_sequence_alignment']:
                print(f"Error: {rec['sequence_id']} range of j_sequence_alignment does not match sequence_alignment")
            end = rec['j_alignment_end']
        else:
            rec['j_alignment_start'] = 0
            rec['j_alignment_end'] = 0
        if c_sequence:
            rec['c_alignment_start'] = end + 1
            rec['c_alignment_end'] = end + len(rec['c_sequence_alignment'])
            if rec['sequence_alignment'][rec['c_alignment_start']-1:rec['c_alignment_end']] != rec['c_sequence_alignment']:
                print(f"Error: {rec['sequence_id']} range of c_sequence_alignment does not match sequence_alignment")

        for field in ('cdr1', 'cdr2'):
            if rec[field]:
                rec[f'{field}_aligned_start'] = find_coord(int(rec[f'{field}_start']), seq_coords) + 1
                rec[f'{field}_aligned_end'] = find_coord(int(rec[f'{field}_end']), seq_coords) + 1
                cdr_seq = rec['sequence_alignment'][rec[f'{field}_aligned_start']-1:rec[f'{field}_aligned_end']]
                if cdr_seq.replace('.', '').replace('-', '') != rec[field]:
                    print(f"Error: {rec['sequence_id']} {field} changed from {rec[field]} to non-matching {cdr_seq}")
                rec[field] = cdr_seq
        
        rec['sequence_alignment'] = rec['sequence_alignment'].replace('-', '.')
        rec['germline_alignment'] = rec['germline_alignment'].replace('-', '.')
        rec['v_sequence_alignment'] = rec['v_sequence_alignment'].replace('-', '.')

        # sanity checks - non gapped characters should match

        if rec['sequence_alignment'].replace('.', '') != vdj_sequence:
            print(f"Error: {rec['sequence_id']} revised sequence_alignment does not match original")
        
        if len(rec['germline_alignment']) != len(rec['sequence_alignment']):
            print(f"Error: {rec['sequence_id']} revised germline_alignment and sequence alignments have different lengths")
        
        if rec['v_sequence_alignment'].replace('.', '') != v_sequence:
            print(f"Error: {rec['sequence_id']} revised v_sequence_alignment does not match original")

        writer.writerow(rec)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Input annotation file')
    parser.add_argument('germline_ref', help='Germline reference file')
    parser.add_argument('outfile', help='Output annotation file')
    args = parser.parse_args()
    main(args.infile, args.outfile, args.germline_ref)