# Align the V-region in an AIRR format arrangement file with an "IMGT gapped" V reference set

# Note - the constant region is not included in the sequence alignment. The MiAIRR spec says "Typically, this will include only the V(D)J region, but that is not a requirement."
# however the Immcantation script CreateGermlines.py errors out if the constant region is included.

import argparse
import csv
from Bio import Align
from receptor_utils import simple_bio_seq as simple

csv.field_size_limit(1024*1024*1024)

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


class report_errors:
    def __init__(self, logfile, echo_errors):
        self.logfile = logfile
        self.log = open(logfile, 'w')
        self.echo_errors = echo_errors
        self.error_count = 0
        self.warning_count = 0

    def report(self, message):
        with open(self.logfile, 'a') as log:
            log.write(message + '\n')
        if self.echo_errors:
            print(message)
        if 'Error' in message:
            self.error_count += 1
        if 'Warning' in message:
            self.warning_count += 1


def main(infile, outfile, germline_ref, logfile, condense_errors, echo_errors, mask_germline_np):
    germlines = simple.read_fasta(germline_ref)
    fi = open(infile, 'r')
    reader = csv.DictReader(fi, delimiter='\t')
    fo = open(outfile, 'w', newline='')
    fieldnames = reader.fieldnames
    fieldnames.extend(['cdr1_aligned_end', 'cdr1_aligned_start', 'cdr2_aligned_start', 'cdr2_aligned_end', 'fwr3_aligned_end', 'fwr4_aligned_start', 'fwr2_aligned_start', 'fwr4_aligned_end', 'fwr1_aligned_end', 'fwr2_aligned_end', 'fwr1_aligned_start', 'fwr3_aligned_start', 'consensus_count', 'duplicate_count'])
    writer = csv.DictWriter(fo, fieldnames, delimiter='\t')
    writer.writeheader()
    reported_fields = []
    reporter = report_errors(logfile, echo_errors)
    rowcount = 0

    for rec in reader:
        rowcount += 1

        for h_c in rec['sequence_id'].split('|'):
            if 'DUPCOUNT=' in h_c:
                rec['duplicate_count'] = h_c.split('=')[1]
            elif 'CONSCOUNT=' in h_c:
                rec['consensus_count'] = h_c.split('=')[1]

        for f in ['stop_codon', 'vj_in_frame', 'v_frameshift', 'productive']:
            if rec[f] == '':
                rec[f] = 'F'

        if not rec['v_call']:
            reporter.report(f"Warning: {rec['sequence_id']}: no v_call. No V-alignment added.")
            continue

        # In the processing below, we leave the alignment spacer as -, while the IMT gap in the reference set is .
        # Once processing is complete, we change all spacers to .

        v_call = rec['v_call'].split(',')[0]
        if v_call not in germlines:
            reporter.report(f"{rec['sequence_id']}: {v_call} not found in germline reference: please provide in the reference all sequences that were used by IgBlast")
            exit(1)

        if not rec['v_sequence_start'] or not rec['v_sequence_end']:
            reporter.report(f"{rec['sequence_id']}: no v-sequence")
            rec['sequence_alignment'] = ''
            rec['germline_alignment'] = ''
            rec['v_sequence_alignment'] = ''
            continue

        v_sequence = rec['sequence'][int(rec['v_sequence_start'])-1:int(rec['v_sequence_end'])]
        end = int(rec['v_sequence_end'])

        if rec['d_sequence_start'] and rec['d_sequence_end']:
            d_sequence = rec['sequence'][int(rec['d_sequence_start'])-1:int(rec['d_sequence_end'])]
            end = int(rec['d_sequence_end'])
        else:
            d_sequence = ''

        if rec['j_sequence_start'] and rec['j_sequence_end']:
            j_sequence = rec['sequence'][int(rec['j_sequence_start'])-1:int(rec['j_sequence_end'])]
            end = int(rec['j_sequence_end'])
        else:
            j_sequence = ''

        # a mis-alignment at the 3' end of the j sequence can cause the j_alignment to be truncated
        fwr4_extra = ''
        if rec['fwr4_end'] and int(rec['fwr4_end']) > end:
            fwr4_extra = rec['sequence'][end:int(rec['fwr4_end'])]
            end = int(rec['fwr4_end'])

        vdj_sequence = rec['sequence'][int(rec['v_sequence_start'])-1:end]

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

        rec['sequence_alignment']
        rec['sequence_alignment'] = aligned_v_gapped + rec['np1'] + rec['d_sequence_alignment'] + rec['np2'] + rec['j_sequence_alignment'] + fwr4_extra

        if mask_germline_np:
            rec['germline_alignment'] = aligned_ref_gapped + 'N'*len(rec['np1']) + rec['d_germline_alignment'] + 'N'*len(rec['np2']) + rec['j_germline_alignment'] + fwr4_extra
        else:
            rec['germline_alignment'] = aligned_ref_gapped + rec['np1'] + rec['d_germline_alignment'] + rec['np2'] + rec['j_germline_alignment'] + fwr4_extra
            
        rec['v_sequence_alignment'] = aligned_v_gapped
        rec['v_germline_alignment'] = aligned_ref_gapped

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
            
        # adjust alignment coords to agree with the new alignment

        rec['v_alignment_start'] = 1
        old_alignment_end = int(rec['v_alignment_end'])
        rec['v_alignment_end'] = int(rec['v_alignment_start']) - 1 + len(rec['v_sequence_alignment'])
        if rec['sequence_alignment'][int(rec['v_alignment_start'])-1:rec['v_alignment_end']] != rec['v_sequence_alignment']:
            reporter.report(f"Error: {rec['sequence_id']} range of v_sequence_alignment does not match sequence_alignment")
        offset = rec['v_alignment_end'] - old_alignment_end

        rec['v_germline_start'] = 1
        rec['v_germline_end'] = int(rec['v_alignment_end'])

        if d_sequence:
            rec['d_alignment_start'] = int(rec['d_alignment_start']) + offset
            rec['d_alignment_end'] = int(rec['d_alignment_end']) + offset
            if rec['sequence_alignment'][rec['d_alignment_start']-1:rec['d_alignment_end']] != rec['d_sequence_alignment']:
                reporter.report(f"Error: {rec['sequence_id']} range of d_sequence_alignment does not match sequence_alignment")
        else:
            rec['d_alignment_start'] = 0
            rec['d_alignment_end'] = 0
        if j_sequence:
            rec['j_alignment_start'] = int(rec['j_alignment_start']) + offset
            rec['j_alignment_end'] = int(rec['j_alignment_end']) + offset
            if rec['sequence_alignment'][rec['j_alignment_start']-1:rec['j_alignment_end']] != rec['j_sequence_alignment']:
                reporter.report(f"Error: {rec['sequence_id']} range of j_sequence_alignment does not match sequence_alignment")
        else:
            rec['j_alignment_start'] = 0
            rec['j_alignment_end'] = 0

        for field in ('cdr1', 'cdr2', 'fwr1', 'fwr2', 'fwr3', 'fwr4'):
            if rec[field]:
                rec[f'{field}_aligned_start'] = find_coord(int(rec[f'{field}_start']), seq_coords) + 1
                rec[f'{field}_aligned_end'] = find_coord(int(rec[f'{field}_end']), seq_coords) + 1
                cdr_seq = rec['sequence_alignment'][rec[f'{field}_aligned_start']-1:rec[f'{field}_aligned_end']]
                if cdr_seq.replace('.', '').replace('-', '') != rec[field]:
                    if int(rec[f'{field}_start']) < int(rec['v_sequence_start']) or int(rec[f'{field}_end']) > int(rec['v_sequence_end']):
                        if rec[field] not in reported_fields:
                            reporter.report(f"Warning: IgBlast {rec['sequence_id']} {field} extended outside the V region (now corrected)")
                            reporter.report(f"v_sequence_start: {rec['v_sequence_start']}, d_sequence_start: {rec['d_sequence_start']}, {field}_start: {rec[f'{field}_start']}, {field}_end: {rec[f'{field}_end']}")
                            reporter.report(f"v_sequence_alignment:\n{rec['v_sequence_alignment']}")
                            if condense_errors:
                                reported_fields.append(rec[field])
                    else:
                        reporter.report(f"Error: {rec['sequence_id']} {field} changed from {rec[field]} to non-matching {cdr_seq}")
                rec[field] = cdr_seq
                rec[field + '_aa'] = simple.translate(cdr_seq)
        
        rec['sequence_alignment'] = rec['sequence_alignment'].replace('-', '.')
        rec['germline_alignment'] = rec['germline_alignment'].replace('-', '.')
        rec['v_sequence_alignment'] = rec['v_sequence_alignment'].replace('-', '.')
        rec['v_germline_alignment'] = rec['v_germline_alignment'].replace('-', '.')

        rec['sequence_alignment_aa'] = simple.translate(rec['sequence_alignment'])
        rec['germline_alignment_aa'] = simple.translate(rec['germline_alignment'])
        rec['v_sequence_alignment_aa'] = simple.translate(rec['v_sequence_alignment'])
        rec['v_germline_alignment_aa'] = simple.translate(rec['v_germline_alignment'])

        # sanity checks - non gapped characters should match

        if rec['sequence_alignment'].replace('.', '') != vdj_sequence:
            reporter.report(f"Error: {rec['sequence_id']} revised sequence_alignment does not match original")
        
        if len(rec['germline_alignment']) != len(rec['sequence_alignment']):
            reporter.report(f"Error: {rec['sequence_id']} revised germline_alignment and sequence alignments have different lengths")
        
        if rec['v_sequence_alignment'].replace('.', '') != v_sequence:
            reporter.report(f"Error: {rec['sequence_id']} revised v_sequence_alignment does not match original")

        writer.writerow(rec)

    reporter.report(f"Processed {rowcount} records")
    reporter.report(f"Errors: {reporter.error_count}, Warnings: {reporter.warning_count}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Input annotation file')
    parser.add_argument('germline_ref', help='Germline reference file')
    parser.add_argument('outfile', help='Output annotation file')
    parser.add_argument('logfile', help='Error log file')
    parser.add_argument('--condense_errors', action='store_true', help='Only report the first error for each similar sequence')
    parser.add_argument('--echo_errors', action='store_true', help='Report errors to stdout as well as in the log file')
    parser.add_argument('--mask_germline_np', action='store_true', help='Mask NP regionss in the germline alignment')

    args = parser.parse_args()
    main(args.infile, args.outfile, args.germline_ref, args.logfile, args.condense_errors, args.echo_errors, args.mask_germline_np)
