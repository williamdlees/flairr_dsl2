# find the number of reads in each UMI using the log produced by mask_primers_extract

from receptor_utils import simple_bio_seq as simple

recs = simple.read_csv('F:/clareo/flairr_dsl2/processed_samples/1001-kinnex/results/reports/MP_m84248_240915_090555_s3_MPE_table.tab', delimiter='\t')

umis = {}

for rec in recs:
    umi = rec['PRIMER']  # it's the second tab, header is wrong
    if umi not in umis:
        umis[umi] = 0
    umis[umi] += 1

# write the umi counts to a csv file, sorted by size

umis = [{'UMI': k, 'Count': v} for k, v in sorted(umis.items(), key=lambda item: item[1], reverse=True)]
simple.write_csv('F:/clareo/flairr_dsl2/processed_samples/1001-kinnex/results/reports/MP_m84248_240915_090555_s3_MPE_table_umi_counts.csv', umis)
