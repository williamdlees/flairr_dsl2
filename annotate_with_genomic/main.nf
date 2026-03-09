// FLAIRR-seq workflow

// these two params must be specified on the command line
params.reads = ""			// FASTA file containing the reads
params.sample_name = ""		// Sample name, to be used in reports and report filenames

// modify these params to meet your requirements
params.species = "Homo_sapiens"
params.locus = "IGH"
params.germline_ref_dir = "$baseDir/../../reference"
params.standard_ref_dir = "$baseDir/../../standard_reference" // The standard allele references (e.g., IMGT/OGRdb)
params.outdir = "$baseDir/../results"
//params.haplotype_genes = "IGHJ6,IGHD2-21,IGHD2-8"
params.haplotype_genes = "IGHJ6"

// these derived params should not need modifying
params.germline_ref = "${params.germline_ref_dir}/${params.species}_${params.locus}"
params.v_ref = "${params.germline_ref}V_gapped.fasta"
params.d_ref = "${params.germline_ref}D.fasta"
params.j_ref = "${params.germline_ref}J.fasta"
params.c_ref = "${params.germline_ref}C.fasta"
params.vdj_ref = "${params.germline_ref}VDJ.fasta"

params.standard_germline_ref = "${params.standard_ref_dir}/${params.species}_${params.locus}"
params.standard_v_ref = "${params.standard_germline_ref}V_gapped.fasta"
params.standard_d_ref = "${params.standard_germline_ref}D.fasta"
params.standard_j_ref = "${params.standard_germline_ref}J.fasta"
params.standard_c_ref = "${params.standard_germline_ref}C.fasta"
params.standard_vdj_ref = "${params.standard_germline_ref}VDJ.fasta"

params.allele_threshold_file = "${params.germline_ref}_allele_thresholds.tsv"
params.aux = "${params.germline_ref}.aux"
params.ndm = "${params.germline_ref}.ndm"
// this might be an issue

params.python_dir = "$baseDir/../python"

include { make_blast_db as make_blast_db_v; make_blast_db as make_blast_db_d; make_blast_db as make_blast_db_j; make_blast_db as make_blast_db_c } from '../modules/make_blast_db'
include { igblast_combo as igblast_combo1} from '../modules/igblast_combo'
include { align_v as align_v1} from '../modules/align_v'
include { define_clones } from '../modules/define_clones'
include { single_clone_representative } from '../modules/single_clone_representative'
include { haplotype_inference_report } from '../modules/haplotype_inference_report'
include { haplotype_const_report } from '../modules/haplotype_const_report'
include { ogrdbstats_report } from '../modules/ogrdbstats_report'

workflow {
	seqs = channel.fromPath(params.reads)
	
	igblast_combo1(seqs, params.v_ref, params.d_ref, params.j_ref, params.c_ref, params.aux, params.ndm)
	align_v1(igblast_combo1.out.output, params.v_ref, 'personalized-with-genomic-seq', params.python_dir)

	define_clones(align_v1.out.annotations, "$baseDir/../python/clonality_threshold.R", "$baseDir/../python/clone_stats.R")
	single_clone_representative(define_clones.out.output)

	haplotype_inference_report(align_v1.out.annotations, params.v_ref, params.d_ref, params.locus, params.haplotype_genes, single_clone_representative.out.ready)
	haplotype_const_report(align_v1.out.annotations, params.vdj_ref, params.python_dir, haplotype_inference_report.out.ready)
	ogrdbstats_report(align_v1.out.annotations, params.standard_v_ref, params.v_ref, params.locus, params.haplotype_genes, params.species, haplotype_const_report.out.ready)
}