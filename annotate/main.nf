// FLAIRR-seq workflow

// these two params must be specified on the command line
params.reads = ""			// FASTA file containing the reads
params.sample_name = ""		// Sample name, to be used in reports and report filenames

// modify these params to meet your requirements
params.species = "Homo_sapiens"
params.locus = "IGH"
params.germline_ref_dir = "$baseDir/../../reference"
params.outdir = "$baseDir/../results"
//params.haplotype_genes = "IGHJ6,IGHD2-21,IGHD2-8"
params.haplotype_genes = "IGHJ6"

// these derived params should not need modifying
params.germline_ref = "${params.germline_ref_dir}/${params.species}_${params.locus}"
params.v_ref = "${params.germline_ref}V_gapped.fasta"
params.d_ref = "${params.germline_ref}D.fasta"
params.j_ref = "${params.germline_ref}J.fasta"
params.c_ref = "${params.germline_ref}C.fasta"
params.aux = "${params.germline_ref}.aux"
params.ndm = "${params.germline_ref}.ndm"

params.python_dir = "$baseDir/../python"

include { make_blast_db as make_blast_db_v; make_blast_db as make_blast_db_d; make_blast_db as make_blast_db_j; make_blast_db as make_blast_db_c } from '../modules/make_blast_db'
include { igblast_combo as igblast_combo1; igblast_combo as igblast_combo2 } from '../modules/igblast_combo'
include { align_v as align_v1; align_v as align_v2 } from '../modules/align_v'
include { TIgGER_bayesian_genotype_Inference as tigger_j_call; TIgGER_bayesian_genotype_Inference as tigger_d_call; TIgGER_bayesian_genotype_Inference as tigger_v_call } from '../modules/tigger_bayesian_genotype_inference'
include { define_clones } from '../modules/define_clones'
include { single_clone_representative } from '../modules/single_clone_representative'
include { haplotype_inference_report } from '../modules/haplotype_inference_report'
include { ogrdbstats_report } from '../modules/ogrdbstats_report'

workflow {
	seqs = channel.fromPath(params.reads)
	
	igblast_combo1(seqs, params.v_ref, params.d_ref, params.j_ref, params.c_ref, params.aux, params.ndm)
	align_v1(igblast_combo1.out.output, params.v_ref, 'non-personalized', params.python_dir)

	tigger_j_call('j_call', 'sequence_alignment', 'false', 'false', align_v1.out.annotations, params.j_ref, "true")
	tigger_d_call('d_call', 'sequence_alignment', 'false', 'false', align_v1.out.annotations, params.d_ref, tigger_j_call.out.ready)	
	tigger_v_call('v_call', 'sequence_alignment', 'false', 'false', align_v1.out.annotations, params.v_ref, tigger_d_call.out.ready)	

	igblast_combo2(seqs, tigger_v_call.out.personal_reference, tigger_d_call.out.personal_reference, tigger_j_call.out.personal_reference, params.c_ref, params.aux, params.ndm)
	align_v2(igblast_combo2.out.output, tigger_v_call.out.personal_reference, 'personalized', params.python_dir)

	define_clones(align_v2.out.annotations, "$baseDir/../python/clonality_threshold.R")
	single_clone_representative(define_clones.out.output)

	haplotype_inference_report(align_v2.out.annotations, tigger_v_call.out.personal_reference, tigger_d_call.out.personal_reference, params.locus, params.haplotype_genes, single_clone_representative.out.ready)
	ogrdbstats_report(align_v2.out.annotations, igblast_combo1.out.consolidated_ref, tigger_v_call.out.personal_reference, params.locus, "", params.species, haplotype_inference_report.out.ready)	
}