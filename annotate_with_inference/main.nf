// FLAIRR-seq annotation workflow, with genotyping and novel allele inference, using Tigger and Piglet

// modify these params to meet your requirements
params.reads = "${baseDir}/../preprocess/results/reads/986-bc1003_1000_consensus-pass_reheader_collapse-unique_atleast-2.fasta"
params.sample_name = "986-bc1003_1000"
params.species = "Homo_sapiens"
params.locus = "IGH"
params.germline_ref_dir = "$baseDir/../../reference"
params.outdir = "$baseDir/../results"

params.Undocumented_Alleles.germline_min = 1  // for test data


// these derived params should not need modifying but please check that the naming of the files matches what's expected here
// and that .ndm and .aux files are present
params.germline_ref = "${params.germline_ref_dir}/${params.species}_${params.locus}"
params.v_ref = "${params.germline_ref}V_gapped.fasta"
params.d_ref = "${params.germline_ref}D.fasta"
params.j_ref = "${params.germline_ref}J.fasta"
params.c_ref = "${params.germline_ref}C.fasta"
params.aux = "${params.germline_ref}.aux"
params.ndm = "${params.germline_ref}.ndm"

include { make_blast_db as make_blast_db_first_v; make_blast_db as make_blast_db_first_d; make_blast_db as make_blast_db_first_j ; make_blast_db as make_blast_db_first_c; make_blast_db as make_blast_db_second_v; make_blast_db as make_blast_db_third_v; make_blast_db as make_blast_db_third_d; make_blast_db as make_blast_db_third_j } from '../modules/make_blast_db'
include { igblast as igblast_first; igblast as igblast_second; igblast as igblast_third; } from '../modules/igblast'
include { makedb as makedb_first; makedb as makedb_second; makedb as makedb_third; } from '../modules/makedb'
include { collapse_annotations as collapse_annotations_first; collapse_annotations as collapse_annotations_second; collapse_annotations as collapse_annotations_third;  } from '../modules/collapse_annotations'
include { Undocumented_Alleles } from '../modules/Undocumented_Alleles'
include { create_germlines as create_germlines_first; create_germlines as create_germlines_second } from '../modules/create_germlines'
include { define_clones } from '../modules/define_clones'
include { single_clone_representative } from '../modules/single_clone_representative'
include { TIgGER_bayesian_genotype_Inference as tigger_j_call; TIgGER_bayesian_genotype_Inference as tigger_d_call; TIgGER_bayesian_genotype_Inference as tigger_v_call } from '../modules/tigger_bayesian_genotype_inference'
include { haplotype_inference_report } from '../modules/haplotype_inference_report'


workflow {
	seqs = channel.fromPath(params.reads)

	make_blast_db_first_v(params.v_ref, true)
	make_blast_db_first_d(params.d_ref, make_blast_db_first_v.out.ready)
	make_blast_db_first_j(params.j_ref, make_blast_db_first_d.out.ready)
	make_blast_db_first_c(params.c_ref, make_blast_db_first_j.out.ready)
	
	igblast_first(seqs, make_blast_db_first_v.out.blastdb, make_blast_db_first_d.out.blastdb, make_blast_db_first_j.out.blastdb, make_blast_db_first_c.out.blastdb, params.aux, params.ndm)
	makedb_first(seqs, igblast_first.out.output, params.v_ref, params.d_ref, params.j_ref, params.c_ref, 'non-personalized')

	collapse_annotations_first(makedb_first.out.annotations, "")
	Undocumented_Alleles(collapse_annotations_first.out.output, params.v_ref)
	make_blast_db_second_v(Undocumented_Alleles.out.novel_germline, true)	

	igblast_second(seqs, make_blast_db_second_v.out.blastdb, make_blast_db_first_d.out.blastdb, make_blast_db_first_j.out.blastdb, make_blast_db_first_c.out.blastdb, params.aux, params.ndm)
	makedb_second(seqs, igblast_second.out.output, Undocumented_Alleles.out.novel_germline, params.d_ref, params.j_ref, params.c_ref, 'v-personalised-1')

	collapse_annotations_second(makedb_second.out.annotations, "")	
	create_germlines_first(collapse_annotations_second.out.output, Undocumented_Alleles.out.novel_germline, params.d_ref, params.j_ref, "false")
	define_clones(create_germlines_first.out.output)
	create_germlines_second(define_clones.out.output, Undocumented_Alleles.out.novel_germline, params.d_ref, params.j_ref, "true")
	single_clone_representative(create_germlines_second.out.output)

	tigger_j_call('j_call', 'sequence_alignment', 'false', 'false', single_clone_representative.out.output, params.j_ref, "true")
	tigger_d_call('d_call', 'sequence_alignment', 'false', 'false', single_clone_representative.out.output, params.d_ref, tigger_j_call.out.ready)	
	tigger_v_call('v_call', 'sequence_alignment', 'false', 'false', single_clone_representative.out.output, params.v_ref, tigger_d_call.out.ready)	
	make_blast_db_third_v(tigger_v_call.out.personal_reference, tigger_v_call.out.ready)
	make_blast_db_third_d(tigger_d_call.out.personal_reference, make_blast_db_third_v.out.ready)
	make_blast_db_third_j(tigger_j_call.out.personal_reference, make_blast_db_third_d.out.ready)
	
	igblast_third(seqs, make_blast_db_third_v.out.blastdb, make_blast_db_third_d.out.blastdb, make_blast_db_third_j.out.blastdb, make_blast_db_first_c.out.blastdb, params.aux, params.ndm)
	makedb_third(seqs, igblast_third.out.output, tigger_v_call.out.personal_reference, tigger_d_call.out.personal_reference, tigger_j_call.out.personal_reference, params.c_ref, 'v-personalised-1')
	collapse_annotations_third(makedb_third.out.annotations, "final")	
	
	haplotype_inference_report(collapse_annotations_third.out.output, tigger_v_call.out.personal_reference, tigger_d_call.out.personal_reference, "false", "true")
	

/*
	create_germlines_pass2(define_clones.out.output, params.v_ref, params.d_ref, params.j_ref, "true")
	ogrdbstats_report(makedb_first, )
	vdjbase_genotype_report(Undocumented_Alleles.out.novel_germline, )

*/
}