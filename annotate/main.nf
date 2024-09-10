// FLAIRR-seq workflow

// modify these params to meet your requirements
params.reads = "${baseDir}/../preprocess/results/reads/986-bc1003_1000_consensus-pass_reheader_collapse-unique_atleast-2.fasta"
params.species = "Homo_sapiens"
params.locus = "IGH"
params.germline_ref_dir = "$baseDir/../../reference"
params.outdir = "$baseDir/../results"

// these derived params should not need modifying
params.germline_ref = "${params.germline_ref_dir}/${params.species}_${params.locus}"
params.v_ref = "${params.germline_ref}V_gapped.fasta"
params.d_ref = "${params.germline_ref}D.fasta"
params.j_ref = "${params.germline_ref}J.fasta"
params.c_ref = "${params.germline_ref}C.fasta"
params.aux = "${params.germline_ref}.aux"
params.ndm = "${params.germline_ref}.ndm"
params.blast_db_dir = channel.fromPath("$baseDir/blast_db", type: 'dir').collect()

include { make_blast_db as make_blast_db_v; make_blast_db as make_blast_db_d; make_blast_db as make_blast_db_j; make_blast_db as make_blast_db_c } from '../modules/make_blast_db'
include { igblast } from '../modules/igblast'
include { makedb } from '../modules/makedb'
include { collapse_annotations } from '../modules/collapse_annotations'
include { create_germlines as create_germlines_pass1; create_germlines as create_germlines_pass2 } from '../modules/create_germlines'
include { define_clones } from '../modules/define_clones'
include { single_clone_representative } from '../modules/single_clone_representative'

workflow {
	make_blast_db_v(params.v_ref, params.blast_db_dir)
	make_blast_db_d(params.d_ref, params.blast_db_dir)
	make_blast_db_j(params.j_ref, params.blast_db_dir)
	make_blast_db_c(params.c_ref, params.blast_db_dir)
	
	seqs = channel.fromPath(params.reads)
	igblast(seqs, make_blast_db_v.out.blastdb, make_blast_db_d.out.blastdb, make_blast_db_j.out.blastdb, make_blast_db_c.out.blastdb, params.aux, params.ndm)
	makedb(seqs, igblast.out.output, params.v_ref, params.d_ref, params.j_ref, params.c_ref, 'non-personalized')
	collapse_annotations(makedb.out.annotations, "")
	create_germlines_pass1(collapse_annotations.out.output, params.v_ref, params.d_ref, params.j_ref)
	define_clones(create_germlines_pass1.out.output)
	create_germlines_pass2(define_clones.out.output, params.v_ref, params.d_ref, params.j_ref)
	single_clone_representative(create_germlines_pass2.out.output)
}