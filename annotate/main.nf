// FLAIRR-seq workflow

// these two params must be specified on the command line
params.reads = ""			// FASTA file containing the reads
params.sample_name = ""		// Sample name, to be used in reports and report filenames

// modify these params to meet your requirements
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

include { make_blast_db as make_blast_db_v; make_blast_db as make_blast_db_d; make_blast_db as make_blast_db_j; make_blast_db as make_blast_db_c } from '../modules/make_blast_db'
include { igblast } from '../modules/igblast'
include { makedb } from '../modules/makedb'
include { collapse_annotations } from '../modules/collapse_annotations'
include { create_germlines as create_germlines_pass1; create_germlines as create_germlines_pass2 } from '../modules/create_germlines'
include { define_clones } from '../modules/define_clones'
include { single_clone_representative } from '../modules/single_clone_representative'

workflow {
	make_blast_db_v(params.v_ref, true)
	make_blast_db_d(params.d_ref, make_blast_db_v.out.ready)
	make_blast_db_j(params.j_ref, make_blast_db_d.out.ready)
	make_blast_db_c(params.c_ref, make_blast_db_j.out.ready)
	
	seqs = channel.fromPath(params.reads)
	igblast(seqs, make_blast_db_v.out.blastdb, make_blast_db_d.out.blastdb, make_blast_db_j.out.blastdb, make_blast_db_c.out.blastdb, params.aux, params.ndm)
	makedb(seqs, igblast.out.output, params.v_ref, params.d_ref, params.j_ref, params.c_ref, 'non-personalized')
	collapse_annotations(makedb.out.annotations, "pass-1")
	create_germlines_pass1(collapse_annotations.out.output, params.v_ref, params.d_ref, params.j_ref, "false", "pass-1")
	define_clones(create_germlines_pass1.out.output)
	create_germlines_pass2(define_clones.out.output, params.v_ref, params.d_ref, params.j_ref, "true", "pass-2")
	single_clone_representative(create_germlines_pass2.out.output)
}