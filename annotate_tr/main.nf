// FLAIRR-seq workflow

// single-sample mode params
params.reads = null			// FASTA file containing the reads
params.sample_name = null		// Sample name, to be used in reports and report filenames
params.input = null			// TSV with sample_name<TAB>ignored_read_path

// modify these params to meet your requirements
params.species = "Homo_sapiens"
params.locus = "IGH"
params.output_locus = null
params.germline_ref_dir = "$baseDir/../../reference"
params.outdir = "${launchDir}/results"
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

include { igblast_combo as igblast_combo1; igblast_combo as igblast_combo2 } from '../modules/igblast_combo'
include { align_v as align_v1; align_v as align_v2 } from '../modules/align_v'
include { TIgGER_bayesian_genotype_Inference as tigger_j_call; TIgGER_bayesian_genotype_Inference as tigger_d_call; TIgGER_bayesian_genotype_Inference as tigger_v_call } from '../modules/tigger_bayesian_genotype_inference'
include { define_clones } from '../modules/define_clones'
include { ogrdbstats_report } from '../modules/ogrdbstats_report'

workflow {
	def output_locus = params.output_locus ?: params.locus
	def single_mode = params.reads && params.sample_name && !params.input
	def table_mode = params.input && !params.reads && !params.sample_name

	if (!single_mode && !table_mode) {
		error "Specify either --reads and --sample_name, or --input (TSV: sample_name<TAB>ignored_read_path)."
	}

	def seqs = single_mode \
		? channel
			.fromPath(params.reads, checkIfExists: true)
			.map { reads -> tuple(params.sample_name.toString(), reads) }
		: channel
			.fromPath(params.input, checkIfExists: true)
			.splitCsv(header: false, sep: '\t')
			.filter { row -> row && row.size() > 0 && row[0].toString().trim() && !row[0].toString().trim().startsWith('#') }
			.map { row ->
				if (row.size() < 2) {
					error "Invalid --input row: ${row}. Expected sample_name<TAB>read_path"
				}
				def sample = row[0].toString().trim()
				def sample_reads = file(params.outdir)
					.resolve(sample)
					.resolve(output_locus)
					.resolve("reads")
					.resolve("${sample}_atleast-2.fasta")
				tuple(sample, file(sample_reads.toString(), checkIfExists: true))
			}

	def ref_v_ch = seqs.map { name, reads_path -> tuple(name, file(params.v_ref, checkIfExists: true)) }
	def ref_d_ch = seqs.map { name, reads_path -> tuple(name, file(params.d_ref, checkIfExists: true)) }
	def ref_j_ch = seqs.map { name, reads_path -> tuple(name, file(params.j_ref, checkIfExists: true)) }
	def first_ready_ch = seqs.map { name, reads_path -> tuple(name, true) }

	igblast_combo1(seqs, ref_v_ch, ref_d_ch, ref_j_ch, params.c_ref, params.aux, params.ndm, params.python_dir)
	align_v1(igblast_combo1.out.output, ref_v_ch, 'non-personalized', params.python_dir)

	tigger_j_call('j_call', 'sequence_alignment', 'false', 'false', align_v1.out.annotations, params.j_ref, first_ready_ch)
	tigger_d_call('d_call', 'sequence_alignment', 'false', 'false', align_v1.out.annotations, params.d_ref, tigger_j_call.out.ready)	
	tigger_v_call('v_call', 'sequence_alignment', 'false', 'false', align_v1.out.annotations, params.v_ref, tigger_d_call.out.ready)	

	igblast_combo2(seqs, tigger_v_call.out.personal_reference, tigger_d_call.out.personal_reference, tigger_j_call.out.personal_reference, params.c_ref, params.aux, params.ndm, params.python_dir)
	align_v2(igblast_combo2.out.output, tigger_v_call.out.personal_reference, 'personalized', params.python_dir)
	define_clones(align_v2.out.annotations, "$baseDir/../python/clonality_threshold.R", "$baseDir/../python/clone_stats.R")
	def define_ready_ch = define_clones.out.output.map { name, clone_file -> tuple(name, true) }

	ogrdbstats_report(define_clones.out.output, igblast_combo1.out.consolidated_ref, tigger_v_call.out.personal_reference, params.locus, "", params.species, define_ready_ch)	
}