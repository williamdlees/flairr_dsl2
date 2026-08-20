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
			.flatMap { row ->
				if (row.size() < 2) {
					error "Invalid --input row: ${row}. Expected sample_name<TAB>read_path"
				}
				def sample = row[0].toString().trim()
				def sample_reads = file(params.outdir)
					.resolve(sample)
					.resolve(output_locus)
					.resolve("reads")
					.resolve("${sample}_atleast-2.fasta")

				(sample_reads.exists() && sample_reads.size() > 0) ? [tuple(sample, sample_reads)] : []
			}

	def ref_v_ch = seqs.map { name, reads_path -> tuple(name, file(params.v_ref, checkIfExists: true)) }
	def ref_d_ch = seqs.map { name, reads_path ->
		def d_ref = file(params.d_ref)
		if (!d_ref.exists()) {
			log.warn "D reference not found: ${params.d_ref}; continuing with a placeholder reference."
		}
		tuple(name, d_ref)
	}
	def ref_j_ch = seqs.map { name, reads_path -> tuple(name, file(params.j_ref, checkIfExists: true)) }
	def first_ready_ch = seqs.map { name, reads_path -> tuple(name, true) }

	def igblast1_input = seqs.join(ref_v_ch).join(ref_d_ch).join(ref_j_ch)
	igblast_combo1(igblast1_input, params.c_ref, params.aux, params.ndm, params.python_dir)
	def align1_input = igblast_combo1.out.output.join(ref_v_ch)
	align_v1(align1_input, 'non-personalized', params.python_dir)

	def d_file_pers = ref_d_ch

	if(params.locus == 'IGH' || params.locus == 'TRB' || params.locus == 'TRD')
	{
		def tigger_j_input = align_v1.out.annotations.join(first_ready_ch)
		tigger_j_call('j_call', 'sequence_alignment', 'false', 'false', tigger_j_input, params.j_ref)
		def tigger_d_input = align_v1.out.annotations.join(tigger_j_call.out.ready)
		tigger_d_call('d_call', 'sequence_alignment', 'false', 'false', tigger_d_input, params.d_ref)	
		def tigger_v_input = align_v1.out.annotations.join(tigger_d_call.out.ready)
		tigger_v_call('v_call', 'sequence_alignment', 'false', 'false', tigger_v_input, params.v_ref)
		d_file_pers = tigger_d_call.out.personal_reference	
	} else {
		def tigger_j_input = align_v1.out.annotations.join(first_ready_ch)
		tigger_j_call('j_call', 'sequence_alignment', 'false', 'false', tigger_j_input, params.j_ref)
		def tigger_v_input = align_v1.out.annotations.join(tigger_j_call.out.ready)
		tigger_v_call('v_call', 'sequence_alignment', 'false', 'false', tigger_v_input, params.v_ref)
	}

	def igblast2_input = seqs.join(tigger_v_call.out.personal_reference).join(d_file_pers).join(tigger_j_call.out.personal_reference)
	igblast_combo2(igblast2_input, params.c_ref, params.aux, params.ndm, params.python_dir)
	def align2_input = igblast_combo2.out.output.join(tigger_v_call.out.personal_reference)
	align_v2(align2_input, 'personalized', params.python_dir)
	define_clones(align_v2.out.annotations, "$baseDir/../python/clonality_threshold.R", "$baseDir/../python/clone_stats.R")

	def define_ready_ch = define_clones.out.output.map { name, clone_file -> tuple(name, true) }
	def ogrdb_input = define_clones.out.output.join(igblast_combo1.out.consolidated_ref).join(tigger_v_call.out.personal_reference).join(define_ready_ch)
	ogrdbstats_report(ogrdb_input, params.locus, "", params.species)	
}