// FLAIRR-seq workflow

params.reads = "../../easton_short_samples/986-bc1003_1000.fastq"
params.outdir = "$baseDir/results"
params.nproc = 20
params.mate = "single"
params.run_FastQC = true
// params.constant_region is used to specify the constant region primer file to match. 
// for example GM matches the file FLAIRRSeq/primers/CPRIMERS_GM.fasta

params.constant_region = "GM"

// MaskPrimers_CPRIMERS
C_R1_primers = "$projectDir/primers/CPRIMERS_${params.constant_region}.fasta"
C_R2_primers = "$projectDir"		// can't be empty

C_MaskPrimers = [method: ["align"]]
C_MaskPrimers.mode = ["cut"]
C_MaskPrimers.start = [0]
C_MaskPrimers.extract_length = [22]
C_MaskPrimers.barcode = ["false"]
C_MaskPrimers.barcode_field = [""]
C_MaskPrimers.primer_field = ["CPRIMER"]
C_MaskPrimers.maxerror = [0.3]
C_MaskPrimers.revpr = ["false"]
C_MaskPrimers.maxlen = [50]
C_MaskPrimers.skiprc = ["false"]
C_MaskPrimers.failed = "true"
C_MaskPrimers.nproc = params.nproc
C_MaskPrimers.umi_length = ["0"]

// MaskPrimers_VPRIMERS
V_R1_primers = "$projectDir/primers/VPRIMERS.fasta"
V_R2_primers = "$projectDir/None.1"		// can't be empty

V_MaskPrimers =  [method: ["align"]]
V_MaskPrimers.method = ["align"]
V_MaskPrimers.mode = ["cut"]
V_MaskPrimers.start = [0]
V_MaskPrimers.extract_length = [22]
V_MaskPrimers.barcode = ["false"]
V_MaskPrimers.barcode_field = [""]
V_MaskPrimers.primer_field = ["VVprIMER"]
V_MaskPrimers.maxerror = [0.3]
V_MaskPrimers.revpr = ["false"]
V_MaskPrimers.maxlen = [50]
V_MaskPrimers.skiprc = ["false"]
V_MaskPrimers.failed = "true"
V_MaskPrimers.nproc = params.nproc
V_MaskPrimers.umi_length = ["0"]

// MaskPrimers_EXTRACT
E_R1_primers = "$projectDir/None.1"		// can't be empty
E_R2_primers = "$projectDir/None.2"

E_MaskPrimers =  [method: ["extract"]]
E_MaskPrimers.method = ["extract"]
E_MaskPrimers.mode = ["cut"]
E_MaskPrimers.start = [0]
E_MaskPrimers.extract_length = [22]
E_MaskPrimers.barcode = ["false"]
E_MaskPrimers.barcode_field = [""]
E_MaskPrimers.primer_field = ["BARCODE"]
E_MaskPrimers.maxerror = [0.3]
E_MaskPrimers.revpr = ["false"]
E_MaskPrimers.maxlen = [50]
E_MaskPrimers.skiprc = ["false"]
E_MaskPrimers.failed = "true"
E_MaskPrimers.nproc = params.nproc
E_MaskPrimers.umi_length = ["0"]

include { filter_seq_quality } from '../modules/filter_seq_quality'
include { parse_log as parse_log_FSQ; parse_log as parse_log_FSL; parse_log as parse_log_MPC; parse_log as parse_log_MPV; parse_log as parse_log_MPE } from '../modules/parse_log'
include { filter_seq_length } from '../modules/filter_seq_length'
include { FastQC } from '../modules/fastqc'
include { MaskPrimers as MaskPrimers_CPRIMERS; MaskPrimers as MaskPrimers_VPRIMERS; MaskPrimers as MaskPrimers_EPRIMERS; } from '../modules/mask_primers'

workflow {
	read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
	
	filter_seq_quality(read_pairs_ch)
	parse_log_FSQ(filter_seq_quality.out.log_file, "FSQ", "ID QUALITY")

	filter_seq_length(filter_seq_quality.out.output)
	parse_log_FSL(filter_seq_length.out.log_file, "FSL", "ID LENGTH")
	
	FastQC(filter_seq_length.out.output)
	
	MaskPrimers_CPRIMERS(filter_seq_length.out.output, C_MaskPrimers, C_R1_primers, C_R2_primers)
	parse_log_MPC(MaskPrimers_CPRIMERS.out.log_file, "MPC", "D PRIMER BARCODE ERROR")

	MaskPrimers_VPRIMERS(MaskPrimers_CPRIMERS.out.output, V_MaskPrimers, V_R1_primers, V_R2_primers)
	parse_log_MPV(MaskPrimers_VPRIMERS.out.log_file, "MPV", "D PRIMER BARCODE ERROR")

	MaskPrimers_EPRIMERS(MaskPrimers_VPRIMERS.out.output, E_MaskPrimers, E_R1_primers, E_R2_primers)
	parse_log_MPE(MaskPrimers_EPRIMERS.out.log_file, "MPE", "D PRIMER BARCODE ERROR")
}

