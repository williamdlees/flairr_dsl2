$HOSTNAME = ""
params.outdir = 'results'  

// filter_seq_quality
params.filter_seq_quality.method = "quality"
params.filter_seq_quality.nproc = params.nproc
params.filter_seq_quality.q = "20"
params.filter_seq_quality.n_length = ""
params.filter_seq_quality.n_missing = ""

// filter_seq_length
params.filter_seq_length.method = "length"
params.filter_seq_length.nproc = params.nproc
params.filter_seq_length.q = ""
params.filter_seq_length.n_length = "750"
params.filter_seq_length.n_missing = ""

// MaskPrimers_CPRIMERS
params.MaskPrimers_CPRIMERS.method = ["align"]
params.MaskPrimers_CPRIMERS.mode = ["cut"]
params.MaskPrimers_CPRIMERS.start = [0]
params.MaskPrimers_CPRIMERS.barcode = ["false"]
params.MaskPrimers_CPRIMERS.barcode_field = [""]
params.MaskPrimers_CPRIMERS.primer_field = ["CPRIMER"]
params.MaskPrimers_CPRIMERS.maxerror = [0.3]
params.MaskPrimers_CPRIMERS.revpr = ["false"]
params.MaskPrimers_CPRIMERS.maxlen = [50]
params.MaskPrimers_CPRIMERS.skiprc = ["false"]
params.MaskPrimers_CPRIMERS.failed = "true"
params.MaskPrimers_CPRIMERS.nproc = params.nproc
params.MaskPrimers_CPRIMERS.R1_primers = "${params.projectDir}/primers/CPRIMERS_${params.constant_region}.fasta"
params.MaskPrimers_CPRIMERS.R2_primers = ""

params.parse_log_MP_CPRIMERS.suffix = "_CPRIMERS"

// MaskPrimers_VPRIMERS
params.MaskPrimers_VPRIMERS.method = ["align"]
params.MaskPrimers_VPRIMERS.mode = ["cut"]
params.MaskPrimers_VPRIMERS.start = [0]
params.MaskPrimers_VPRIMERS.barcode = ["false"]
params.MaskPrimers_VPRIMERS.barcode_field = [""]
params.MaskPrimers_VPRIMERS.primer_field = ["VPRIMER"]
params.MaskPrimers_VPRIMERS.maxerror = [0.3]
params.MaskPrimers_VPRIMERS.revpr = ["false"]
params.MaskPrimers_VPRIMERS.maxlen = [50]
params.MaskPrimers_VPRIMERS.skiprc = ["false"]
params.MaskPrimers_VPRIMERS.failed = "true"
params.MaskPrimers_VPRIMERS.nproc = params.nproc
params.MaskPrimers_VPRIMERS.R1_primers = "${params.projectDir}/primers/VPRIMERS.fasta"
params.MaskPrimers_VPRIMERS.R2_primers = ""

params.parse_log_MP_VPRIMERS.suffix = "_VPRIMERS"
// MaskPrimers_EXTRACT
params.MaskPrimers_EXTRACT.method = ["extract"]
params.MaskPrimers_EXTRACT.mode = ["cut"]
params.MaskPrimers_EXTRACT.start = [0]
params.MaskPrimers_EXTRACT.extract_length = [22]
params.MaskPrimers_EXTRACT.barcode = ["false"]
params.MaskPrimers_EXTRACT.barcode_field = [""]
params.MaskPrimers_EXTRACT.primer_field = ["BARCODE"]
params.MaskPrimers_EXTRACT.maxerror = [0.3]
params.MaskPrimers_EXTRACT.revpr = ["false"]
params.MaskPrimers_EXTRACT.maxlen = [50]
params.MaskPrimers_EXTRACT.skiprc = ["false"]
params.MaskPrimers_EXTRACT.failed = "true"
params.MaskPrimers_EXTRACT.nproc = params.nproc
params.MaskPrimers_EXTRACT.R1_primers = ""
params.MaskPrimers_EXTRACT.R2_primers = ""

params.parse_log_MP_EXTRACT.suffix = "_EXTRACT"

params.check_for_seqs.primers_file = "${params.projectDir}/primers/VPRIMERS.fasta"

// align_sets
params.align_sets.method = "muscle"
params.align_sets.bf = "BARCODE"
params.align_sets.div = "fasle"
params.align_sets.failed = "false"
params.align_sets.nproc = params.nproc
params.align_sets.muscle_exec = "/usr/local/bin/muscle"
params.align_sets.offset_table = ""
params.align_sets.pf = ""
params.align_sets.mode = ""
params.align_sets.primer_file = ""
params.align_sets.reverse = "false"

// cluster_sets
params.cluster_sets.method = ["set"]
params.cluster_sets.failed = ["false"]
params.cluster_sets.nproc = params.nproc
params.cluster_sets.cluster_field = ["CLUSTER"]
params.cluster_sets.cluster_tool = ["usearch"]
params.cluster_sets.cluster_exec = ["/usr/local/bin/usearch"]
params.cluster_sets.set_field = ["BARCODE"]
params.cluster_sets.barcode_field = ["BARCODE"]

// parse_headers_copy
params.parse_headers_copy.method = "copy"
params.parse_headers_copy.act = "cat"
params.parse_headers_copy.args = "-f BARCODE -k CLUSTER"

// build_consensus
params.build_consensus.failed = "true"
params.build_consensus.nproc = params.nproc
params.build_consensus.barcode_field = ["CLUSTER"]
params.build_consensus.primer_field = ["CPRIMER"]
params.build_consensus.act = ["none"]
params.build_consensus.copy_field = [""]
params.build_consensus.maxerror = ["0.1"]
params.build_consensus.prcons = ["0.6"]
params.build_consensus.maxgap = ["0.5"]
params.build_consensus.maxdiv = ["none"]
params.build_consensus.dep = ["false"]

// parse_headers_copy
params.parse_headers_collapse.method = "collapse"
params.parse_headers_collapse.act = "min"
params.parse_headers_collapse.args = "-f CONSCOUNT"

// collapse_seq
params.collapse_seq.max_missing = 20 
params.collapse_seq.inner = "true" 
params.collapse_seq.fasta = "true" 
params.collapse_seq.act = "sum" 
params.collapse_seq.uf = "CREGION" 
params.collapse_seq.cf = "CONSCOUNT" 
params.collapse_seq.nproc = params.nproc
params.collapse_seq.failed = "true"

// split_seq
params.split_seq.field = "CONSCOUNT"
params.split_seq.num = 2
params.split_seq.fasta = "true"

// parse_headers_table
params.parse_headers_table.method = "table"
params.parse_headers_table.args = "-f ID PRCONS CONSCOUNT DUPCOUNT"

params.Alignment_FLAIRSEQ_IgBlastn.num_threads = params.nproc
params.Alignment_FLAIRSEQ_IgBlastn.ig_seqtype = "Ig"
params.Alignment_FLAIRSEQ_IgBlastn.outfmt = "MakeDb"
params.Alignment_FLAIRSEQ_IgBlastn.num_alignments_V = "10"
params.Alignment_FLAIRSEQ_IgBlastn.domain_system = "imgt"


params.Alignment_FLAIRSEQ_MakeDb.failed = "true"
params.Alignment_FLAIRSEQ_MakeDb.format = "airr"
params.Alignment_FLAIRSEQ_MakeDb.regions = "default"
params.Alignment_FLAIRSEQ_MakeDb.extended = "true"
params.Alignment_FLAIRSEQ_MakeDb.asisid = "false"
params.Alignment_FLAIRSEQ_MakeDb.asiscalls = "false"
params.Alignment_FLAIRSEQ_MakeDb.inferjunction = "false"
params.Alignment_FLAIRSEQ_MakeDb.partial = "false"
params.Alignment_FLAIRSEQ_MakeDb.name_alignment = "_Alignment_FLAIRSEQ"

// Process Parameters for Alignment_FLAIRSEQ_Collapse_AIRRseq:
params.Alignment_FLAIRSEQ_Collapse_AIRRseq.conscount_min = 2
params.Alignment_FLAIRSEQ_Collapse_AIRRseq.n_max = 10
params.Alignment_FLAIRSEQ_Collapse_AIRRseq.name_alignment = "_Alignment_FLAIRSEQ"

// Process Parameters for Clone_AIRRseq_First_CreateGermlines:
params.Clone_AIRRseq_First_CreateGermlines.failed = "false"
params.Clone_AIRRseq_First_CreateGermlines.format = "airr"
params.Clone_AIRRseq_First_CreateGermlines.g = "dmask"
params.Clone_AIRRseq_First_CreateGermlines.cloned = "false"
params.Clone_AIRRseq_First_CreateGermlines.seq_field = ""
params.Clone_AIRRseq_First_CreateGermlines.v_field = ""
params.Clone_AIRRseq_First_CreateGermlines.d_field = ""
params.Clone_AIRRseq_First_CreateGermlines.j_field = ""
params.Clone_AIRRseq_First_CreateGermlines.clone_field = ""

params.Clone_AIRRseq_DefineClones.failed = "false"
params.Clone_AIRRseq_DefineClones.format = "airr"
params.Clone_AIRRseq_DefineClones.seq_field = ""
params.Clone_AIRRseq_DefineClones.v_field = ""
params.Clone_AIRRseq_DefineClones.d_field = ""
params.Clone_AIRRseq_DefineClones.j_field = ""
params.Clone_AIRRseq_DefineClones.group_fields =  ""
params.Clone_AIRRseq_DefineClones.mode = "gene"
params.Clone_AIRRseq_DefineClones.dist = "0.2"
params.Clone_AIRRseq_DefineClones.norm = "len"
params.Clone_AIRRseq_DefineClones.act = "set"
params.Clone_AIRRseq_DefineClones.model = "hh_s5f"
params.Clone_AIRRseq_DefineClones.sym = "min"
params.Clone_AIRRseq_DefineClones.link = "single"
params.Clone_AIRRseq_DefineClones.maxmiss = "0"

// Process Parameters for Clone_AIRRseq_Second_CreateGermlines:
params.Clone_AIRRseq_Second_CreateGermlines.failed = "false"
params.Clone_AIRRseq_Second_CreateGermlines.format = "airr"
params.Clone_AIRRseq_Second_CreateGermlines.g = "dmask"
params.Clone_AIRRseq_Second_CreateGermlines.cloned = "true"
params.Clone_AIRRseq_Second_CreateGermlines.seq_field = ""
params.Clone_AIRRseq_Second_CreateGermlines.v_field = ""
params.Clone_AIRRseq_Second_CreateGermlines.d_field = ""
params.Clone_AIRRseq_Second_CreateGermlines.j_field = ""
params.Clone_AIRRseq_Second_CreateGermlines.clone_field = ""



if (!params.mate){params.mate = ""} 
if (!params.reads){params.reads = ""} 
if (!params.d_germline){params.d_germline = ""} 
if (!params.j_germline){params.j_germline = ""} 
if (!params.v_germline){params.v_germline = ""} 
if (!params.alignment_mate){params.alignment_mate = ""} 
if (!params.auxiliary_data){params.auxiliary_data = ""} 
if (!params.custom_internal_data){params.custom_internal_data = ""} 
if (!params.c_germline){params.c_germline = ""} 
// Stage empty file to be used as an optional input where required
ch_empty_file_1 = file("$baseDir/.emptyfiles/NO_FILE_1", hidden:true)
ch_empty_file_2 = file("$baseDir/.emptyfiles/NO_FILE_2", hidden:true)
ch_empty_file_3 = file("$baseDir/.emptyfiles/NO_FILE_3", hidden:true)
