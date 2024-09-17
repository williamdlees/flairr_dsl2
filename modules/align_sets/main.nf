// Process Parameters for align_sets:


process align_sets {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /AS_.*$/) "reports/$filename"}

	input:
		tuple val(name),file(reads)

	output:
		tuple val(name), file("*_align-pass.fastq"), emit: output
		tuple val(name), file("AS_*"), emit: log_file
		tuple val(name), file("*_align-fail.fastq") optional true
		tuple val(name), file("out*") optional true

	script:
		method = params.align_sets.method
		bf = params.align_sets.bf
		div = params.align_sets.div
		failed = params.align_sets.failed
		nproc = params.align_sets.nproc
		mate = params.mate

		muscle_exec = params.align_sets.muscle_exec

		offset_table = params.align_sets.offset_table
		pf = params.align_sets.pf
		mode = params.align_sets.mode

		primer_file = params.align_sets.primer_file
		reverse = params.align_sets.reverse

		readArray = reads.toString().split(' ')	

		reverse_arg = (reverse=="false") ? "" : "--reverse"
		div_arg = (div=="false") ? "" : "--div"
		failed_arg = (failed=="true") ? "--failed" : "" 
		bf = "--bf ${bf}"

		primer_file_argv = ""

		if(method=="offset"){
			pf = "--pf ${pf}"
			mode = "--mode ${mode}"
			offset_table_argv = "-d ${offset_table}"
			muscle_exec_argv = ""
		}else{
			pf = ""
			mode = ""
			offset_table_argv = ""
			muscle_exec_argv = "--exec ${muscle_exec}"
			
			if(method=="table"){
				primer_file_argv = "-p ${primer_file}"
			}
		}

		if(mate=="pair"){
			R1 = readArray[0]
			R2 = readArray[1]
			
			
			"""
			AlignSets.py ${method} -s ${R1} ${bf} ${muscle_exec_argv} ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R1_${name}.log --nproc ${nproc} >> out_${R1}_AS.log
			AlignSets.py ${method} -s ${R2} ${bf} ${muscle_exec_argv} ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_R2_${name}.log --nproc ${nproc} >> out_${R1}_AS.log
			"""
			
		}else{
			R1 = readArray[0]
			"""
			AlignSets.py ${method} -s ${R1} ${bf} ${muscle_exec_argv} ${div_arg} ${reverse_arg} ${failed_arg} ${pf} ${offset_table_argv} ${mode} ${primer_file_argv} --log AS_${name}.log --nproc ${nproc} >> out_${R1}_AS.log
			"""
		}

}
