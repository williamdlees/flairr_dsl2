

process MaskPrimers {

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_primers-fail.fastq$/) "failed_reads/$filename"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /MP_.*$/) "reports/$filename"}

	input:
		tuple val(name), path(reads)
		val mask_params
		path R1_primers
		path R2_primers
		val ready

	output:
		tuple val(name), path("*_primers-pass.fastq"), emit: output
		tuple val(name), path("*_primers-fail.fastq") optional true
		tuple val(name), path("MP_*"), emit: log_file						// to parse_log
		tuple val(name), path("out*")

	script:
		mate = mask_params.mate
		method = mask_params.method
		barcode_field = mask_params.barcode_field
		primer_field = mask_params.primer_field
		barcode = mask_params.barcode
		revpr = mask_params.revpr
		mode = mask_params.mode
		failed = mask_params.failed
		nproc = mask_params.nproc
		maxerror = mask_params.maxerror
		umi_length = mask_params.umi_length
		start = mask_params.start
		extract_length = mask_params.extract_length
		maxlen = mask_params.maxlen
		skiprc = mask_params.skiprc
		//R1_primers = mask_params.R1_primers
		//R2_primers = mask_params.R2_primers

		method = (method.collect().size==2) ? method : [method[0],method[0]]
		barcode_field = (barcode_field.collect().size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]
		primer_field = (primer_field.collect().size==2) ? primer_field : [primer_field[0],primer_field[0]]
		barcode = (barcode.collect().size==2) ? barcode : [barcode[0],barcode[0]]
		revpr = (revpr.collect().size==2) ? revpr : [revpr[0],revpr[0]]
		mode = (mode.collect().size==2) ? mode : [mode[0],mode[0]]
		maxerror = (maxerror.collect().size==2) ? maxerror : [maxerror[0],maxerror[0]]
		umi_length = (umi_length.collect().size==2) ? umi_length : [umi_length[0],umi_length[0]]
		start = (start.collect().size==2) ? start : [start[0],start[0]]
		extract_length = (extract_length.collect().size==2) ? extract_length : [extract_length[0],extract_length[0]]
		maxlen = (maxlen.collect().size==2) ? maxlen : [maxlen[0],maxlen[0]]
		skiprc = (skiprc.collect().size==2) ? skiprc : [skiprc[0],skiprc[0]]
		failed = (failed=="true") ? "--failed" : ""

		def args_values = [];
		[method,barcode_field,primer_field,barcode,revpr,mode,maxerror,umi_length,start,extract_length,maxlen,skiprc].transpose().each { m,bf,pf,bc,rp,md,mr,ul,s,el,ml,sk -> {
			
			if(m=="align"){
				s = ""
			}else{
				if(bc=="false"){
					s = "--start ${s}"
				}else{
					s = s + ul
					s = "--start ${s}"
				}
			}
			
			el = (m=="extract") ? "--len ${el}" : ""
			mr = (m=="extract") ? "" : "--maxerror ${mr}" 
			ml = (m=="align") ? "--maxlen ${ml}" : "" 
			sk = (m=="align" && sk=="true") ? "--skiprc" : "" 
			
			PRIMER_FIELD = "${pf}"
			
			// all
			bf = (bf=="") ? "" : "--bf ${bf}"
			pf = (pf=="") ? "" : "--pf ${pf}"
			bc = (bc=="false") ? "" : "--barcode"
			rp = (rp=="false") ? "" : "--revpr"
			args_values.add("${m} --mode ${md} ${bf} ${pf} ${bc} ${rp} ${mr} ${s} ${el} ${ml} ${sk}")
			
			
		}}

		readArray = reads.toString().split(' ')
		if(mate=="pair"){
			args_1 = args_values[0]
			args_2 = args_values[1]
			
		  


			R1 = readArray[0]
			R2 = readArray[1]
			
			R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
			R2_primers = (method[1]=="extract") ? "" : "-p ${R2_primers}"
			
			
			"""
			
			MaskPrimers.py ${args_1} -s ${R1} ${R1_primers} --log MP_R1_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
			MaskPrimers.py ${args_2} -s ${R2} ${R2_primers} --log MP_R2_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
			"""
		}else{
			args_1 = args_values[0]
			
			R1_primers = (method[0]=="extract") ? "" : "-p ${R1_primers}"
			
			R1 = readArray[0]

			"""
			echo -e "Assuming inputs for R1\n"
			
			MaskPrimers.py ${args_1} -s ${reads} ${R1_primers} --log MP_${name}.log  --nproc ${nproc} ${failed} >> out_${R1}_MP.log
			"""
		}

}
