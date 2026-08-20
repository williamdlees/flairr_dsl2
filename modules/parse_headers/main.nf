

process parse_headers {
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> (filename.endsWith("_reheader.fastq") || filename.endsWith(".tab")) ? "${name}/${params.output_locus ?: params.locus}/reports/${filename}" : null}
	
	input:
		tuple val(name), path(reads), val(ready)
		val prefix
		val method
		val act
		val args

	output:
		tuple val(name), path("{*.tab,*_reheader.fastq}"), emit: output
		tuple val(name), val(true), emit: ready

	script:
		outname = name + "_" + prefix
		if(method=="collapse" || method=="copy" || method=="rename" || method=="merge"){
			out="_reheader.fastq"
			act = (act=="none") ? "" : "--act ${act}"
			"""
			if [ ! -s ${reads} ]; then
				touch ${outname}${out}
			else
				ParseHeaders.py  ${method} --outname ${outname} -s ${reads} ${args} ${act}
			fi
			"""
		}else{
			if(method=="table"){
					out=".tab"
					"""
					if [ ! -s ${reads} ]; then
						touch ${outname}${out}
					else
						ParseHeaders.py ${method} --outname ${outname} -s ${reads} ${args}
					fi
					"""	
			}else{
				out="_reheader.fastq"
				"""
				if [ ! -s ${reads} ]; then
					touch ${outname}${out}
				else
					ParseHeaders.py ${method} --outname ${outname} -s ${reads} ${args}
				fi
				"""		
			}
		}


}
