

process parse_headers {
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*${out}$/) "reports/$filename"}
	
	input:
		tuple val(name), path(reads)
		val prefix
		val method
		val act
		val args
		val ready

	output:
		tuple val(name), path("*${out}"), emit: output
		val(true), emit: ready

	script:
		outname = prefix + "_" + name
		if(method=="collapse" || method=="copy" || method=="rename" || method=="merge"){
			out="_reheader.fastq"
			act = (act=="none") ? "" : "--act ${act}"
			"""
			ParseHeaders.py  ${method} --outname ${outname} -s ${reads} ${args} ${act}
			"""
		}else{
			if(method=="table"){
					out=".tab"
					"""
					ParseHeaders.py ${method} --outname ${outname} -s ${reads} ${args}
					"""	
			}else{
				out="_reheader.fastq"
				"""
				ParseHeaders.py ${method} --outname ${outname} -s ${reads} ${args}
				"""		
			}
		}


}
