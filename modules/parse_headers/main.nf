

process parse_headers {
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*${out}$/) "reports/$filename"}
	
	input:
		tuple val(name), file(reads)
		val method
		val act
		val args

	output:
		tuple val(name), file("*${out}"), emit: output

	script:
		if(method=="collapse" || method=="copy" || method=="rename" || method=="merge"){
			out="_reheader.fastq"
			act = (act=="none") ? "" : "--act ${act}"
			"""
			ParseHeaders.py  ${method} -s ${reads} ${args} ${act}
			"""
		}else{
			if(method=="table"){
					out=".tab"
					"""
					ParseHeaders.py ${method} -s ${reads} ${args}
					"""	
			}else{
				out="_reheader.fastq"
				"""
				ParseHeaders.py ${method} -s ${reads} ${args}
				"""		
			}
		}


}
