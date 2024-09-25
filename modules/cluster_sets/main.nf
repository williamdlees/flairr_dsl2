

process cluster_sets {

	input:
		tuple val(name),path(reads)
		val ready

	output:
		tuple val(name),path("*_cluster-pass.fastq"), emit: output
		tuple val(name),path("*_cluster-fail.fastq") optional true

	script:
		name = params.sample_name
		method = params.cluster_sets.method
		failed = params.cluster_sets.failed
		nproc = params.cluster_sets.nproc
		cluster_field = params.cluster_sets.cluster_field
		ident = params.cluster_sets.ident
		length = params.cluster_sets.length
		prefix = params.cluster_sets.prefix
		cluster_tool = params.cluster_sets.cluster_tool
		cluster_exec = params.cluster_sets.cluster_exec
		set_field = params.cluster_sets.set_field
		start = params.cluster_sets.start
		end = params.cluster_sets.end
		barcode_field = params.cluster_sets.barcode_field
		mate = params.mate

		method = (method.size==2) ? method : [method[0],method[0]]
		failed = (failed.size==2) ? failed : [failed[0],failed[0]]
		cluster_field = (cluster_field.size==2) ? cluster_field : [cluster_field[0],cluster_field[0]]
		ident = (ident.size==2) ? ident : [ident[0],ident[0]]
		length = (length.size==2) ? length : [length[0],length[0]]
		prefix = (prefix.size==2) ? prefix : [prefix[0],prefix[0]]
		cluster_tool = (cluster_tool.size==2) ? cluster_tool : [cluster_tool[0],cluster_tool[0]]
		cluster_exec = (cluster_exec.size==2) ? cluster_exec : [cluster_exec[0],cluster_exec[0]]
		set_field = (set_field.size==2) ? set_field : [set_field[0],set_field[0]]
		start = (start.size==2) ? start : [start[0],start[0]]
		end = (end.size==2) ? end : [end[0],end[0]]
		barcode_field = (barcode_field.size==2) ? barcode_field : [barcode_field[0],barcode_field[0]]

		def args_values = [];
		[method, failed, cluster_field, ident, length, prefix, cluster_tool, cluster_exec, set_field, start, end, barcode_field].transpose().each { m, f, cf, i, l, p, ct, ce, sf, s, e, bf -> {
			f = (f=="true") ? "--failed" : ""
			p = (p=="") ? "" : "--prefix ${p}" 
			ce = (ce=="") ? "" : "--exec ${ce}" 
			sf = (m=="set") ? "-f ${sf}" : ""
			s = (m=="barcode") ? "" : "--start ${s}" 
			e = (m=="barcode") ? "" : (e=="") ? "" : "--end ${e}" 
			bf = (m=="barcode") ? "-f ${bf}" : ""
			args_values.add("${m} ${f} -k ${cf} --ident ${i} --length ${l} ${p} --cluster ${ct} ${ce} ${sf} ${s} ${e} ${bf}")
		}}


		if(mate=="pair"){
			// files
			readArray = reads.toString().split(' ')	
			R1 = readArray.grep(~/.*R1.*/)[0]
			R2 = readArray.grep(~/.*R2.*/)[0]
			
			args_1 = args_values[0]
			args_2 = args_values[1]
			
			"""
			ClusterSets.py ${args_1} -s $R1  --nproc ${nproc}
			ClusterSets.py ${args_2} -s $R2  --nproc ${nproc}
			"""
		}else{
			args_1 = args_values[0]
			"""
			ClusterSets.py ${args_1} -s $reads --nproc ${nproc}
			"""
		}


}
