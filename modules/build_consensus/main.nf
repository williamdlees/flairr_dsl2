

process build_consensus {
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /BC.*$/) "reports/$filename"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_consensus-fail.fastq$/) "failed_reads/$filename"}
	
	input:
		tuple val(name), path(reads)

	output:
		tuple val(name), path("*_consensus-pass.fastq"), emit: output
		tuple val(name), path("BC*"), emit: log_file
		tuple val(name), path("*_consensus-fail.fastq") optional true

	script:
		name = params.sample_name
		mate = params.mate
		failed = params.build_consensus.failed
		nproc = params.build_consensus.nproc
		barcode_field = params.build_consensus.barcode_field
		primer_field = params.build_consensus.primer_field
		act = params.build_consensus.act
		copy_field = params.build_consensus.copy_field
		mincount = params.build_consensus.mincount
		minqual = params.build_consensus.minqual
		minfreq = params.build_consensus.minfreq
		maxerror = params.build_consensus.maxerror
		prcons = params.build_consensus.prcons
		maxgap = params.build_consensus.maxgap
		maxdiv = params.build_consensus.maxdiv
		dep = params.build_consensus.dep
		//* @style @condition:{act="none",},{act="min",copy_field},{act="max",copy_field},{act="sum",copy_field},{act="set",copy_field},{act="majority",copy_field} @array:{barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep} @multicolumn:{failed,nproc},{barcode_field,primer_field,act,copy_field}, {mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep}

		args_values_bc = args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep)

		// args 
		if(isCollectionOrArray_bc(args_values_bc)){
			args_1 = args_values_bc[0]
			args_2 = args_values_bc[1]
		}else{
			args_1 = args_values_bc
			args_2 = args_values_bc
		}

		failed = (failed=="true") ? "--failed" : "" 


		if(mate=="pair"){
			// files
			readArray = reads.toString().split(' ')	
			R1 = readArray.grep(~/.*R1.*/)[0]
			R2 = readArray.grep(~/.*R2.*/)[0]
			
			"""
			BuildConsensus.py -s $R1 ${args_1} --outname ${name}_R1 --log BC_${name}_R1.log ${failed} --nproc ${nproc}
			BuildConsensus.py -s $R2 ${args_2} --outname ${name}_R2 --log BC_${name}_R2.log ${failed} --nproc ${nproc}
			"""
		}else{
			"""
			BuildConsensus.py -s $reads ${args_1} --outname ${name} --log BC_${name}.log ${failed} --nproc ${nproc}
			"""
		}
}


boolean isCollectionOrArray_bc(object) {    
    [Collection, Object[]].any { it.isAssignableFrom(object.getClass()) }
}

def args_creator_bc(barcode_field, primer_field, act, copy_field, mincount, minqual, minfreq, maxerror, prcons, maxgap, maxdiv, dep){
	def args_values;
    if(isCollectionOrArray_bc(barcode_field) || isCollectionOrArray_bc(primer_field) || isCollectionOrArray_bc(copy_field) || isCollectionOrArray_bc(mincount) || isCollectionOrArray_bc(minqual) || isCollectionOrArray_bc(minfreq) || isCollectionOrArray_bc(maxerror) || isCollectionOrArray_bc(prcons) || isCollectionOrArray_bc(maxgap) || isCollectionOrArray_bc(maxdiv) || isCollectionOrArray_bc(dep)){
    	primer_field = (isCollectionOrArray_bc(primer_field)) ? primer_field : [primer_field,primer_field]
    	act = (isCollectionOrArray_bc(act)) ? act : [act,act]
    	copy_field = (isCollectionOrArray_bc(copy_field)) ? copy_field : [copy_field,copy_field]
    	mincount = (isCollectionOrArray_bc(mincount)) ? mincount : [mincount,mincount]
    	minqual = (isCollectionOrArray_bc(minqual)) ? minqual : [minqual,minqual]
    	minfreq = (isCollectionOrArray_bc(minfreq)) ? minfreq : [minfreq,minfreq]
    	maxerror = (isCollectionOrArray_bc(maxerror)) ? maxerror : [maxerror,maxerror]
    	prcons = (isCollectionOrArray_bc(prcons)) ? prcons : [prcons,prcons]
    	maxgap = (isCollectionOrArray_bc(maxgap)) ? maxgap : [maxgap,maxgap]
    	maxdiv = (isCollectionOrArray_bc(maxdiv)) ? maxdiv : [maxdiv,maxdiv]
    	dep = (isCollectionOrArray_bc(dep)) ? dep : [dep,dep]
    	args_values = []
        [barcode_field,primer_field,act,copy_field,mincount,minqual,minfreq,maxerror,prcons,maxgap,maxdiv,dep].transpose().each { bf,pf,a,cf,mc,mq,mf,mr,pc,mg,md,d -> {
            bf = (bf=="") ? "" : "--bf ${bf}"
            pf = (pf=="") ? "" : "--pf ${pf}" 
            a = (a=="none") ? "" : "--act ${a}" 
            cf = (cf=="") ? "" : "--cf ${cf}" 
            mr = (mr=="none") ? "" : "--maxerror ${mr}" 
            pc = (pc=="none") ? "" : "--prcons ${pc}" 
            mg = (mg=="none") ? "" : "--maxgap ${mg}" 
            md = (md=="none") ? "" : "--maxdiv ${md}" 
            d = (d=="true") ? "--dep" : "" 
            args_values.add("${bf} ${pf} ${a} ${cf} -n ${mc} -q ${mq} --freq ${mf} ${mr} ${pc} ${mg} ${md} ${d}")
        }}
    }else{
        barcode_field = (barcode_field=="") ? "" : "--bf ${barcode_field}"
        primer_field = (primer_field=="") ? "" : "--pf ${primer_field}" 
        act = (act=="none") ? "" : "--act ${act}" 
        copy_field = (copy_field=="") ? "" : "--cf ${copy_field}" 
        maxerror = (maxerror=="none") ? "" : "--maxerror ${maxerror}" 
        prcons = (prcons=="none") ? "" : "--prcons ${prcons}" 
        maxgap = (maxgap=="none") ? "" : "--maxgap ${maxgap}" 
        maxdiv = (maxdiv=="none") ? "" : "--maxdiv ${maxdiv}" 
        dep = (dep=="true") ? "--dep" : "" 
        args_values = "${barcode_field} ${primer_field} ${act} ${copy_field} -n ${mincount} -q ${minqual} --freq ${minfreq} ${maxerror} ${prcons} ${maxgap} ${maxdiv} ${dep}"
    }
    return args_values
}


