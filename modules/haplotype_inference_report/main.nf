process haplotype_inference_report {
	// from rabhit: "Error in names(x) <- value : 'names' attribute [15] must be the same length as the vector [0]"
	// means that there aren't any single-assigned records in the set that can be genotyped with the siigned gene
	// to avoid this problem, don't do D haplotying with small test datasets

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_haplotype.tsv$/) "genotype_report/${name}_haplotype.tsv"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_binomDel.tsv$/) "genotype_report/${name}_binomDel.tsv"}
	
	input:
		path(airrFile)
		path(v_germline)
		path(d_germline)
		val(locus)
		val(haplotype_genes)
		val(ready)

	output:
		path("*_haplotype.tsv") optional true
		path("*_binomDel.tsv"), emit: deletions optional true
		val(true), emit: ready		

	script:
		name = params.sample_name
		v_germline = v_germline.name.startsWith('NO_FILE') ? "" : "${v_germline}"
		d_germline = d_germline.name.startsWith('NO_FILE') ? "" : "${d_germline}"
		outname = airrFile.name.toString().substring(0, airrFile.name.toString().indexOf("_db-pass"))
			
		template "haplotype_inference_report.R"
}