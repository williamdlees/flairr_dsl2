process haplotype_inference_report {
	// from rabhit: "Error in names(x) <- value : 'names' attribute [15] must be the same length as the vector [0]"
	// means that there aren't any single-assigned records in the set that can be genotyped with the siigned gene
	// to avoid this problem, don't do D haplotying with small test datasets

	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_haplotype.tsv$/) "${name}/${params.output_locus ?: params.locus}/genotype_report/${name}_haplotype.tsv"}
	publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*_binomDel.tsv$/) "${name}/${params.output_locus ?: params.locus}/genotype_report/${name}_binomDel.tsv"}
	
	input:
		tuple val(name), path(airrFile)
		tuple val(name_v), path(v_germline)
		tuple val(name_d), path(d_germline)
		val(locus)
		val(haplotype_genes)
		tuple val(name_ready), val(ready)

	output:
		tuple val(name), path("*_haplotype.tsv"), optional: true
		tuple val(name), path("*_binomDel.tsv"), emit: deletions, optional: true
		tuple val(name), val(true), emit: ready		

	script:
		v_germline = v_germline.name.startsWith('NO_FILE') ? "" : "${v_germline}"
		d_germline = d_germline.name.startsWith('NO_FILE') ? "" : "${d_germline}"
		outname = airrFile.getBaseName() + "_haplotype.tsv"
			
		template "haplotype_inference_report.R"
}