# Build project-level table summarising allele calls for each sample

from receptor_utils import simple_bio_seq as simple
from argparse import ArgumentParser
import glob


def process_tigger_genotype(project_name, sample_name, genotype_file):
    genotype = simple.read_csv(genotype_file, delimiter='\t')
    results = []

    total_seqs_all_genes = {}
    for row in genotype:
        gene = row["gene"]
        gene_type = gene[3]
        locus = gene[:3]

        kdiff = row["k_diff"]
        # allele counts come from un-genotyped column
        allele_counts = {}
        ac = str(row["counts"]).split(',')
        total_count = sum([int(x) if x.isnumeric() else 0 for x in ac])
        for index, allele in enumerate(str(row["alleles"]).split(",")):
            if not ((allele == "Unk") or ("NR" in allele) or ("del" in allele.lower())) and index < len(ac):
                allele_counts[allele] = ac[index]
        for index, allele in enumerate(str(row["genotyped_alleles"]).split(",")):
            if (allele == "Unk") or ("NR" in allele):
                continue
            elif len(str(allele)) == 1:
                allele = "0" + str(allele)
            elif ("del" in allele.lower()):
                allele = "Del"

            # check if the allele exists in the genotype
            count = 0
            if allele != "Del":
                count = int(allele_counts[allele]) if allele in allele_counts else 0

            results.append({
                "project": project_name,
                "subject": sample_name,
                "sample_name": sample_name,
                "num_alleles": len(row["genotyped_alleles"].split(",")),
                "locus": locus,
                "gene": gene,
                "gene_type": gene_type,
                "vdjbase_allele": gene + '*' + allele,
                "seqs": count,
                "gene_seqs": total_count,
                "frequency": 0,
                "gene_frequency": 0,
                "kdiff": kdiff
            })

            if gene_type not in total_seqs_all_genes:
                total_seqs_all_genes[gene_type] = 0
            total_seqs_all_genes[gene_type] += count

    for result in results:
        result["frequency"] = f"{result['seqs'] / total_seqs_all_genes[result['gene_type']]:.5f}"
        result["gene_frequency"] = f"{result['gene_seqs'] / total_seqs_all_genes[result['gene_type']]:.5f}"

    return results


args = ArgumentParser()
args.add_argument("project", help="Project name")
args.add_argument("indir", help="Input directory")
args.add_argument("outfile", help="Output file")
args = args.parse_args()

genotype_files = glob.glob("**/*genotype_report*.tsv", recursive=True, root_dir=args.indir)

results = []

for genotype_file in genotype_files:
    sample_name = genotype_file.replace("\\", "/").split("/")[0]
    results.extend(process_tigger_genotype(args.project, sample_name, args.indir + '/' + genotype_file))

simple.write_csv(args.outfile, results)
