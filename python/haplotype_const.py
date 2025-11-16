import argparse
import csv
from receptor_utils import simple_bio_seq as simple
import os.path
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Allele Threshold Table taken from PIgLET - see https://bitbucket.org/yaarilab/piglet/src/master/ for atrribution and terms


def run_haplotype(haplotype_gene, annot_file):
    haplotype_data = {}

    with open(annot_file, 'r') as f_annot:
        reader = csv.DictReader(f_annot, delimiter='\t')
        for rec in reader:
            if ',' in rec['c_call'] or not rec['c_call'] or not rec['c_sequence_start'] or not rec['c_sequence_end']:
                continue

            c_gene = rec['c_call'].split('*')[0]

            if c_gene != haplotype_gene:
                continue

            if rec['c_call'] not in haplotype_data:
                haplotype_data[rec['c_call']] = {
                    'total': 1,
                    'v_calls': {},
                    'j_calls': {},
                    'd_calls': {}
                }

            found_calls = False
            for gtype in ['v', 'd', 'j']:        
                if rec[f'{gtype}_call'] and ',' not in rec[f'{gtype}_call']: 
                    found_calls = True
                    gene = rec[f'{gtype}_call'].split('*')[0]
                    if gene not in haplotype_data[rec['c_call']][f'{gtype}_calls']:
                        haplotype_data[rec['c_call']][f'{gtype}_calls'][gene] = {}
                    if rec[f'{gtype}_call'] not in haplotype_data[rec['c_call']][f'{gtype}_calls'][gene]:
                        haplotype_data[rec['c_call']][f'{gtype}_calls'][gene][rec[f'{gtype}_call']] = 0
                    haplotype_data[rec['c_call']][f'{gtype}_calls'][gene][rec[f'{gtype}_call']] += 1

            if found_calls:
                haplotype_data[rec['c_call']]['total'] += 1

    return haplotype_data


ANCHOR_FILTER_THRESHOLD = 0.6  # 60% of max


def filter_haplotype_data(haplotype_data):
    # FInd the c_call with the highest count. Then remove any c_calls with less than ANCHOR_FILTER_THRESHOLD of that count.
    if not haplotype_data:
        return haplotype_data
    
    max_count = max(data['total'] for data in haplotype_data.values())
    threshold = max_count * ANCHOR_FILTER_THRESHOLD
    filtered_data = {c_call: data for c_call, data in haplotype_data.items() if data['total'] >= threshold}
    return filtered_data


DEFAULT_THRESHOLD = 1e-04       # default threshold if allele is not in the table


def z_score(allele_name, allele_count, repertoire_depth, allele_thresholds):
    # Calculate z-score for an allele count given the repertoire depth and allele thresholds
    threshold = allele_thresholds.get(allele_name, DEFAULT_THRESHOLD)
    return (allele_count - (repertoire_depth * threshold)) / (repertoire_depth * threshold * (1 - threshold))**0.5


def write_haplotype_report(haplotype_data, anchors, gene_list, allele_thresholds, report_file):
    # for each gene type v, d, j in turn, write a csv file with one line per v,d,j gene
    # the columns should be: gene, alleles, total_count, count_per_allele where alleles and count_per_allele are comma-separated lists
    # geneses within each v,d,j type should be sorted alphabetically

    rows = []
    hap_counts = [0] * len(anchors)

    for gtype in ['v', 'd', 'j']:
        for gene in sorted(gene_list[gtype]):
            total_count = 0
            for i in range(len(anchors)):
                if gene in haplotype_data[anchors[i]][f'{gtype}_calls']:
                    allele_counts = haplotype_data[anchors[i]][f'{gtype}_calls'][gene]
                    gene_count = sum(allele_counts.values())
                    hap_counts[i] += gene_count
                    total_count += gene_count

            row = {
                'gene': gene,
                'total': total_count
            }

            for i in range(len(anchors)):
                if gene in haplotype_data[anchors[i]][f'{gtype}_calls']:
                    allele_counts = haplotype_data[anchors[i]][f'{gtype}_calls'][gene]
                    alleles = ';'.join(sorted(allele_counts.keys()))
                    counts = ';'.join(str(allele_counts[allele]) for allele in sorted(allele_counts.keys()))
                    gene_count = sum(allele_counts.values())
                    row[anchors[i]] = alleles
                    row[f'{anchors[i]} count'] = counts
                    row[f'{anchors[i]} %'] = f'{(gene_count / total_count * 100):.2f}' if total_count > 0 else '0.00'
                else:
                    row[anchors[i]] = ''
                    row[f'{anchors[i]} count'] = ''
                    row[f'{anchors[i]} %'] = ''

            rows.append(row)

    # Calculate z scores and haplotype assignments
    for row in rows:
        for i in range(len(anchors)):
            z_scores = []
            hap_alleles = []
            if row[f'{anchors[i]} count']:
                for allele_count, allele_name in zip(row[f'{anchors[i]} count'].split(';'), row[anchors[i]].split(';')):
                    z = z_score(allele_name, int(allele_count), hap_counts[i], allele_thresholds)
                    z_scores.append(f'{z:.2f}')
                    if z > 0:
                        hap_alleles.append(allele_name)
            row[f'{anchors[i]} z'] = ';'.join(z_scores)
            row[f'{anchors[i]} hap'] = ';'.join(hap_alleles) if hap_alleles else ''

    with open(report_file, 'w', newline='') as f_report:
        fieldnames = ['gene', 'total', anchors[0], anchors[1], f'{anchors[0]} count', f'{anchors[1]} count', f'{anchors[0]} %', f'{anchors[1]} %', f'{anchors[0]} z', f'{anchors[1]} z', f'{anchors[0]} hap', f'{anchors[1]} hap']
        writer = csv.DictWriter(f_report, fieldnames=fieldnames, delimiter=',')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# scan the annotation file for possible haplotype genes. Return a set of gene names and their counts.
# ambiguous calls are ignored
def scan_for_genes(annot_file):
    gene_counts = {}
    with open(annot_file, 'r') as f_annot:
        reader = csv.DictReader(f_annot, delimiter='\t')
        for rec in reader:
            if rec['c_call'] and ',' not in rec['c_call']:
                gene = rec['c_call'].split('*')[0]
                if gene not in gene_counts:
                    gene_counts[gene] = 0
                gene_counts[gene] += 1
    return gene_counts


def make_gene_list(vdj_ref_file):
    gene_list = {}
    refs = simple.read_fasta(vdj_ref_file)
    for allele_name in refs.keys():
        gene_name = allele_name.split('*')[0]
        gene_type = gene_name[3].lower()
        if gene_type not in gene_list:
            gene_list[gene_type] = set()
        if gene_name not in gene_list[gene_type]:
            gene_list[gene_type].add(gene_name)

    return gene_list


def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.isfile(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x


def plot_haplotype_data(csv_file, output_file=None):
    """
    Create a Plotly visualization of haplotype data with two columns (one per anchor gene).
    Each column contains horizontal bars for V, D, J genes, colored by haplotyped alleles.
    
    Args:
        csv_file: Path to the haplotype report CSV file
        output_file: Optional path to save the HTML plot (if None, uses csv_file name with .html)
    """
    # Read the CSV file
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    if not rows:
        print(f"No data found in {csv_file}")
        return
    
    # Extract anchor names from column headers (columns ending with ' hap')
    anchors = []
    for col in rows[0].keys():
        if col.endswith(' hap'):
            anchors.append(col.replace(' hap', ''))
    
    if len(anchors) != 2:
        print(f"Expected 2 anchor genes, found {len(anchors)}: {anchors}")
        return
    
    # Create subplots with 2 columns
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=[f"<b>{anchors[0]}</b>", f"<b>{anchors[1]}</b>"],
        horizontal_spacing=0.15
    )
    
    # Color palette for alleles
    color_palette = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
        '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5'
    ]
    
    # Collect all unique alleles across both anchors to assign consistent colors
    all_alleles = set()
    for row in rows:
        for anchor in anchors:
            hap_col = f'{anchor} hap'
            if row[hap_col]:
                alleles = row[hap_col].split(';')
                all_alleles.update(alleles)
    
    # Extract unique allele identifiers (part after *) for color assignment
    allele_ids = set()
    for allele in all_alleles:
        if '*' in allele:
            allele_ids.add(allele.split('*')[1])
        else:
            allele_ids.add(allele)
    
    # Sort allele identifiers and assign colors
    sorted_allele_ids = sorted(allele_ids)
    allele_id_colors = {allele_id: color_palette[i % len(color_palette)]
                        for i, allele_id in enumerate(sorted_allele_ids)}
    
    # Track which allele identifiers have been added to legend (in sorted order)
    allele_ids_in_legend = []
    
    # Process data for each anchor
    for anchor_idx, anchor in enumerate(anchors):
        col_num = anchor_idx + 1
        hap_col = f'{anchor} hap'
        count_col = f'{anchor} count'
        z_col = f'{anchor} z'
        
        # Collect ALL genes (not just those with haplotype data)
        genes_with_data = []
        for row in rows:
            genes_with_data.append(row['gene'])
        
        if not genes_with_data:
            print(f"No data for anchor {anchor}")
            continue
        
        # Reverse order so first gene appears at top
        genes_with_data = genes_with_data[::-1]
        
        # Create bars for each gene
        for gene in genes_with_data:
            # Find the row for this gene
            row = next(r for r in rows if r['gene'] == gene)
            alleles = row[hap_col].split(';') if row[hap_col] else []
            
            if not alleles:
                # No alleles assigned - create a white bar
                fig.add_trace(
                    go.Bar(
                        x=[1.0],
                        y=[gene],
                        name='No assignment',
                        orientation='h',
                        marker=dict(color='white', line=dict(color='lightgray', width=1)),
                        text='',
                        hovertemplate=f"<b>{gene}</b><br>No alleles assigned<br><extra></extra>",
                        showlegend=False
                    ),
                    row=1, col=col_num
                )
                continue
            
            # Get counts and z-scores for this gene
            # Note: these correspond to ALL alleles in the <anchor> column, not just haplotyped ones
            all_alleles_for_gene = row[anchor].split(';') if row[anchor] else []
            counts = row[count_col].split(';') if row[count_col] else []
            z_scores = row[z_col].split(';') if row[z_col] else []
            
            # Create a mapping from allele name to its count and z-score
            allele_info = {}
            for i, allele_name in enumerate(all_alleles_for_gene):
                allele_info[allele_name] = {
                    'count': counts[i] if i < len(counts) else '0',
                    'z_score': z_scores[i] if i < len(z_scores) else 'N/A'
                }
            
            # Calculate total reads for this anchor/gene (from haplotyped alleles only)
            total_reads = sum(int(allele_info[allele]['count']) for allele in alleles if allele in allele_info)
            
            # Calculate segment width (equal division)
            segment_width = 1.0 / len(alleles)
            
            # Sort alleles by their identifier to maintain consistent order
            alleles_sorted = sorted(alleles, key=lambda x: x.split('*')[1] if '*' in x else x)
            
            # Add a bar segment for each allele
            for allele_idx, allele in enumerate(alleles_sorted):
                # Get count and z-score for this specific allele
                count = allele_info[allele]['count'] if allele in allele_info else '0'
                z_score = allele_info[allele]['z_score'] if allele in allele_info else 'N/A'
                
                # Extract allele identifier (part after *)
                allele_id = allele.split('*')[1] if '*' in allele else allele
                
                # Get consistent color for this allele identifier
                color = allele_id_colors[allele_id]
                
                # Calculate x position (start of this segment)
                x_start = allele_idx * segment_width
                
                # Show in legend only once per allele identifier (in sorted order)
                # Legend order is determined by the order in sorted_allele_ids
                show_legend = allele_id not in allele_ids_in_legend
                if show_legend:
                    allele_ids_in_legend.append(allele_id)
                
                # Find the position in sorted list for legend ordering
                legend_rank = sorted_allele_ids.index(allele_id)
                
                # Build hover text with read counts and z-score
                hover_text = f"<b>{gene}</b><br>{allele}<br>{count}/{total_reads} reads<br>z-score: {z_score}<extra></extra>"
                
                fig.add_trace(
                    go.Bar(
                        x=[segment_width],
                        y=[gene],
                        name=f"*{allele_id}",  # Legend shows just the allele identifier
                        orientation='h',
                        marker=dict(color=color),
                        text=allele_id,  # Show allele number
                        textposition='inside',
                        insidetextanchor='middle',  # Center text horizontally within the bar segment
                        textangle=0,  # Horizontal text (no rotation)
                        textfont=dict(color='white', size=10),
                        base=x_start,
                        hovertemplate=hover_text,
                        showlegend=show_legend,
                        legendgroup=allele_id,
                        legendrank=legend_rank
                    ),
                    row=1, col=col_num
                )
    
    # Update layout
    fig.update_xaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        range=[0, 1]
    )
    
    fig.update_yaxes(
        title_text="Gene",
        tickfont=dict(size=10)
    )
    
    # Extract gene name from CSV filename for title
    gene_name = os.path.basename(csv_file).split('_')[-1].replace('.csv', '')
    
    fig.update_layout(
        title=dict(
            text=f"<b>Haplotype Analysis: {gene_name}</b>",
            x=0.5,
            xanchor='center',
            font=dict(size=16)
        ),
        height=max(400, len(genes_with_data) * 25),  # Dynamic height based on number of genes
        barmode='stack',
        plot_bgcolor='white',
        margin=dict(l=100, r=200, t=100, b=50),
        legend=dict(
            title=dict(text="<b>Alleles</b>"),
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            font=dict(size=10),
            traceorder="normal"
        )
    )
    
    # Determine output filename
    if output_file is None:
        output_file = csv_file.rsplit('.', 1)[0] + '.html'
    
    # Save the plot
    fig.write_html(output_file)
    print(f"Plot saved to {output_file}")
    
    return fig


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('annot_file', type=extant_file, help='Input annotation file')
    parser.add_argument('vdj_ref_file', type=extant_file, help='VDJ reference file')
    parser.add_argument('allele_threshold_table', type=extant_file, help='Allele threshold table file')
    parser.add_argument('report_file', help='Output haplotype report')
    parser.add_argument('--plot', action='store_true', help='Generate Plotly visualization for each haplotype report')
    args = parser.parse_args()

    gene_list = make_gene_list(args.vdj_ref_file)
    gene_counts = scan_for_genes(args.annot_file)

    allele_thresholds = simple.read_csv(args.allele_threshold_table, delimiter='\t')
    allele_thresholds = {row['allele']: float(row['threshold']) for row in allele_thresholds}

    # Qualifying genes must have > 1000 counts and must represent at least 10% of total c_calls
    total_counts = sum(gene_counts.values())
    qualifying_genes = [gene for gene, count in gene_counts.items() if count > 1000 and count / total_counts >= 0.1]
    if qualifying_genes:
        print(f'Qualifying haplotype genes: {qualifying_genes}')
    else:
        print('No qualifying haplotype genes found')
        exit(0)

    for haplotype_gene in qualifying_genes:
        print(f'Processing haplotype gene: {haplotype_gene}')
        fs = args.report_file.rsplit(".", 1)
        name = fs[0]
        ext = fs[1] if len(fs) > 1 else ''
        filename = f'{name}_{haplotype_gene}.{ext}' if ext else f'{name}_{haplotype_gene}'
        haplotype_data = run_haplotype(haplotype_gene, args.annot_file)
        filtered_data = filter_haplotype_data(haplotype_data)

        if len(filtered_data) != 2:
            print('Gene unsuitable for haplotyping. Anchor totals:')
            for c_call, data in haplotype_data.items():
                print(f'{c_call}{"*" if c_call in filtered_data else ""}: {data["total"]}')
            if len(filtered_data) < 2:
                print('Only one dominant anchor after filtering.')
            else:
                print('More than two dominant anchors after filtering.')
            continue

        write_haplotype_report(haplotype_data, list(filtered_data.keys()), gene_list, allele_thresholds, filename)
        
        if args.plot:
            plot_haplotype_data(filename)
