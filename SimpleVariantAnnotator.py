import argparse
import os
import sys
import gffutils
from pyfaidx import Fasta
from Bio.Seq import Seq
import re
import pandas as pd
from cyvcf2 import VCF, Writer
import bisect

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = "A simple tool for variant annotation. Based on a sequence (FASTA format) and its annotations (GFF format), it analyzes a list of variants (VCF format), determines the region in which they occur (intergenic/intron/exon), and identifies how they modify the protein sequence.",
        epilog = "\nSimpleVariantAnnotator v0.2.1\nAuthor: Daniel M. Duszczak\nLicense: MIT\n\nIf you find this tool useful, please consider citing it.\nSource code and documentation are available at: https://github.com/maybedannyornot/SVA",
        formatter_class = argparse.RawTextHelpFormatter
    )
    parser.add_argument('--version', action = 'version', version = 'SVA v0.2.1', help = 'Information about current version')
    group = parser.add_argument_group('Basic input file')
    group.add_argument('-v', '--vcf', required = True, help = 'Path to input variant file (VCF)')
    group.add_argument('-f', '--fasta', required = True, help = 'Path to input sequence file (FASTA)')
    group = parser.add_argument_group('Annotation file', description = 'As long as you do not specify the path to a database (--db) previously created from a GFF file, the GFF file (--gff) is absolutely required')
    group.add_argument('-g', '--gff', required = False, help = 'Path to input annotation file (GFF)')
    group.add_argument('--db', default = ':memory:', help = 'Path to the file for writing to/reading from a genome database built from a GFF file (Default: storing the database in RAM)')
    group = parser.add_argument_group('Output settings')
    group.add_argument('-o', '--outdir', required = False, default = 'SVA_results', help = 'Output directory (Default: ./SVA_results/)')
    group.add_argument('-b', '--basename', type = str, help = 'Optional prefix for output files (Default: VCF filename)')
    group.add_argument('--force', action = 'store_true', help = 'Allow overwriting of output files')
    return parser.parse_args()

def gene_intervals(db):
    """Creates a sorted list of genes (based on gffutils db) for quickly searching for intergenic regions"""
    genes = []
    for gene in db.features_of_type('gene'):
        genes.append((gene.chrom, gene.start, gene.end, gene.id, gene.strand))
    return sorted(genes, key = lambda x: (x[0], x[1]))  # Sorted by localization (chromosome and start position)

def get_feature_id(feature,real_number):
    """Attempts to extract a readable exon/CDS name from the attributes"""
    if "Name" in feature.attributes:
        feature_id = feature.attributes['Name'][0]
    elif feature.id:
        feature_id = feature.id
    else:
        feature_id = "Exon"
        
    feature_id = re.sub(r'_\d+$', '', feature_id)
    return f"{feature_id}_{real_number}"

def extract_transcript_info(db, transcript, fasta):
    """Combines CDS sequences, translates to protein and creates a position map"""
    cds_list = list(db.children(transcript, featuretype = 'CDS', order_by='start')) # Extraction of CDSs
    if not cds_list:
        cds_list = list(db.children(transcript, featuretype = 'exon', order_by='start')) # No CDS in annotation, let's try with exons
    if not cds_list:
        return None
    
    strand = cds_list[0].strand
    if strand == '-':
        cds_list = cds_list[::-1] # Reversed exon order
        
    spliced_seq = ""
    genomic_to_local = {}
    
    current = 0 # Current local position base
    exon_number = 1
    exon_boundaries = [] 
    
    for cds in cds_list:
        seq = fasta[cds.chrom][cds.start - 1 : cds.end].seq
        if strand == '-':
            seq = str(Seq(seq).reverse_complement()) # Reverse complement (from Bio.Seq)
        
        spliced_seq += seq
        
        feature_id = get_feature_id(cds, exon_number)
        exon_boundaries.append((cds.start, cds.end, feature_id))
        
        # Creating a position map
        if strand == '+':
            for i in range(cds.start, cds.end + 1):
                genomic_to_local[i] = current + (i - cds.start)
        else:
            for i in range(cds.start, cds.end + 1):
                genomic_to_local[i] = current + (cds.end - i)
                
        current += len(seq)
        exon_number += 1
        
    protein_seq = str(Seq(spliced_seq).translate(to_stop = False))
    
    return {
        'strand': strand,
        'spliced_seq': spliced_seq,
        'protein': protein_seq,
        'mapping': genomic_to_local,
        'exons': sorted(exon_boundaries, key = lambda x: x[0]) # Sorted by start position
    }
    
def format_distance(pos, feature_start, feature_end, strand, is_left):
    """
    Calculates the distance and assigns a sign based on the direction (upstream or downstream):\n
    Upstream (5') = minus (-)\n
    Downstream (3') = plus (+)
    """
    if is_left:
        dist = pos - feature_end
        if dist < 0: dist = 0
        return f'+{dist}' if strand == '+' else f'-{dist}'
    else:
        dist = feature_start - pos
        if dist < 0: dist = 0
        return f'-{dist}' if strand == '+' else f'+{dist}'
    
def get_arrow(strand):
    """Returns a graphical representation of the strand direction"""
    return "->" if strand == '+' else "<-"
    
def get_intergenic_parts(chrom, pos, sorted_genes, fasta):
    """Returns the formatted distance and gene names with arrows"""
    chrom_genes = [gene for gene in sorted_genes if gene[0] == chrom]
    if not chrom_genes:
        if chrom in fasta:
            contig_len = len(fasta[chrom])
            dist_start = pos
            dist_end = contig_len - pos
            return f'(+{dist_start}bp / -{dist_end}bp)', 'Start / End'
        else:
            return '(Unknown / Unknown)', 'Start / End'
    
    starts = [gene[1] for gene in chrom_genes]
    idx = bisect.bisect_right(starts, pos)
    
    if idx > 0:
        left_gene = chrom_genes[idx-1]
        l_name = left_gene[3]
        l_strand = left_gene[4]
        l_dist = format_distance(pos, left_gene[1], left_gene[2], l_strand, True)
        l_part = f'{l_name} {get_arrow(l_strand)}'
    else:
        l_dist = f'+{pos}'
        l_part = 'Start'
    
    if idx < len(chrom_genes):
        right_gene = chrom_genes[idx]
        r_name = right_gene[3]
        r_strand = right_gene[4]
        r_dist = format_distance(pos, right_gene[1], right_gene[2], r_strand, False)
        r_part = f'{get_arrow(r_strand)} {r_name}'
    else:
        if chrom in fasta:
            r_dist = f'-{len(fasta[chrom]) - pos}'
        else:
            r_dist = 'Unknown'
        r_part = 'End'
        
    return f'({l_dist}bp / {r_dist}bp)', f'{l_part} / {r_part}'

def main():
    # Parsing arguments and setting output filenames
    arguments = parse_arguments()
    
    db_path = arguments.db
    gff_path = arguments.gff
    if db_path == ':memory:' or not os.path.exists(db_path):
        if not gff_path:
            print('\nERROR! Missing GFF file!')
            if db_path == ':memory:':
                print('If you want to store database in RAM, you must specify the GFF file.')
            else:
                print(f'To generate file {db_path}, you must specify the GFF file.')
            sys.exit(1)
    
    fasta_path = arguments.fasta
    vcf_path = arguments.vcf
    output_dir = arguments.outdir
    
    os.makedirs(output_dir, exist_ok = True)
    
    if arguments.basename:
        basename_prefix = arguments.basename
    else:
        basename_prefix = os.path.basename(vcf_path).replace('.vcf.gz', '').replace('.vcf', '')
        
    out_proteins = os.path.join(output_dir, f'{basename_prefix}_proteins.fasta')
    out_genes_report = os.path.join(output_dir, f'{basename_prefix}_genes_report.csv')
    out_report = os.path.join(output_dir, f'{basename_prefix}_variants_report.csv')
    out_summary = os.path.join(output_dir, f'{basename_prefix}_short_summary.csv')
    out_vcf = os.path.join(output_dir, f'{basename_prefix}_annotated.vcf.gz')
    
    # Overwrite check
    output_files = [out_proteins, out_genes_report, out_report, out_summary, out_vcf]
    
    if not arguments.force:
        for path in output_files:
            if os.path.exists(path):
                print(f'\nERROR! Output file already exists: {path}')
                print("Use --force to overwrite files or change filenames using --basename.")
                sys.exit(1)
                
    
    # Building genome database from GFF 
    if db_path == ':memory:':
        print('Building a database from GFF file in RAM...')
        db = gffutils.create_db(gff_path, dbfn = db_path, force = True, keep_order = True, merge_strategy = 'merge')
    else:
        if os.path.exists(db_path):
            print(f'Loading existing genome database from file {db_path}...')
            db = gffutils.FeatureDB(db_path, keep_order = True)
        else:
            print('Building a database from GFF file...')
            db = gffutils.create_db(gff_path, dbfn = db_path, force = True, keep_order = True, merge_strategy = 'merge')
            print(f'Database saved to file {db_path}')
    
    # Creating a sorted list of genes
    sorted_genes = gene_intervals(db)
    
    # Reading FASTA file
    print("Reading a FASTA file...")
    fasta = Fasta(fasta_path)
    
    # Genes report and proteins translation
    print("Generating a gene raport, exon splicing, and protein translation...")
    transcripts_data = {}
    genes_report = []
    
    with open(out_proteins, 'w') as file_prot:
        for gene in db.features_of_type('gene'):
            gene_length = gene.end - gene.start + 1
            
            transcripts = list(db.children(gene, featuretype = ['mRNA', 'transcript']))
            if not transcripts:
                genes_report.append({
                    'Gene': gene.id,
                    'Total_Length_bp': gene_length,
                    'Exon count': 0,
                    'Exon_Lengths_bp': "-"
                })
                continue
            
            transcript = transcripts[0]
            t_info = extract_transcript_info(db, transcript, fasta)
            
            if t_info:
                transcripts_data[transcript.id] = t_info
                file_prot.write(f'>{gene.id}_{transcript.id}\n{t_info["protein"]}\n')
                
                exons = t_info['exons']
                exon_count = len(exons)
                exon_lengths = [str(abs(exon[1] - exon[0]) + 1) for exon in exons]
                if t_info['strand'] == '-':
                    exon_lengths.reverse()
                
                genes_report.append({
                    'Gene': gene.id,
                    'Total_Length_bp': gene_length,
                    'Exon count': exon_count,
                    'Exon_Lengths_bp': ', '.join(exon_lengths)
                })
            else:
                genes_report.append({
                    'Gene': gene.id,
                    'Total_Length_bp': gene_length,
                    'Exon count': 0,
                    'Exon_Lengths_bp': "-"
                })
    print(f'Protein sequences saved to file {out_proteins}')
    pd.DataFrame(genes_report).to_csv(out_genes_report, index = False)
    print(f'Gene report saved to file {out_genes_report}')
    
    print("Reading VCF file...")
    vcf = VCF(vcf_path)
    
    print("Analyzing variants...")
    # Preprocessing of variants
    mutations = {}
    
    for variant in vcf:
        chrom, pos = variant.CHROM, variant.POS
        alt = variant.ALT[0]
        
        overlapping_genes = list(db.region(seqid = chrom, start = pos, end = pos, featuretype = 'gene')) # Where
        if overlapping_genes:
            gene = overlapping_genes[0]
            transcripts = list(db.children(gene, featuretype = ['mRNA', 'transcript']))
            if transcripts and transcripts[0].id in transcripts_data:
                t_id = transcripts[0].id
                t_data = transcripts_data[t_id]
                if pos in t_data['mapping']:
                    t_pos = t_data['mapping'][pos]
                    strand = t_data['strand']
                    eff_alt = alt if strand == '+' else str(Seq(alt).complement())
                    
                    if t_id not in mutations:
                        mutations[t_id] = {}
                    mutations[t_id][t_pos] = eff_alt # Save the mutation under a specific index
    vcf.close()
    # Annotation and report generation
    print("Creating annotation and reports...")
    vcf = VCF(vcf_path)
    vcf.add_info_to_header({'ID': 'MUT_INFO', 'Description': 'Custom Annotation: Type, Gene(s), Detail', 'Type': 'String', 'Number': '1'})
    out_vcf_file = Writer(out_vcf, vcf)
    
    report = []
    summary = {}
    global_count = 0
    intergenic_count = 0
    
    for variant in vcf:
        global_count += 1
        
        chrom = variant.CHROM
        pos = variant.POS
        ref = variant.REF
        alt = variant.ALT[0]
        mut_info = ''
        
        overlapping_genes = list(db.region(seqid = chrom, start = pos, end = pos, featuretype = 'gene')) # Where
        
        if not overlapping_genes: # Intergenic
            intergenic_count += 1
            dist, genes = get_intergenic_parts(chrom, pos, sorted_genes, fasta)
            detail = f'{ref}>{alt}'
            mut_info = f'Intergenic {dist}, {genes}'
            
            report.append({
                'Chrom': chrom,
                'Pos': pos,
                'Type': f'Intergenic {dist}',
                'Gene': genes,
                'Mutation': detail
            })
            
        else: # Intragenic
            gene = overlapping_genes[0]
            gene_id = gene.id
            
            if gene_id not in summary:
                summary[gene_id] = {'total': 0, 'exonic': 0, 'intronic': 0, 'syn': 0, 'nonsyn': 0, 'frameshift': 0}
            summary[gene_id]['total'] += 1
            
            transcripts = list(db.children(gene, featuretype = ['mRNA', 'transcript']))
            if not transcripts or transcripts[0].id not in transcripts_data:
                continue
            
            t_id = transcripts[0].id
            t_data = transcripts_data[t_id]
            strand = t_data['strand']
            arrow = get_arrow(strand)
            
            if pos in t_data['mapping']: # Position is mapped -> Exon
                summary[gene_id]['exonic'] += 1
                
                t_pos = t_data['mapping'][pos]
                codon_start = t_pos - (t_pos % 3)
                
                origin_codon = list(t_data['spliced_seq'][codon_start:codon_start+3])
                mut_codon = list(origin_codon)
                
                # Looking for all mutations in the codon
                for i in range(3):
                    current_pos = codon_start + i
                    if current_pos in mutations.get(t_id, {}):
                        mut_codon[i] = mutations[t_id][current_pos]
                        
                origin_codon = "".join(origin_codon)
                mut_codon = "".join(mut_codon)
                
                detail = f'{ref}>{alt}'
                if len(origin_codon) == 3 and len(mut_codon) == 3:
                    origin_aa = str(Seq(origin_codon).translate())
                    mut_aa = str(Seq(mut_codon).translate())
                    aa_pos = (codon_start // 3) + 1
                    detail2 = f'{origin_aa}{aa_pos}{mut_aa} ({origin_codon}>{mut_codon})'
                    if origin_aa == mut_aa:
                        summary[gene_id]['syn'] += 1
                    else:
                        summary[gene_id]['nonsyn'] += 1
                else:
                    detail2 = 'Frameshift'
                    summary[gene_id]['frameshift'] += 1
                
                combined_detail = detail + ', ' + detail2
                gene = f'{gene_id} {arrow}'
                mut_info = f'Exon, {gene}, {detail2}'
                
                report.append({
                    'Chrom': chrom,
                    'Pos': pos,
                    'Type': f'Exon',
                    'Gene': gene,
                    'Mutation': combined_detail
                })
                
            else: # Position is not mapped -> Intron
                summary[gene_id]['intronic'] += 1
                exons = t_data['exons']
                
                for i in range(len(exons) - 1):
                    if exons[i][1] < pos < exons[i+1][0]:
                        left_exon = exons[i]
                        right_exon = exons[i+1]
                        
                        l_dist = format_distance(pos, left_exon[0], left_exon[1], strand, True)
                        r_dist = format_distance(pos, right_exon[0], right_exon[1], strand, False)
                        
                        dist = f'({l_dist}bp / {r_dist}bp)'
                        genes = f'{left_exon[2]} {arrow} / {arrow} {right_exon[2]}'
                        break
                
                detail = f'{ref}>{alt}'
                mut_info = f'Intron {dist}, {genes}'
                
                report.append({
                    'Chrom': chrom,
                    'Pos': pos,
                    'Type': f'Intron {dist}',
                    'Gene': genes,
                    'Mutation': detail
                })
                
        variant.INFO['MUT_INFO'] = mut_info
        out_vcf_file.write_record(variant)
    
    out_vcf_file.close()
    vcf.close()
    print(f'Annotated variants saved to file {out_vcf}')
    
    print('Saving additional report files...')
    df_report = pd.DataFrame(report)
    df_report.to_csv(out_report, index = False)
    print(f'Mutation report saved to file {out_report}')
    
    summary_list = []
    global_stats = {'total': 0, 'exonic': 0, 'intronic': 0, 'syn': 0, 'nonsyn': 0, 'frameshift': 0}
    for gene, stats in summary.items():
        stats['Gene'] = gene # Adds key name to data
        summary_list.append(stats)
        for k in global_stats: # Adds gene stats to global stats
            global_stats[k] += stats[k]
            
    df_summary = pd.DataFrame(summary_list)
    if not df_summary.empty:
        df_summary = df_summary[['Gene','total','exonic','intronic','syn','nonsyn','frameshift']]
    df_summary.to_csv(out_summary, index = False)
    print(f'Summary report saved to file {out_summary}')
    print('Generating a final summary...')
    print('\n --- A FINAL SUMMARY ---')
    print(f'Number of identified mutations: {global_count}')
    print(f' Intergenic: {intergenic_count}')
    print(f' Intragenic: {global_stats["total"]}')
    print(f'    in CDSs: {global_stats["exonic"]}  ({global_stats["frameshift"]} frameshift(s))')
    print(f'        Synonymous: {global_stats["syn"]}')
    print(f'        Non-synonymous: {global_stats["nonsyn"]}')
    print(f'    in introns: {global_stats["intronic"]}')
    
if __name__ == '__main__':
    main()