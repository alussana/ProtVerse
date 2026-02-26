#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
download Christopher et al 2025 dataset of LOPIT-DC spatial proteomics
profile changes for A549 cells perturbed with 6 Gy X-ray exposure

both total proteomics abundance and fractionation profiles are provided

Ref <https://doi.org/10.1016/j.mcpro.2024.100888>
*/
process dl_lopit2025 {

    publishDir "${out_dir}",
            pattern: 'databases/lopit2025/1-s2.0-S1535947624001786-mmc1.xlsx',
            mode: 'copy'

    output:
        path 'databases/lopit2025/1-s2.0-S1535947624001786-mmc1.xlsx'
    script:
    """
    mkdir -p databases/lopit2025

    wget \
        -O databases/lopit2025/1-s2.0-S1535947624001786-mmc1.xlsx \
        ${params.url_lopit2025}        
    """

}


/*
build dataset of protein localization changes

there are 3 replicates and 8 fractions (spin speeds) per replicate:
['1000xg',
 '3000xg',
 '5000xg',
 '9000xg',
 '15000xg',
 '79000xg',
 '120000xg',
 'Supernatant']

when considering same replicate and spin speed,
expr and ctrl values have a very high correlation 
PearsonRResult(statistic=0.9354219430500028, pvalue=0.0)
therefore are assumed to be usable for fold change calculation

.META:
1   Accession   A0AV96
2   Replicate   1
3   Spin speed  1000xg
4   Fold Change 1.2345
*/
process lopit2025_prot_loc_changes {

    publishDir "${out_dir}",
            pattern: "databases/lopit2025/lopit2025_prot_loc_changes.tsv",
            mode: 'copy'

    input:
        path "input/input.xlsx"

    output:
        path "databases/lopit2025/lopit2025_prot_loc_changes.tsv"

    script:
    """
    mkdir -p databases/lopit2025

    lopit2025_prot_loc_changes.py \
        input/input.xlsx \
        databases/lopit2025/lopit2025_prot_loc_changes.tsv
    """

}


/*
build list of genes appearing in the dataset
*/
process get_lopit2025_genes {

    publishDir "${out_dir}",
            pattern: "databases/lopit2025/gene_list.txt",
            mode: 'copy'

    input:
        path "input/table.tsv"

    output:
        path "databases/lopit2025/gene_list.txt"

    script:
    """
    mkdir -p databases/lopit2025

    cat input/table.tsv \
        | cut -f1 | sed '1d' | sort | uniq \
        > databases/lopit2025/gene_list.txt
    """
}


/*
Build the partial input vector from lopit2025 data

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2

.META: lopit2025_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2
... lopit2025 features
*/
process lopit2025_featvec {

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/lopit2025_table.tsv')

    output:
        path 'lopit2025_input_vector.tsv'

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi
    
    lopit2025_input_vector.py  \
        input/gene_pairs.tsv \
        input/lopit2025_table.tsv \
        > lopit2025_input_vector.tsv
    """ 

}