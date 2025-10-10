#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
download Mitchell et al 2022 dataset of proteome-wide profile changes for 875
compounds in a human cancer cell line, measurement taken
after 24h from perturbation

Ref <https://doi.org/10.1038/s41587-022-01539-0>
*/
process dl_mitchell2023 {

    publishDir "${out_dir}",
            pattern: 'databases/mitchell2023/41587_2022_1539_MOESM3_ESM.csv',
            mode: 'copy'

    output:
        path 'databases/mitchell2023/41587_2022_1539_MOESM3_ESM.csv'

    script:
    """
    mkdir -p databases/mitchell2023

    wget \
        -O databases/mitchell2023/41587_2022_1539_MOESM3_ESM.csv \
        ${params.url_mitchell2023}        
    """

}

/*
get the list of gene names to be translated
*/
process get_mitchell2023_genes {

    publishDir "${out_dir}",
            pattern: "databases/mitchell2023/genes.txt",
            mode: 'copy'

    input:
        path 'input/matrix.tsv'

    output:
        path 'databases/mitchell2023/genes.txt'

    script:
    """
    mkdir -p databases/mitchell2023

    cat input/matrix.tsv \
        | cut -f1 \
        | sed '1d' | awk 'NF' \
        | sort | uniq \
        > databases/mitchell2023/genes.txt
    """

}

/*
keep only genes having a defined Gene Name

NOTE: some gene names are duplicated due to different protein isoforms
      9714 unique gene names, 9950 rows

.META: 
1   Gene_names
... perturbation
*/
process parse_mitchell2023 {

    publishDir "${out_dir}",
                pattern: "databases/mitchell2023/matrix.tsv",
                mode: 'copy'

    input:
        path 'input/table_S1.csv'

    output:
        path 'databases/mitchell2023/matrix.tsv'

    script:
    """
    mkdir -p databases/mitchell2023

    parse_mitchell2023.py \
        input/table_S1.csv \
        > databases/mitchell2023/matrix.tsv
    """

}

/*
Build the partial input vector from mitchell2023 data

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2

.META: dependency_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2
... mitchell2023 features
*/
process mitchell2023_featvec {

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/mitchell2023_table.tsv')

    output:
        path 'mitchell2023_input_vector.tsv'

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi
    
    mitchell2023_input_vector.py  \
        input/gene_pairs.tsv \
        input/mitchell2023_table.tsv \
        > mitchell2023_input_vector.tsv
    """ 

}
