#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
download 294 SILAC ratios for 10,323 proteins that make up ProteomeHD

<https://www.proteomehd.net/download>

.META: databases/proteomeHD/table_S1.csv, 10323 proteins x 294 perturbations
1  Majority_protein_IDs
2  Simplified_protein_ID
3  Protein_names
4  Gene_names
5  RatioHL_${perturbation_name} 
*/
process dl_proteomehd {

    publishDir "${out_dir}",
                pattern: "databases/proteomeHD/table_S1.csv",
                mode: 'copy'

    output:
        path 'databases/proteomeHD/table_S1.csv'

    script:
    """
    mkdir -p databases/proteomeHD

    wget -O databases/proteomeHD/table_S1.csv '${params.url_proteomeHD}'
    """

}

/*
keep only genes having a defined Gene Name

when multiple gene names are associated to the same record,
explode the record (e.g. 3 gene names will result in a row duplicated 3 times,
each mapped to a different individual gene name)

.META: 
1  Gene_names
2  RatioHL_${perturbation_name}
*/
process parse_proteomehd {

    publishDir "${out_dir}",
                pattern: "databases/proteomeHD/SILAC_ratios.tsv",
                mode: 'copy'

    input:
        path 'input/table_S1.csv'

    output:
        path 'databases/proteomeHD/SILAC_ratios.tsv'

    script:
    """
    mkdir -p databases/proteomeHD

    parse_proteomehd.py \
        input/table_S1.csv \
        > databases/proteomeHD/SILAC_ratios.tsv
    """

}

/*
get the list of gene names in SILAC ratios matrix to be translated
*/
process get_proteomehd_genes {

    publishDir "${out_dir}",
            pattern: "databases/proteomeHD/genes.txt",
            mode: 'copy'

    input:
        path 'input/SILAC_ratios.tsv'

    output:
        path 'databases/proteomeHD/genes.txt'

    script:
    """
    mkdir -p databases/proteomeHD

    cat input/SILAC_ratios.tsv \
        | cut -f1 | sed '1d' | awk 'NF' | sort | uniq \
        > databases/proteomeHD/genes.txt
    """

}

/*
Build the partial input vector from proteomeHD data

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
... proteomehd features
*/
process proteomehd_featvec {

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/proteomehd_table.tsv')

    output:
        path 'proteomehd_input_vector.tsv'

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi
    
    proteomehd_input_vector.py  \
        input/gene_pairs.tsv \
        input/proteomehd_table.tsv \
        > proteomehd_input_vector.tsv
    """ 

}