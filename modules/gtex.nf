#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Download median transcripts per million in GTEx tissues
56200 ENSG (rows)
54 tissues (columns)
*/
process download_median_tpm {

    publishDir "${out_dir}",
                pattern: "databases/gtex/*gz",
                mode: 'copy'

    output:
        path 'databases/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz'

    script:
    """
    mkdir -p databases/gtex
    wget -P databases/gtex ${params.url_GTEx_medianTPM}
    """

}

/*
Generate a genes (rows) x samples (columns) TPM matrix
Note: In the ENSG names the version (/\..*$/) is truncated
*/
process parse_median_tpm {

    publishDir "${out_dir}",
                pattern: "databases/gtex/median_tpm.tsv",
                mode: 'copy'

    input:
        path 'input/gtex_median_tpm,gct.gz'

    output:
        path 'databases/gtex/median_tpm.tsv'

    script:
    """
    mkdir -p databases/gtex
    zcat input/gtex_median_tpm,gct.gz | sed '1,2d' | cut -f1,3- \
        | sed -r 's/(ENSG[0-9][0-9]*)\\.[0-9][0-9]*/\\1/' \
        > databases/gtex/median_tpm.tsv
    """

}

/*
Get the list of gene names in the median TPM matrix to be translated
*/
process get_gtex_genes {

    publishDir "${out_dir}",
                pattern: "databases/gtex/genes.txt",
                mode: 'copy'

    input:
        path 'input/median_tpm.tsv'

    output:
        path 'databases/gtex/genes.txt'

    script:
    """
    mkdir -p databases/gtex
    
    cat input/median_tpm.tsv | cut -f1 | sed '1d' | awk 'NF' \
        | sort | uniq \
        > databases/gtex/genes.txt
    """

}

/*
Build the partial input vector from GTEx data

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2

.META: gtex_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2
... gtex features
*/
process gtex_featvec {

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/gene_tpm.tsv')

    output:
        path 'gtex_input_vector.tsv'

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi
    
    gtex_input_vector.py  \
        input/gene_pairs.tsv \
        input/gene_tpm.tsv \
        > gtex_input_vector.tsv
    """ 

}