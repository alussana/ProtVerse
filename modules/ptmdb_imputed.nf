#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
PTMDB phosphosites across samples were imputed using a variational
autoencoder with SCVI <https://scvi-tools.org>
*/
process parse_ptmdb_imputed {

    publishDir "${out_dir}", 
                pattern: "databases/ptmdb_imputed/ptmdb_imputed.tsv",
                mode: 'copy'

    input:
        path 'input/ptmdb_imputed.tsv'

    output:
        path "databases/ptmdb_imputed/ptmdb_imputed.tsv"

    script:
    """
    mkdir -p databases/ptmdb_imputed

    ptmdb_imputed_phos.py \
        input/ptmdb_imputed.tsv \
        databases/ptmdb_imputed/ptmdb_imputed.tsv    
    """

}

/*
Get the list of gene names in the ptmdb to be translated
*/
process get_ptmdb_imputed_genes {

    publishDir "${out_dir}",
                pattern: 'databases/ptmdb_imputed/genes.txt',
                mode: 'copy'

    input:
        path 'input/ptmdb_matrix.tsv'

    output:
        path 'databases/ptmdb_imputed/genes.txt'

    script:
    """
    mkdir -p databases/ptmdb_imputed

    cat input/ptmdb_matrix.tsv \
        | cut -f1 | awk 'NF' \
        | sort | uniq \
        > databases/ptmdb_imputed/genes.txt
    """

}

/*
Build the partial input vector from PTMDB data
*/
process ptmdb_imputed_featvec {

    input:
        tuple val(gene_pairs),
              file('input/ptmdb_imputed.tsv')

    output:
        path 'ptmdb_imputed_input_vector.tsv'

    script:
    """
    echo -e "${gene_pairs}" | awk 'NF' > gene_pairs.tsv

    ptmdb_imputed_input_vector.py  \
        gene_pairs.tsv \
        input/ptmdb_imputed.tsv \
        > ptmdb_imputed_input_vector.tsv
    """ 

}