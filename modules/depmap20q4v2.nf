#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
download Achilles_gene_effect.csv version 20Q4v2
from the DepMap portal <https://depmap.org/portal/download/all/>

it consists of CERES scores for 789 cell lines x 18119 genes, where the
principle components strongly related to known batch effects removed,
then shifted and scaled per cell line so the median nonessential KO effect
is 0 and the median essential KO effect is -1
*/
process dl_depmap20q4v2 {

    publishDir "${out_dir}",
                pattern: "databases/DepMap20Q4v2/25770029",
                mode: 'copy'

    output:
        path 'databases/DepMap20Q4v2/25770029' 
    
    script:
    """
    mkdir -p databases/DepMap20Q4v2

    wget -P databases/DepMap20Q4v2/ \
        ${params.url_DepMap20Q4v2_scaled}
    """
}

/*
parse Achilles_gene_effect.csv version 20Q4v2
from the DepMap portal <https://depmap.org/portal/download/all/>
making in a tab-separated file with a header of HUGO genes and
index (i.e. DepMap_ID of cell lines) on the first column
*/
process parse_depmap20q4v2 {

    publishDir "${out_dir}",
                pattern: "databases/DepMap20Q4v2/Achilles_gene_effect.tsv",
                mode: 'copy'

    input:
        path 'input/Achilles_gene_effect.csv'

    output:
        path 'databases/DepMap20Q4v2/Achilles_gene_effect.tsv' 
    
    script:
    """
    mkdir -p databases/DepMap20Q4v2

    cat \
        <(cat input/Achilles_gene_effect.csv \
            | sed -n '1p' | tr ',' '\\t' | sed -r 's/ \\([0-9]+\\)//g' \
        ) \
        <(cat input/Achilles_gene_effect.csv \
            | sed '1d' | tr ',' '\\t' \
        ) \
    > databases/DepMap20Q4v2/Achilles_gene_effect.tsv   
    """
}

/*
get genes names reported in Achilles_gene_effect.csv version 20Q4v2
from the DepMap portal <https://depmap.org/portal/download/all/>
*/
process get_depmap20q4v2_genes {

    input:
        path 'input/Achilles_gene_effect.tsv'

    output:
        path 'genes.txt'

    script:
    """
    cat input/Achilles_gene_effect.tsv \
        | sed -n '1p' | tr '\\t' '\\n' | sed '1d' \
        >  genes.txt
    """

}

/*
process the Ceres scores as described in Gheorghe and Hart 2022
<https://doi.org/10.1101/2022.08.03.502694>

perform covariance normalization (whitening) using PCA

plot PC1 and PC2 for the PCA-transformed matrix, w/ and w/o whitening

NOTE: in this step cell lines with missing data are discarded,
      and as a result 730 cell lines are kept
*/
process covar_norm {

    publishDir "${out_dir}",
                pattern: "databases/DepMap20Q4v2/*",
                mode: 'copy'

    input:
        path 'input/Achilles_gene_effect.tsv'

    output:
        path 'databases/DepMap20Q4v2/Achilles_gene_effect_white.tsv', emit: X_white
        path 'databases/DepMap20Q4v2/*.png'
    
    script:
    """
    mkdir -p databases/DepMap20Q4v2

    depmap20q4v2_covar_norm.py \
        input/Achilles_gene_effect.tsv \
        databases/DepMap20Q4v2/Achilles_gene_effect_pcs.png \
        databases/DepMap20Q4v2/Achilles_gene_effect_white_pcs.png \
        databases/DepMap20Q4v2/Achilles_gene_effect_white.tsv
    """

}