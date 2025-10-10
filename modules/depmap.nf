#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Get gene dependencies fold changes computed in
<https://doi.org/10.1038/s41467-021-21898-7>

908 cell lines (columns) (model_id from <https://cellmodelpassports.sanger.ac.uk/>)
17486 genes (rows) (HUGO gene symbols)
*/
process download_essentiality_matrices {

    publishDir "${out_dir}",
                pattern: "databases/depMap/*txt",
                mode: 'copy'

    output:
        path 'databases/depMap/CERES_FC.txt', 
            emit: CERES_FC
        path 'databases/depMap/CRISPRcleanR_FC.txt',
            emit: CRISPRcleanR_FC

    script:
    """
    mkdir -p databases/depMap/
    wget --user-agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36" \\
         --referer="https://example.com/" \\
         --content-on-error --content-disposition -O databases/depMap/26705222.zip \\
         "${params.url_depMap}"

    unzip databases/depMap/26705222.zip

    mv integrated_Sanger_Broad_essentiality_matrices_20201201/CERES_FC.txt \\
        databases/depMap/CERES_FC.txt

    mv integrated_Sanger_Broad_essentiality_matrices_20201201/CRISPRcleanR_FC.txt \\
        databases/depMap/CRISPRcleanR_FC.txt

    rm -fr integrated_Sanger_Broad_essentiality_matrices_20201201
    """

}

/*
Get the list of gene names in the essentiality matrices to be translated
*/
process get_depmap_genes {

    publishDir "${out_dir}",
                pattern: "databases/depMap/genes.txt",
                mode: 'copy'

    input:
        path 'input/gene_dependencies.txt'

    output:
        path 'databases/depMap/genes.txt'

    script:
    """
    mkdir -p databases/depMap
    
    cat input/gene_dependencies.txt | cut -f1 | awk 'NF' \
        | sort | uniq \
        > databases/depMap/genes.txt
    """

}

/*
Build the partial input vector from gene dependency data

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
... depmap features
*/
process dependency_featvec {

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/dependency_table.tsv')

    output:
        path 'dependency_input_vector.tsv'

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi
    
    dependency_input_vector.py  \
        input/gene_pairs.tsv \
        input/dependency_table.tsv \
        > dependency_input_vector.tsv
    """ 

}