#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Download hu.MAP 3.0 data, containing probability scores for
25992007 protein pairs

.META:
1   gene_1 UniProtAC
2   gene_2 UniProtAC
3   score 
*/
process download_humap3_net {

    publishDir "${out_dir}",
                pattern: 'databases/huMAP3/source.tsv.gz',
                mode: 'copy'

    output:
        path 'databases/huMAP3/source.tsv.gz'

    script:
    """
    mkdir -p databases/huMAP3
    wget -O databases/huMAP3/source.tsv.gz -P databases/huMAP3 ${params.url_humap3_network}
    """

}

/*
filter for humap probability score > ${params.humap3_min_score}

round all probability scores to 4 decimal digits, save in tsv format

.META:
1   gene_a
2   gene_b
4   prob 
*/
process humap3_table {

    publishDir "${out_dir}",
                pattern: 'databases/huMAP3/huMAP3_table.tsv',
                mode: 'copy'

    input:
        path 'input/net.tsv.gz'

    output:
        path 'databases/huMAP3/huMAP3_table.tsv'

    script:
    """
    mkdir -p databases/huMAP3

    zcat input/net.tsv.gz | awk '\$3>${params.humap3_min_score}' \
        > humap_filtered.tsv

    humap.py \
        humap_filtered.tsv \
        > databases/huMAP3/huMAP3_table.tsv
    """

}

/*
List all genes found in the network
*/
process get_humap3_genes {

    publishDir "${out_dir}",
                pattern: 'databases/huMAP3/genes.txt',
                mode: 'copy'

    input:
        path 'input/net.tsv.gz'

    output:
        path 'databases/huMAP3/genes.txt'

    script:
    """
    mkdir -p databases/huMAP3

    zcat input/net.tsv.gz | cut -f1,2 | tr '\t' '\n' \
        | sort | uniq \
        > databases/huMAP3/genes.txt
    """
}

/*
Build the partial input vector from hu.MAP data (probability score)

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2

.META: humap_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2
5   humap score
*/
process humap3_score {

    memory "16G"

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/humap_table.tsv')

    output:
        path 'humap_input_vector.tsv'

    script:
    """
    humap3_input_vector.py  \
        input/gene_pairs.tsv \
        input/humap_table.tsv \
        > humap_input_vector.tsv
    """ 

}