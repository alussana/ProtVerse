#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Get human genes phylogenetic profiling data from
<https://doi.org/10.1016/j.celrep.2015.01.025>

Information for 31406 orthogroups across 177 species including 19974 human
protein-coding genes
*/
process download_orthogroups {

    publishDir "${out_dir}",
                pattern: 'databases/orthogroups/hOP_occurrence_matrix.tsv',
                mode: 'copy'

    publishDir "${out_dir}",
                pattern: 'databases/orthogroups/gene_orthogroup_key.tsv',
                mode: 'copy'

    publishDir "${out_dir}",
                pattern: 'databases/orthogroups/PCS_hOPMatrix.tsv',
                mode: 'copy'

    output:
        path 'databases/orthogroups/hOP_occurrence_matrix.tsv', emit: occurrence_matrix
        path 'databases/orthogroups/gene_orthogroup_key.tsv' , emit: gene2orthogroup
        path 'databases/orthogroups/PCS_hOPMatrix.tsv', emit: PCS

    script:
    """
    mkdir -p databases/orthogroups/

    wget -O databases/orthogroups/1-s2.0-S2211124715000509-mmc6.zip \
        ${params.url_orthogroups}

    unzip -o databases/orthogroups/1-s2.0-S2211124715000509-mmc6.zip

    cat \
        <(paste \
            <(echo "Orthogroup_Label" \
                ) \
            <(cat File3_Species_Names_Ordered.txt \
                | cut -f2 | sed '1d' | transpose.awk \
                | perl -p -i -e "s/\r//g" \
                ) \
            ) \
        <(paste \
            <(cat File2_Orthogroup_Labels.txt | cut -f2 | sed '1d' \
                | perl -p -i -e "s/\r//g" \
                ) \
            <(cat File4_hOPMAP.txt | tr ' ' '\t' \
                ) \
            ) \
        > databases/orthogroups/hOP_occurrence_matrix.tsv

    paste \
        <(cat File1_Gene_Orthogroup_Key.txt | cut -f1-3) \
        <(cat File1_Gene_Orthogroup_Key.txt | cut -f4- | tr "\\t" "," | sed 's/\\(.*\\),/\\1/') \
        > databases/orthogroups/gene_orthogroup_key.tsv

    mv File5_PCS_hOPMatrix.txt databases/orthogroups/PCS_hOPMatrix.tsv
    """

}

/*
Get the list of reported gene names to be translated
*/
process get_orthogroups_genes {

    publishDir "${out_dir}",
                pattern: 'databases/orthogroups/genes.txt',
                mode: 'copy'

    input:
        path 'input/gene2orthogroup.tsv'

    output:
        path 'databases/orthogroups/genes.txt'

    script:
    """
    mkdir -p databases/orthogroups

    cat input/gene2orthogroup.tsv \
        | cut -f1 | sed '1d' | awk 'NF' \
        | sort | uniq \
        > databases/orthogroups/genes.txt
    """

}

/*
Build the partial input vector from orthogroup data

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2

.META: orthogroup_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2
5   highest orthogroup pcs for the gene pair
*/
process orthogroup_featvec {

    debug "false"

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/orthogroup_table.tsv'),
              file('input/orthogroup_pcs.tsv')

    output:
        path 'orthogroup_input_vector.tsv'

    script:
    """
    orthogroup_input_vector.py  \
        input/gene_pairs.tsv \
        input/orthogroup_table.tsv \
        input/orthogroup_pcs.tsv \
        > orthogroup_input_vector.tsv
    """ 

}