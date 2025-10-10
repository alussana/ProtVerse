#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
.META:
1   module name (int starting from 0)
2   tab-separated list of genes
*/
process run_blant {

    container "${projectDir}/envs/blant-clusters.sif"
    cpus "${params.clusterone_n_cores}"
    memory '20G'

    publishDir "${out_dir}", pattern: "modules/${net}/blant/*", mode: 'copy'

    input:
        tuple val(net),
              file('input/edges.tsv'),
              val(min_weight)

    output:
        tuple val("${net}/blant"),
              file("modules/${net}/blant/clusters_minW${min_weight}.tsv"),
              val(min_weight)

    script:
    """
    mkdir -p modules/${net}/blant

    cat input/edges.tsv \\
        | awk '\$3>${min_weight}' \\
        > edges_minW${min_weight}.elw

    /scripts/blant-clusters.sh \\
        -1 /blant 4 ${min_weight} edges_minW${min_weight}.elw \\
        > modules/${net}/blant/clusters_minW${min_weight}.tsv
    """

}