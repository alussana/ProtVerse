#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
apply disparity filter to prune ProtVerse given a maximum alpha

p_{ij} = \frac{w_{ij}}{s_{i}}

a_{ij} = (1 - p_{ij})^{(k - 1)}

ref:
M. Á. Serranoa , M. Boguñáb, and A. Vespignani (2009)
"Extracting the multiscale backbone of complex weighted networks"
<https://doi.org/10.1073/pnas.0808904106>
*/
process disparity_filter_prune {

    memory "${params.adj_dask_memory}"
    cpus "${params.adj_dask_n_proc}"

    input:
        path 'input/edges.parquet.tar.gz'

    output:
        path "pruned_${params.disparity_filter_alpha}.tsv.gz"

    script:
    """
    tar -zxvf input/edges.parquet.tar.gz

    disparity_filter.py \\
        edges.parquet \\
        pruned_${params.disparity_filter_alpha}.tsv \\
        ${params.disparity_filter_alpha}

    gzip pruned_${params.disparity_filter_alpha}.tsv \\
    """

}


/*
plot:
* degree distribution on histogram with logarithmic axis
* degree complementary cumulative density function
*/
process plot_parquet_net_stats {

    memory "${params.adj_dask_memory}"
    cpus "${params.adj_dask_n_proc}"

	publishDir "${out_dir}",
        pattern: "${id}/*.pdf",
        mode: 'copy'

    input:
        path 'input/edges.parquet.tar.gz'
        val id

    output:
         path "${id}/*.pdf"

    script:
    """
    mkdir -p ${id} edges.parquet

    tar \\
        -zxvf input/edges.parquet.tar.gz \\
        -C edges.parquet \\
        --strip-components=1

    plot_parquet_net_stats.py \\
        edges.parquet \\
        ${id}/pruned_${params.disparity_filter_alpha}_
    """

}


/*
TODO

plot:
* degree distribution on histogram with logarithmic axis
* degree complementary cumulative density function
*/
process plot_tsv_net_stats {

    memory "${params.adj_dask_memory}"
    cpus "${params.adj_dask_n_proc}"

	publishDir "${out_dir}",
        pattern: "${id}/*.pdf",
        mode: 'copy'

    input:
        path 'input/edges.parquet.tar.gz'
        val id

    output:
         path "${id}/*.pdf"

    script:
    """
    mkdir -p ${id} edges.parquet

    tar \\
        -zxvf input/edges.parquet.tar.gz \\
        -C edges.parquet \\
        --strip-components=1

    plot_tsv_net_stats.py \\
        edges.parquet \\
        ${id}/pruned_${params.disparity_filter_alpha}_
    """

}


/*
NOT USED

merge two parquet graphs
*/
process cat_pq_graphs {

    input:
        path 'input/edges1.parquet.tar.gz'
        path 'input/edges2.parquet.tar.gz'

    output:
        path "ProtVerse/edge_alpha${params.disparity_filter_alpha}.parquet.tar.gz"

    script:
    """
    mkdir -p ProtVerse

    tar -zxOf input/edges1.parquet.tar.gz > edges1.parquet
    tar -zxOf input/edges2.parquet.tar.gz > edges2.parquet

    concat_parquet_graphs.py \\
        edges1.parquet \\
        edges2.parquet \\
        ProtVerse/edges_alpha${params.disparity_filter_alpha}.parquet

    tar -zcvf \\
        ProtVerse/edges_${params.disparity_filter_alpha}.parquet.tar.gz \\
        ProtVerse/edges_${params.disparity_filter_alpha}.parquet
    """

}


/*
merge two tsv graphs
*/
process cat_tsv_graphs {

    publishDir "${out_dir}",
        pattern: "ProtVerse/*tsv",
        mode: 'copy'

    input:
        path 'input/edges1.tsv.gz'
        path 'input/edges2.tsv.gz'

    output:
        path "ProtVerse/edges_alpha${params.disparity_filter_alpha}.tsv.gz"

    script:
    """
    mkdir -p ProtVerse

    concat_tsv_graphs.py \\
        input/edges1.tsv.gz \\
        input/edges2.tsv.gz \\
    | gzip > ProtVerse/edges_alpha${params.disparity_filter_alpha}.tsv.gz
    """

}


/*
run network propagation (random walk with restart), given a list of edges and
a list of seed nodes
[...]
*/
process rwr_lof {

    input:
        path 'input/edges.tsv.gz'
        tuple val(id),
              file('input/seeds.txt')

    output:
        path "rwr_lof/${id}.tsv.gz"

    script:
    """
    mkdir -p rwr_lof

    rwr.py \\
        input/edges.tsv.gz \\
        input/seeds.txt \\
    | gzip > rwr_lof/${id}.tsv.gz
    """

}
