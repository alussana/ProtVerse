#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Use the singularity image to run ClusterONE

.META:
1   module name (int starting from 0)
2   tab-separated list of genes
*/
process run_clusterone {

    container "${projectDir}/envs/cl1.sif"
    cpus "${params.clusterone_n_cores}"
    memory '20G'

    publishDir "${out_dir}", pattern: "modules/${net}/clusterone/*", mode: 'copy'

    input:
        tuple val(net),
              file('input/edges.tsv'),
              val(min_density)

    output:
        tuple val("${net}/clusterone"),
              file("modules/${net}/clusterone/clusters_minD${min_density}.tsv"),
              val(min_density)

    script:
    """
    mkdir -p modules/${net}/clusterone

    java -jar /cluster_one-1.2.jar \
        --min-density "${min_density}" \
        --max-overlap "${params.cl1_max_overlap_match_coeff}" \
        --min-size "${params.cl1_min_size}" \
        input/edges.tsv \
        | awk '{print (NR-1)"\\t"\$0}' \
        > modules/${net}/clusterone/clusters_minD${min_density}.tsv
    """

}


process cl1_modules_stats {

    publishDir "${out_dir}",
               pattern: "modules/${net}/*.pdf",
               mode: 'copy'

    input:
        tuple val(net),
              file('input/modules.tsv'),
              val(min_density)

    output:
        path "modules/${net}/*.pdf"

    script:
    """
    mkdir -p modules/${net}/clusterone/

    cat input/modules.tsv | cut -f2- > modules.tsv

    modules_stats.py \
        modules.tsv \
        modules/${net}/logSizeDistrib_minD${min_density}.pdf \
        modules/${net}/overlapRatioDistrib_minD${min_density}.pdf \
        modules/${net}/pleiotropyDistrib_minD${min_density}.pdf \
        modules/${net}/logPleiotropyDistrib_minD${min_density}.pdf \
        modules/${net}/sizeDistrib_minD${min_density}.pdf \
        modules/${net}/jointSizeOverlapRatio_minD${min_density}.pdf \
        modules/${net}/jointLogSizeOverlapRatio_minD${min_density}.pdf \
        ${net} \
        ${params.cl1_min_size} \
        ${params.max_module_size}
    """

}


process cl1_joined_modules_stats {

    publishDir "${out_dir}",
               pattern: "modules/${net}/*.pdf",
               mode: 'copy'

    input:
        tuple val(net),
              file('input/modules.tsv')

    output:
        path "modules/${net}/*.pdf"

    script:
    """
    mkdir -p modules/${net}/clusterone/

    cat input/modules.tsv | cut -f2- > modules.tsv

    modules_stats.py \
        modules.tsv \
        modules/${net}/logSizeDistrib_joined.pdf \
        modules/${net}/overlapRatioDistrib_joined.pdf \
        modules/${net}/pleiotropyDistrib_joined.pdf \
        modules/${net}/logPleiotropyDistrib.pdf \
        modules/${net}/sizeDistrib_joined.pdf \
        modules/${net}/jointSizeOverlapRatio_joined.pdf \
        modules/${net}/jointLogSizeOverlapRatio_joined.pdf \
        ${net} \
        ${params.cl1_min_size} \
        ${params.max_module_size}
    """

}


/*
join different modules sets obtained with different min-density parameters
into a single modules set

modules are renamed accordingly to the min-density value used to detect them

filter out modules where the number of member genes is greater than 
${params.max_module_size}
*/
process join_clusters {

    publishDir "${out_dir}",
               pattern: "modules/${net}/clusters_joined.tsv",
               mode: 'copy'

    input:
        val net
        path "input/*"

    output:
        tuple val("${net}_joined"),
              file("modules/${net}/clusters_joined.tsv")

    script:
    """
    mkdir -p modules/${net}/

    for file in \$(ls input/); do \
        d=\$(echo "\$file" | sed -r 's/clusters_minD(.\\...)\\.tsv/\\1/'); \
        cat input/\$file \
        | awk -v d=\${d} '{print "c"\$1"_d"d"\\t"\$0}' \
        | cut -f1,3-; \
    done \
    > modules/${net}/clusters_unfiltered.tsv
    
    
    cat modules/${net}/clusters_unfiltered.tsv \
        | awk -F "\\t" 'NF<=${params.max_module_size}+1' \
        > modules/${net}/clusters_joined.tsv
    """

}

/*
== TODO ==

Build vector of cluster membership for each gene from the ClusterONE output
*/
process CLUSTER_MEMBERSHIPS {

    publishDir "${out_dir}", pattern: "modules/${net}/clusterone/*", mode: 'copy'

    input:
        tuple val(net),
              file('input/clusters.tsv')

    output:
        tuple val(net),
              file("modules/${net}/clusterone/vectors.tsv")

    script:
    """
    mkdir -p modules/${net}

    cat input/clusters.tsv | tr '\t' '\n' | sort | uniq \
        > nodes.txt

    clusters2vec.py \
        nodes.txt \
        input/clusters.tsv \
        > modules/${net}/clusterone/vectors.tsv
    """

}


/*
join different modules sets obtained with different min-density parameters
into a single modules set

modules are renamed accordingly to the min-density value used to detect them

filter out modules where the number of member genes is greater than 
${params.max_module_size}
*/
process merge_clusters {

    publishDir "${out_dir}",
               pattern: "modules/${id1}__${id2}/clusters_merged.tsv",
               mode: 'copy'

    input:
        tuple val(id1),
              file("input/clusters1.tsv")
        tuple val(id2),
              file("input/clusters2.tsv")

    output:
        tuple val("${id1}__${id2}_merged"),
              file("modules/${id1}__${id2}/clusters_merged.tsv")

    script:
    """
    mkdir -p modules/${id1}__${id2}/
  
    cat \\
        <(awk 'BEGIN{OFS="\t"} {printf "${id1}__"\$1"\t"; for (i=2; i<=NF; i++) printf "%s%s", \$i, (i<NF ? OFS : ORS)}' input/clusters1.tsv) \\
        <(awk 'BEGIN{OFS="\t"} {printf "${id2}__"\$1"\t"; for (i=2; i<=NF; i++) printf "%s%s", \$i, (i<NF ? OFS : ORS)}' input/clusters2.tsv) \\
        > modules/${id1}__${id2}/clusters_merged.tsv
    """

}