#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { translate } from './modules/utils'
include { translate_no_pub } from './modules/utils'
include { translate_no_pub_rand_id } from './modules/utils'
include { concat } from './modules/utils'
include { concat_w_id } from './modules/utils'
include { translatepy } from './modules/utils'
include { translatepy_no_pub } from './modules/utils'
include { translatepy_no_pub_rand_id } from './modules/utils'
include { translate_expand_matrix } from './modules/utils'
include { translate_expand_pairs } from './modules/utils'
include { translate_expand_list } from './modules/utils'
include { publish } from './modules/utils'
include { split } from './modules/utils'
include { concatenate } from './modules/utils'

include { dl_pc12_reactome } from './modules/pathwaycommons'
include { filter_sig_interactions } from './modules/pathwaycommons'
include { get_reactome_genes } from './modules/pathwaycommons'
include { get_reactome_sig_edges } from './modules/pathwaycommons'

include { run_clusterone } from './modules/clusterone'
//include { cluster_memberships } from './modules/clusterone'
include { cl1_modules_stats } from './modules/clusterone'
include { cl1_joined_modules_stats } from './modules/clusterone'
include { join_clusters } from './modules/clusterone'
include { merge_clusters } from './modules/clusterone'

include { run_blant } from './modules/blant'



// ===== //

workflow REACTOME {

    main:
        // get reactome signalling interactions for human genes
        reactome_sif = dl_pc12_reactome()
        sif = filter_sig_interactions( reactome_sif ).sif
        genes = get_reactome_genes( sif )
        net = get_reactome_sig_edges( sif )
        
    emit:
        sif
        genes
        net

}


workflow CLUSTERONE_WCSN {

    take:
        net

    main:
        id = Channel.from("wcsn_${params.wholecellnet_edge_min_threshold}")
        cl1_min_density = Channel.from(params.cl1_protverse_min_density_list)
        run_clusterone_input = id.combine(net)
                                 .combine(cl1_min_density)
        clusters = run_clusterone( run_clusterone_input )
        join_custers_input = clusters.map{ it -> it[1] }.collect()
        joined_clusters = join_clusters( id, join_custers_input )
        cl1_modules_stats( clusters )
        cl1_joined_modules_stats( joined_clusters )

    emit:
        joined_clusters

}


workflow CLUSTERONE_WCMN {

    take:
        net

    main:
        id = Channel.from("wcmn_${params.wholecellnet_edge_min_threshold}")
        cl1_min_density = Channel.from(params.cl1_protverse_min_density_list)
        run_clusterone_input = id.combine(net)
                                 .combine(cl1_min_density)
        clusters = run_clusterone( run_clusterone_input )
        join_custers_input = clusters.map{ it -> it[1] }.collect()
        joined_clusters = join_clusters( id, join_custers_input )
        cl1_modules_stats( clusters )
        cl1_joined_modules_stats( joined_clusters )

    emit:
        joined_clusters

}


workflow CLUSTERONE_WCSN_REACTOME {

    take:
        net

    main:
        id = Channel.from("wcsn_reactome_${params.wholecellnet_edge_min_threshold}")
        cl1_min_density = Channel.from(params.cl1_protverse_min_density_list)
        run_clusterone_input = id.combine(net)
                                 .combine(cl1_min_density)
        clusters = run_clusterone( run_clusterone_input )
        join_custers_input = clusters.map{ it -> it[1] }.collect()
        joined_clusters = join_clusters( id, join_custers_input )
        cl1_modules_stats( clusters )
        cl1_joined_modules_stats( joined_clusters )

    emit:
        joined_clusters

}


workflow BLANT_WCSN {

    take:
        net

    main:
        id = Channel.from("wcsn_${params.wholecellnet_edge_min_threshold}")
        blant_minW_list = Channel.from(params.blant_minW_list)
        run_blant_input = id.combine(net)
                            .combine(blant_minW_list)
        clusters = run_blant( run_blant_input )

    emit:
        clusters

}


workflow WCSN_WCMN_CLUSTERS {

    take:
        wcsn_clusters
        wcmn_clusters

    main:
        merged_clusters = merge_clusters( wcsn_clusters, wcmn_clusters )
        cl1_joined_modules_stats( merged_clusters )

    emit:
        merged_clusters

}


workflow PUBLISH_CONFIG {

    main:
        config_ch = Channel.fromPath("${projectDir}/extract_modules.config")
        val_ch = Channel.of('workflow/extract_modules.config')
        
        publish( config_ch, val_ch )

}


workflow {

    // Signalling model graph
    wcsn = Channel.fromPath("${wcsn_dir}/edges_${params.wholecellnet_edge_min_threshold}minScore_transformed.tsv")


    // Signalling model graph (Reactome genes only)
    //wcsn_reactome = Channel.fromPath("${wcsn_dir}/edges_reactomeGenes_${params.wholecellnet_edge_min_threshold}minScore.tsv")


    // Metabolism model graph
    wcmn = Channel.fromPath("${wcmn_dir}/edges_${params.wholecellnet_edge_min_threshold}minScore_transformed.tsv")
    //wcmn = Channel.fromPath("${wcmn_dir}/edges_${params.wholecellnet_edge_min_threshold}minScore.tsv")


    // Reactome signalling-specific graph from PathwayCommons
    reactome = REACTOME()


    // Signalling modules
    wcsn_clusters = CLUSTERONE_WCSN( wcsn )
    //blant_clusters = BLANT_WCSN( wcsn )


    // Signalling modules (Reactome genes only)
    //wcsn_reactome_clusters = CLUSTERONE_WCSN_REACTOME( wcsn_reactome )


    // metabolic modules
    wcmn_clusters = CLUSTERONE_WCMN( wcmn )


    // join together signalling model- and metabolic model- based clusters
    wcsn_wcmn_clusters = WCSN_WCMN_CLUSTERS( wcsn_clusters, wcmn_clusters )


    // save nextflow.config
    PUBLISH_CONFIG()

}
