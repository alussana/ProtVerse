#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { translate } from './modules/utils'
include { translate_no_pub } from './modules/utils'
include { translate_no_pub_rand_id } from './modules/utils'
include { filter_data_matrix_rows } from './modules/utils'
include { filter_data_matrix_rows_w_id } from './modules/utils'
include { concat } from './modules/utils'
include { cat_and_split } from './modules/utils'
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

include { dl_pc13_intact } from './modules/pathwaycommons'
include { unweighted_intact } from './modules/pathwaycommons'
include { dl_pc13_biogrid } from './modules/pathwaycommons'
include { unweighted_biogrid } from './modules/pathwaycommons'

include { dl_reactome } from './modules/reactome'
include { sig_react_ids } from './modules/reactome'
include { immune_system_react_ids } from './modules/reactome'
include { ensembl_react_lists } from './modules/reactome'
include { pos_examples } from './modules/reactome'
include { get_signalling_genes } from './modules/reactome'
include { list_terms } from './modules/reactome'
include { translate_tab } from './modules/reactome'
include { hs_react_sets } from './modules/reactome'
include { translate_tab_all } from './modules/reactome'
include { hs_gene_list } from './modules/reactome'
include { pos_ex_unique_edges } from './modules/reactome'
include { dl_reactome2022_enrichr } from './modules/reactome'
include { modules_stats } from './modules/reactome'
include { signal_transduction_reactome_leaves } from './modules/reactome'

include { metabolism_react_ids } from './modules/reactome'
include { metabolism_ensembl_react_lists } from './modules/reactome'
include { get_metabolism_genes } from './modules/reactome'
include { metabolism_list_terms } from './modules/reactome'

include { download_median_tpm } from './modules/gtex'
include { parse_median_tpm } from './modules/gtex'
include { get_gtex_genes } from './modules/gtex'
include { gtex_featvec } from './modules/gtex'

include { download_eprot } from './modules/eprot'
include { parse_eprot } from './modules/eprot'
include { merge_eprot_exp } from './modules/eprot'
include { get_eprot_genes } from './modules/eprot'
include { eprot_featvec } from './modules/eprot'

include { dl_proteomehd } from './modules/proteomehd'
include { parse_proteomehd } from './modules/proteomehd'
include { get_proteomehd_genes } from './modules/proteomehd'
include { proteomehd_featvec  } from './modules/proteomehd'

include { dl_mitchell2023 } from './modules/mitchell2023'
include { get_mitchell2023_genes } from './modules/mitchell2023'
include { parse_mitchell2023 } from './modules/mitchell2023'
include { mitchell2023_featvec } from './modules/mitchell2023'

include { open_ptmdb } from './modules/ptmdb'
include { parse_ptmdb } from './modules/ptmdb'
include { merge_ptmdb_exp } from './modules/ptmdb'
include { get_ptmdb_genes } from './modules/ptmdb'
include { ptmdb_featvec } from './modules/ptmdb'

include { dl_lopit2025 } from './modules/lopit2025'
include { lopit2025_prot_loc_changes } from './modules/lopit2025'
include { get_lopit2025_genes } from './modules/lopit2025'
include { lopit2025_featvec } from './modules/lopit2025'

/*include { dl_tahoe_100m } from './modules/tahoe'
include { dl_tahoe_100m_metadata } from './modules/tahoe'
include { dl_tahoe_100m_hf } from './modules/tahoe'*/

include { parse_ptmdb_imputed } from './modules/ptmdb_imputed'
include { get_ptmdb_imputed_genes } from './modules/ptmdb_imputed'
include { ptmdb_imputed_featvec } from './modules/ptmdb_imputed'

include { download_essentiality_matrices } from './modules/depmap'
include { get_depmap_genes } from './modules/depmap'
include { dependency_featvec } from './modules/depmap'

include { dl_depmap20q4v2 } from './modules/depmap20q4v2'
include { parse_depmap20q4v2 } from './modules/depmap20q4v2'
include { covar_norm } from './modules/depmap20q4v2'
include { get_depmap20q4v2_genes } from './modules/depmap20q4v2'

include { download_orthogroups } from './modules/orthogroups2015'
include { get_orthogroups_genes } from './modules/orthogroups2015'
include { orthogroup_featvec } from './modules/orthogroups2015'

include { download_humap3_net } from './modules/humap3'
include { humap3_table } from './modules/humap3'
include { get_humap3_genes } from './modules/humap3'
include { humap3_score } from './modules/humap3'

include { get_ubiquitination } from './modules/ubiquitination'
include { parse_ubiquitination } from './modules/ubiquitination'
include { get_ubiquitination_genes } from './modules/ubiquitination'
include { ubiquitination_featvec } from './modules/ubiquitination'

include { dl_ref_proteome } from './modules/uniprot'
include { topo_domain } from './modules/uniprot'
include { interpro } from './modules/uniprot'
include { length } from './modules/uniprot'
include { dl_human_ids } from './modules/uniprot'
include { filter_idmapping } from './modules/uniprot'
include { get_proteome_genes } from './modules/uniprot'
include { translate_ids; translate_ids as tr_ids } from './modules/uniprot'
include { filter_id_dict } from './modules/uniprot'
include { uniprot2gene_name_and_synonym } from './modules/uniprot'
include { translate_ac_to_gene_name } from './modules/uniprot'
include { viz_ref_genes_coverage } from './modules/uniprot'

include { dl_bernett_2024 } from './modules/bernett2024'
include { parse_and_translate_bernett_2024 } from './modules/bernett2024'
include { bernett2024_conform_cat_and_split } from './modules/bernett2024'
include { bernett2024_features_tables } from './modules/bernett2024'
include { bernett2024_features_pca } from './modules/bernett2024'
include { bernett2024_make_data_splits } from './modules/bernett2024'
include { bernett2024_features_distrib } from './modules/bernett2024'
include { bernett2024_train_valid_rf } from './modules/bernett2024'
include { bernett2024_train_valid_rf_reduced } from './modules/bernett2024'
include { bernett2024_train_valid_xgb } from './modules/bernett2024'
include { bernett2024_xgb_single_source } from './modules/bernett2024'

include { neg_examples } from './modules/training'
include { features_tables } from './modules/training'
include { features_tables_metabolism } from './modules/training'
include { features_pca } from './modules/training'
include { features_pca_metabolism } from './modules/training'
include { split_dataset } from './modules/training'
include { split_dataset_metabolism } from './modules/training'
include { inspect_examples } from './modules/training'
include { shap } from './modules/training'
include { features_distrib } from './modules/training'
include { features_distrib_metabolism } from './modules/training'
include { calibrate } from './modules/training'
include { signalling_train_valid_xgb } from './modules/training'
include { signalling_xgb_single_source } from './modules/training'
include { metabolism_train_valid_xgb } from './modules/training'
include { metabolism_xgb_single_source } from './modules/training'

include { add_examples } from './modules/rf_test'

include { gene_pairs } from './modules/build_network'
include { eprot_features } from './modules/build_network'
include { proteomehd_features } from './modules/build_network'
include { mitchell2023_features } from './modules/build_network'
include { gtex_features } from './modules/build_network'
include { ptmdb_features } from './modules/build_network'
include { ubiquitination_features } from './modules/build_network'
include { dep_features } from './modules/build_network'
include { orthogroup_features } from './modules/build_network'
include { humap_features } from './modules/build_network'
include { lopit2025_features } from './modules/build_network'
include { compute_edge_features } from './modules/build_network'
include { compute_edge_weights; compute_edge_weights as metabolism_compute_edge_weights } from './modules/build_network'
include { concat_edges; metabolism_concat_edges } from './modules/build_network'
include { edges_tsv_2_parquet } from './modules/build_network'
include { add_edges } from './modules/build_network'
include { filter_edge_score; metabolism_filter_edge_score } from './modules/build_network'
include { transform_edge_score; metabolism_transform_edge_score } from './modules/build_network'
include { get_ref_genes } from './modules/build_network'
include { gene_adjacency_vector } from './modules/build_network'
include { correlate_adjacencies } from './modules/build_network'
include { w_net_stats } from './modules/build_network'
include { filter_reactome_genes } from './modules/build_network'


workflow REACTOME {

    main:
        // Reactome signal transduction data
        reactome = dl_reactome()

        sig_transd_reactome = sig_react_ids( reactome )
        
        imm_reactome = immune_system_react_ids( reactome )
        sig_reactome = concatenate( sig_transd_reactome.concat(imm_reactome).collect() )
        //sig_reactome = concatenate( sig_transd_reactome.collect() )

        cofun_uniprot_ids = ensembl_react_lists( sig_reactome,
                                                reactome.uniprot_ids )
        cofun_uniprot_ids_collected = cofun_uniprot_ids.collect()
        gene_list = get_signalling_genes( cofun_uniprot_ids_collected )
        mapped_cofun_ids = cofun_uniprot_ids
                                .flatMap()
                                .map{file -> tuple( file.baseName, file )}
        clusters = list_terms( cofun_uniprot_ids_collected )
        hs_sets = hs_react_sets(reactome.pathway_ids, reactome.uniprot_ids)
        hs_gene_list = hs_gene_list( hs_sets )

        // Reactome metabolism data
        metabolism_reactome = metabolism_react_ids( reactome )
        metabolism_cofun_uniprot_ids = metabolism_ensembl_react_lists( metabolism_reactome,
                                                                       reactome.uniprot_ids )
        metabolism_cofun_uniprot_ids_collected = metabolism_cofun_uniprot_ids.collect()
        metabolism_gene_list = get_metabolism_genes( metabolism_cofun_uniprot_ids_collected )
        metabolism_mapped_cofun_ids = metabolism_cofun_uniprot_ids
                                            .flatMap()
                                            .map{file -> tuple( file.baseName, file )}
        metabolism_clusters = metabolism_list_terms( metabolism_cofun_uniprot_ids_collected )

    emit:
        mapped_cofun_ids
        gene_list
        clusters
        metabolism_mapped_cofun_ids
        metabolism_gene_list
        metabolism_clusters
        hs_sets
        hs_gene_list
                
}


workflow PROTEOMEHD {

    main: 
        table_S1 = dl_proteomehd()
        matrix = parse_proteomehd( table_S1 )
        gene_list = get_proteomehd_genes( matrix )

    emit:
        matrix
        gene_list

}


workflow MITCHELL2023 {

    main: 
        table_S1 = dl_mitchell2023()
        matrix = parse_mitchell2023( table_S1 )
        gene_list = get_mitchell2023_genes( matrix )

    emit:
        matrix
        gene_list

}


workflow GTEX {
    
    main:
        // gtex transcript expression data
        median_tpm = download_median_tpm()
        matrix = parse_median_tpm( median_tpm )
        gene_list = get_gtex_genes( matrix )

    emit:
        matrix
        gene_list

}


workflow EPROT {

    main:
        // Expression Atlas protein abundance data
        exp_ids = channel.fromList( params.eprot_ids )
        eprot_tables = download_eprot( exp_ids )
                            .flatMap()
                            .map{ file -> tuple( file.baseName, file ) }
        eprot_parsed = parse_eprot( eprot_tables ) 
                            .collect()
        matrix = merge_eprot_exp( eprot_parsed )
        gene_list = get_eprot_genes( matrix )

    emit:
        matrix
        gene_list

}


workflow PTMDB {

    main:
        // ptmdb phosphoproteomics data
        zip = channel.fromPath( ptmdb_text_format )
        phos_exp = open_ptmdb( zip )
                        .flatMap()
                        .map{ file -> tuple( file.baseName, file ) }
        phos_logfc = parse_ptmdb( phos_exp )
                        .collect()
        matrix = merge_ptmdb_exp( phos_logfc )
        gene_list = get_ptmdb_genes( matrix )

    emit:
        matrix
        gene_list

}


/*workflow TAHOE {

    main:
        // download from google cloud
        //h5ad_paths = dl_tahoe_100m()
        //metadata = dl_tahoe_100m_metadata()

        // download from hugging face
        dl_tahoe_100m_hf()

        
}*/


workflow LOPIT2025 {

    main:
        // download lopit2025 dataset
        //lopit2025 = dl_lopit2025()

        // build protein localization changes dataset
        lopit2025_xlsx = Channel.fromPath( "${lopit2025_xlsx_path}" )
        table = lopit2025_prot_loc_changes( lopit2025_xlsx )

        // get list of entries (accessions)
        gene_list = get_lopit2025_genes( table )


    emit:
        table
        gene_list

}


workflow PTMDB_IMPUTED {

    main:
        // ptmdb phosphoproteomics data
        h5_file = Channel.fromPath( "${ptmdb_imputed_h5}" )
        matrix = parse_ptmdb_imputed( h5_file )
        gene_list = get_ptmdb_imputed_genes( matrix )

    emit:
        matrix
        gene_list

}


workflow DEPENDENCY {

    main:
        // Integrated Broad and Sanger CRISPR-KO fold changes
        dependency = download_essentiality_matrices()
        gene_list = get_depmap_genes( dependency.CRISPRcleanR_FC )
        matrix = dependency.CRISPRcleanR_FC
    emit:
        matrix
        gene_list

}


workflow DEPMAP20Q2V2 {

    main:
        ceres_scores = dl_DepMap20Q4v2()
        ceres_scores = parse_DepMap20Q4v2( ceres_scores )
        gene_list = get_DepMap20Q4v2_genes( ceres_scores )
        matrix = covar_norm(ceres_scores).X_white

    emit:
        matrix
        gene_list

}


workflow ORTHOGROUPS {

    main:
        // Phylogenetic Co-occurrence Score and orthogroups
        orthogroups = download_orthogroups()
        pcs = orthogroups.PCS
        occurrence_matrix = orthogroups.occurrence_matrix
        gene2orthogroup = orthogroups.gene2orthogroup
        gene_list = get_orthogroups_genes( orthogroups.gene2orthogroup )

    emit:
        pcs
        gene2orthogroup
        occurrence_matrix
        gene_list
}


workflow HUMAP3 {

    main:
        // Get hu.MAP 3.0 data
        net = download_humap3_net()
        table = humap3_table( net )
        gene_list = get_humap3_genes( net )
    
    emit:
        table
        gene_list

}


workflow UNIPROT {

    main:
        ref_prot = dl_ref_proteome()
        idmapping_all = dl_human_ids()
        idmapping = filter_idmapping( idmapping_all )
        uniprot2gene_name_and_synonym_dict = uniprot2gene_name_and_synonym( idmapping_all )
        gene_list = get_proteome_genes( ref_prot )

    emit:
        idmapping
        gene_list
        uniprot2gene_name_and_synonym_dict
}


workflow UBIQUITINATION {

    // Get ubiquitination data from the source file and parse it

    main:
        source = get_ubiquitination( ubiquitination_data_source )
        table = parse_ubiquitination( source )
        gene_list = get_ubiquitination_genes( table )

    emit:
        table
        gene_list    

}


workflow MAP_IDS {

    // generate dictionary table assigning to each reference proteome gene name
    // a corresponding translation in each of the datasets

    take:
        ref_genes
        idmapping
        signalling_genes
        metabolism_genes
        gtex_genes
        eprot_genes
        orthogroups_genes
        depmap_genes
        ubiquitination_genes
        humap3_genes
        ptmdb_genes
        proteomehd_genes
        mitchell2023_genes
        lopit2025_genes
        all_reactome_genes

    main:
        // only use reactome signal trasduction leaves genes
        unfiltered = translate_ids( ref_genes,
                                    idmapping,
                                    signalling_genes,
                                    metabolism_genes,
                                    gtex_genes,
                                    eprot_genes,
                                    orthogroups_genes,
                                    depmap_genes,
                                    ubiquitination_genes,
                                    humap3_genes,
                                    ptmdb_genes,
                                    proteomehd_genes,
                                    mitchell2023_genes,
                                    lopit2025_genes,
                                    all_reactome_genes )

        viz_ref_genes_coverage( unfiltered )

        // only keep the genes found in all the data sources
        //filtered = filter_id_dict( unfiltered )
    
    emit:
        unfiltered
        //filtered

}


workflow FILTER_GTEX {
    
    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows(matrix, dict, dict_col)
        id = Channel.from('gtex')
        tr_input = id.combine(target_matrix)
        f_matrix = translate_expand_matrix( tr_input, dict, dict_col )

    emit:
        f_matrix

}


workflow FILTER_EPROT {
    
    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows(matrix, dict, dict_col)
        id = Channel.from('eprot')
        tr_input = id.combine(target_matrix)
        f_matrix = translate_expand_matrix(tr_input, dict, dict_col )

    emit:
        f_matrix

}


workflow FILTER_PROTEOMEHD {
    
    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows(matrix, dict, dict_col)
        id = Channel.from('proteomehd')
        tr_input = id.combine(target_matrix)
        f_matrix = translate_expand_matrix(tr_input, dict, dict_col )

    emit:
        f_matrix

}


workflow FILTER_MITCHELL2023 {

    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows(matrix, dict, dict_col)
        id = Channel.from('mitchell2023')
        tr_input = id.combine(target_matrix)
        f_matrix = translate_expand_matrix(tr_input, dict, dict_col )

    emit:
        f_matrix

}


workflow FILTER_LOPIT2025 {

    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows(matrix, dict, dict_col)
        id = Channel.from('lopit2025')
        tr_input = id.combine(target_matrix)
        f_matrix = translatepy(tr_input, dict, dict_col )

    emit:
        f_matrix

}


workflow FILTER_REACTOME {
    
    take:
        gene_lists__dict__col
        dict
        dict_col

    main:
        target_gene_lists = filter_data_matrix_rows_w_id( gene_lists__dict__col )
        target_lists__dict__col = target_gene_lists
                                        .combine(dict)
                                        .combine(dict_col)
        tr_gene_lists = translate_expand_list( target_lists__dict__col )
                            .flatMap()
                            .map{file -> tuple( file.baseName, file )}

    emit:
        tr_gene_lists

}


workflow FILTER_METABOLISM {

    take:
        gene_lists__dict__col
        dict
        dict_col

    main:
        target_gene_lists = filter_data_matrix_rows_w_id( gene_lists__dict__col )
        target_lists__dict__col = target_gene_lists
                                        .combine(dict)
                                        .combine(dict_col)
        tr_gene_lists = translate_expand_list( target_lists__dict__col )
                            .flatMap()
                            .map{file -> tuple( file.baseName, file )}

    emit:
        tr_gene_lists

}


workflow FILTER_DEPENDENCY {
    
    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows(matrix, dict, dict_col)
        id = Channel.from('dependency')
        tr_input = id.combine(target_matrix)
        f_matrix = translate_expand_matrix(tr_input, dict, dict_col )

    emit:
        f_matrix

}


workflow FILTER_PTMDB {
    
    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows(matrix, dict, dict_col)
        id = Channel.from('ptmdb')
        tr_input = id.combine(target_matrix)
        f_matrix = translate_expand_matrix(tr_input, dict, dict_col )

    emit:
        f_matrix

}


workflow FILTER_ORTHOGROUPS {
    
    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows( matrix, dict, dict_col )
        id = Channel.from( 'orthogroups' )
        tr_input = id.combine( target_matrix )
        f_matrix = translate_expand_matrix( tr_input, dict, dict_col )

    emit:
        f_matrix

}


workflow FILTER_UBIQUITINATION {
    
    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows(matrix, dict, dict_col)
        id = Channel.from('ubiquitination')
        tr_input = id.combine(target_matrix)
        f_matrix = translate_expand_matrix(tr_input, dict, dict_col )

    emit:
        f_matrix

}


workflow FILTER_IVKAPHE {

    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows(matrix, dict, dict_col)
        id = Channel.from('ivkaphe')
        target_matrix_chunks = target_matrix
                                    .splitText( by: 10000 )
        tr_input_chunks = id.combine( target_matrix_chunks )
                            .combine( dict )
                            .combine( dict_col )
        f_matrix_chunks = translate_expand_pairs( tr_input_chunks )
                            .collect()
        f_matrix_name = Channel.from('ivkaphe.tsv')
        header = Channel.from('kinase\ttarget\tpos\tivkaphe_score')
        f_matrix = CONCAT_W_ID( f_matrix_name, f_matrix_chunks, header )
        kinases = get_ivkaphe_kinases( f_matrix )

    emit:
        f_matrix
        kinases

}


workflow FILTER_HUMAP3 {
    
    take:
        matrix
        dict
        dict_col

    main:
        target_matrix = filter_data_matrix_rows(matrix, dict, dict_col)
        id = Channel.from('humap3')
        target_matrix_chunks = target_matrix
                                    .splitText( by: 50000 )
        tr_input_chunks = id.combine( target_matrix_chunks )
                            .combine( dict )
                            .combine( dict_col )
        f_matrix_chunks = translate_expand_pairs( tr_input_chunks )
                            .collect()
        f_matrix_name = Channel.from('humap3.tsv')
        header = Channel.from('gene_1\tgene_2\thumap_score')
        f_matrix = concat_w_id( f_matrix_name, f_matrix_chunks, header )

    emit:
        f_matrix

}


workflow BERNETT2024 {

    take:
        uniprot2gene_name_and_synonym_dict

    main:
        bernett_2024_zip = dl_bernett_2024()

        bernett_2024 = parse_and_translate_bernett_2024( bernett_2024_zip,
                                                         uniprot2gene_name_and_synonym_dict )

        train_ids = bernett_2024.train_ids
        valid_ids = bernett_2024.valid_ids
        test_ids = bernett_2024.test_ids

    emit:
        train_ids
        valid_ids
        test_ids
}


// ===== //

workflow GENERATE_EXAMPLES {

    take:
        f_reactome_lists
        id_dict
        f_gtex_m
        f_eprot_m
        f_proteomehd_m
        f_mitchell2023_m
        f_dependency_m
        f_orthogroup_m
        orthogroup_pcs
        f_ptmdb_m
        f_ubiquitination_m
        f_humap3_m
        f_lopit2025_m

    main:
        // generate training data
        pos_examples = pos_examples( f_reactome_lists ).collect()
        pos_examples_uniq = pos_ex_unique_edges( pos_examples )
        n = pos_examples_uniq.n
        pos_examples = pos_examples_uniq.edges
        neg_examples = neg_examples( id_dict, n )
        concat_input = pos_examples
                            .concat( neg_examples )
                            .collect()
        
        // split gene pairs training examples in chunks
        examples = cat_and_split( concat_input, 10000 )
                        .flatten()

        // build features from eprot
        eprot_featvec_input = examples.combine( f_eprot_m )
        eprot_featvec_ch = eprot_featvec( eprot_featvec_input ).collect()
        
        // build features from proteomehd
        proteomehd_featvec_input = examples.combine( f_proteomehd_m )
        proteomehd_featvec_ch = proteomehd_featvec( proteomehd_featvec_input )
                                    .collect()

        // build features from mitechell2023
        mitchell2023_featvec_input = examples.combine( f_mitchell2023_m )
        mitchell2023_featvec_ch = mitchell2023_featvec( mitchell2023_featvec_input )
                                    .collect()
        
        // build features from gtex
        gtex_featvec_input = examples.combine( f_gtex_m )
        gtex_featvec_ch = gtex_featvec( gtex_featvec_input )
                            .collect()

        // build features from ptmdb
        ptmdb_featvec_input = examples.combine( f_ptmdb_m )
        ptmdb_featvec_ch = ptmdb_featvec( ptmdb_featvec_input )
                                .collect()

        // build features from ubiquitination
        ubiquitination_featvec_input = examples.combine( f_ubiquitination_m )
        ubiquitination_featvec_ch = ubiquitination_featvec( ubiquitination_featvec_input )
                                        .collect()

        // build features from depmap
        dependency_featvec_input = examples.combine( f_dependency_m )
        dependency_featvec_ch = dependency_featvec( dependency_featvec_input )
                                    .collect()

        // build features from orthology information
        orthogroup_featvec_input = examples
                                    .combine( f_orthogroup_m )
                                    .combine( orthogroup_pcs )
        orthogroup_featvec_ch = orthogroup_featvec( orthogroup_featvec_input )
                                    .collect()

        // build features from humap3
        humap3_score_input = examples.combine( f_humap3_m )
        humap3_score_ch = humap3_score( humap3_score_input )
                            .collect()
        // build features from lopit`2025
        lopit2025_featvec_input = examples.combine( f_lopit2025_m )
        lopit2025_featvec_ch = lopit2025_featvec( lopit2025_featvec_input )
                                    .collect()
        // join features together
        dataset = features_tables( eprot_featvec_ch,
                                   proteomehd_featvec_ch,
                                   mitchell2023_featvec_ch,
                                   gtex_featvec_ch,
                                   ptmdb_featvec_ch,
                                   dependency_featvec_ch,
                                   ubiquitination_featvec_ch,
                                   orthogroup_featvec_ch,
                                   humap3_score_ch,
                                   lopit2025_featvec_ch )/*,
                                   ivkaphe_featvec_ch*/

        features_pca( dataset )
        
        // get balanced random splits for training and testing
        splits = split_dataset( dataset )
        //inspect_examples( splits.X_train,
        //                  splits.X_test,
        //                  splits.y_train,
        //                  splits.y_test )

        features_distrib( dataset )

        X_train = splits.X_train
        X_test = splits.X_test
        y_train = splits.y_train
        y_test = splits.y_test
        pos_examples = pos_examples.collect()
    
    emit:
        X_train
        X_test
        y_train
        y_test
        pos_examples
        n
}


workflow GENERATE_EXAMPLES_METABOLISM {

    take:
        f_reactome_lists
        id_dict
        f_gtex_m
        f_eprot_m
        f_proteomehd_m
        f_mitchell2023_m
        f_dependency_m
        f_orthogroup_m
        orthogroup_pcs
        f_ptmdb_m
        f_ubiquitination_m
        f_humap3_m
        f_lopit2025_m
        

    main:
        // generate training data
        pos_examples = pos_examples( f_reactome_lists ).collect()
        pos_examples_uniq = pos_ex_unique_edges( pos_examples )
        n = pos_examples_uniq.n
        pos_examples = pos_examples_uniq.edges
        neg_examples = neg_examples( id_dict, n )
        concat_input = pos_examples
                            .concat( neg_examples )
                            .collect()
        
        // split gene pairs training examples in chunks
        examples = cat_and_split( concat_input, 10000 )
                        .flatten()

        // build features from eprot
        eprot_featvec_input = examples.combine( f_eprot_m )
        eprot_featvec_ch = eprot_featvec( eprot_featvec_input ).collect()
        
        // build features from proteomehd
        proteomehd_featvec_input = examples.combine( f_proteomehd_m )
        proteomehd_featvec_ch = proteomehd_featvec( proteomehd_featvec_input )
                                    .collect()

        // build features from mitechell2023
        mitchell2023_featvec_input = examples.combine( f_mitchell2023_m )
        mitchell2023_featvec_ch = mitchell2023_featvec( mitchell2023_featvec_input )
                                    .collect()
        
        // build features from gtex
        gtex_featvec_input = examples.combine( f_gtex_m )
        gtex_featvec_ch = gtex_featvec( gtex_featvec_input )
                            .collect()

        // build features from ptmdb
        ptmdb_featvec_input = examples.combine( f_ptmdb_m )
        ptmdb_featvec_ch = ptmdb_featvec( ptmdb_featvec_input )
                                .collect()

        // build features from ubiquitination
        ubiquitination_featvec_input = examples.combine( f_ubiquitination_m )
        ubiquitination_featvec_ch = ubiquitination_featvec( ubiquitination_featvec_input )
                                        .collect()

        // build features from depmap
        dependency_featvec_input = examples.combine( f_dependency_m )
        dependency_featvec_ch = dependency_featvec( dependency_featvec_input )
                                    .collect()

        // build features from orthology information
        orthogroup_featvec_input = examples
                                    .combine( f_orthogroup_m )
                                    .combine( orthogroup_pcs )
        orthogroup_featvec_ch = orthogroup_featvec( orthogroup_featvec_input )
                                    .collect()

        // build features from humap3
        humap3_score_input = examples.combine( f_humap3_m )
        humap3_score_ch = humap3_score( humap3_score_input )
                            .collect()
        
        // build features from lopit2025
        lopit2025_featvec_input = examples.combine( f_lopit2025_m )
        lopit2025_featvec_ch = lopit2025_featvec( lopit2025_featvec_input )
                                    .collect()

        // join features together
        dataset = features_tables_metabolism( eprot_featvec_ch,
                                              proteomehd_featvec_ch,
                                              mitchell2023_featvec_ch,
                                              gtex_featvec_ch,
                                              ptmdb_featvec_ch,
                                              dependency_featvec_ch,
                                              ubiquitination_featvec_ch,
                                              orthogroup_featvec_ch,
                                              humap3_score_ch,
                                              lopit2025_featvec_ch )

        features_pca_metabolism( dataset )
        
        // get balanced random splits for training and testing
        splits = split_dataset_metabolism( dataset )
        //inspect_examples( splits.X_train,
        //                  splits.X_test,
        //                  splits.y_train,
        //                  splits.y_test )

        features_distrib_metabolism( dataset )

        X_train = splits.X_train
        X_test = splits.X_test
        y_train = splits.y_train
        y_test = splits.y_test
        pos_examples = pos_examples.collect()
    
    emit:
        X_train
        X_test
        y_train
        y_test
        pos_examples
        n
}


workflow GENERATE_EXAMPLES_BERNETT2024 {

    take:
        train_ids
        valid_ids
        test_ids
        id_dict
        f_gtex_m
        f_eprot_m
        f_proteomehd_m
        f_mitchell2023_m
        f_dependency_m
        f_orthogroup_m
        orthogroup_pcs
        f_ptmdb_m
        f_ubiquitination_m
        f_humap3_m
        f_lopit2025_m

    main:
        // concatenate all gene pairs and split them in chunks
        examples = bernett2024_conform_cat_and_split(
                            train_ids,
                            valid_ids,
                            test_ids,
                            10000 ).flatten()

        // build features from eprot
        eprot_featvec_input = examples.combine( f_eprot_m )
        eprot_featvec_ch = eprot_featvec( eprot_featvec_input ).collect()
        
        // build features from proteomehd
        proteomehd_featvec_input = examples.combine( f_proteomehd_m )
        proteomehd_featvec_ch = proteomehd_featvec( proteomehd_featvec_input )
                                    .collect()

        // build features from mitechell2023
        mitchell2023_featvec_input = examples.combine( f_mitchell2023_m )
        mitchell2023_featvec_ch = mitchell2023_featvec( mitchell2023_featvec_input )
                                    .collect()
        
        // build features from gtex
        gtex_featvec_input = examples.combine( f_gtex_m )
        gtex_featvec_ch = gtex_featvec( gtex_featvec_input )
                            .collect()

        // build features from ptmdb
        ptmdb_featvec_input = examples.combine( f_ptmdb_m )
        ptmdb_featvec_ch = ptmdb_featvec( ptmdb_featvec_input )
                                .collect()

        // build features from ubiquitination
        ubiquitination_featvec_input = examples.combine( f_ubiquitination_m )
        ubiquitination_featvec_ch = ubiquitination_featvec( ubiquitination_featvec_input )
                                        .collect()

        // build features from depmap
        dependency_featvec_input = examples.combine( f_dependency_m )
        dependency_featvec_ch = dependency_featvec( dependency_featvec_input )
                                    .collect()

        // build features from orthology information
        orthogroup_featvec_input = examples
                                    .combine( f_orthogroup_m )
                                    .combine( orthogroup_pcs )
        orthogroup_featvec_ch = orthogroup_featvec( orthogroup_featvec_input )
                                    .collect()

        // build features from humap3
        humap3_score_input = examples.combine( f_humap3_m )
        humap3_score_ch = humap3_score( humap3_score_input )
                            .collect()

        // build features from lopit2025
        lopit2025_featvec_input = examples.combine( f_lopit2025_m )
        lopit2025_featvec_ch = lopit2025_featvec( lopit2025_featvec_input )
                                    .collect()
        
        // join features together
        dataset = bernett2024_features_tables( eprot_featvec_ch,
                                               proteomehd_featvec_ch,
                                               mitchell2023_featvec_ch,
                                               gtex_featvec_ch,
                                               ptmdb_featvec_ch,
                                               ubiquitination_featvec_ch,
                                               dependency_featvec_ch,
                                               orthogroup_featvec_ch,
                                               humap3_score_ch,
                                               lopit2025_featvec_ch )

        bernett2024_features_pca( dataset )


        splits = bernett2024_make_data_splits( dataset, 
                                               train_ids,
                                               valid_ids,
                                               test_ids )

        //inspect_examples( splits.X_train,
        //                  splits.X_test,
        //                  splits.y_train,
        //                  splits.y_test )*/

        bernett2024_features_distrib( dataset )

        X_train = splits.X_train
        X_valid = splits.X_valid
        X_test = splits.X_test
        y_train = splits.y_train
        y_valid = splits.y_valid
        y_test = splits.y_test

    emit:
        X_train
        X_valid
        X_test
        y_train
        y_valid
        y_test
    
}


workflow FIT_SIGNALLING {

    take:
        X_train
        X_valid
        y_train
        y_valid

    main:
        xgb = signalling_train_valid_xgb( X_train, 
                                          X_valid,
                                          y_train,
                                          y_valid ).model
        omics_model = xgb

    emit:
        omics_model

}


workflow FIT_SIGNALLING_SINGLE_SOURCE {

    take:
        X_train
        X_valid
        y_train
        y_valid

    main:
        xgb = signalling_xgb_single_source( X_train, 
                                            X_valid,
                                            y_train,
                                            y_valid )

}


workflow FIT_METABOLISM {

    take:
        X_train
        X_valid
        y_train
        y_valid

    main:
        xgb = metabolism_train_valid_xgb( X_train, 
                                          X_valid,
                                          y_train,
                                          y_valid ).model
        omics_model = xgb

    emit:
        omics_model

}


workflow FIT_METABOLISM_SINGLE_SOURCE {

    take:
        X_train
        X_valid
        y_train
        y_valid

    main:
        xgb = metabolism_xgb_single_source( X_train, 
                                            X_valid,
                                            y_train,
                                            y_valid )

}


workflow FIT_BERNETT2024 {

    take:
        X_train
        X_valid
        X_test
        y_train
        y_valid
        y_test

    main:

        xgb = bernett2024_train_valid_xgb( X_train, 
                                           X_test,
                                           y_train,
                                           y_test ).model
                                                  
}


workflow FIT_BERNETT2024_SINGLE_SOURCE {

    take:
        X_train
        X_valid
        X_test
        y_train
        y_valid
        y_test

    main:

        xgb = bernett2024_xgb_single_source( X_train, 
                                             X_test,
                                             y_train,
                                             y_test )
                                                  
}


workflow ADD_NEG_EXAMPLES {

    take:
        id_dict
        f_gtex_m
        f_eprot_m
        f_proteomehd_m
        f_mitchell2023_m
        f_dependency_m
        f_orthogroup_m
        orthogroup_pcs
        f_ptmdb_m
        f_ivkaphe_m
        ivkaphe_kinases
        f_humap_m
        n
        X_test
        y_test

    main:
        neg_examples = neg_examples( id_dict, n )
        concat_input = neg_examples.collect()

        // split gene pairs in chunks
        examples = cat_and_split( concat_input, 10000 )
                        .flatten()

        // build features from eprot
        eprot_featvec_input = examples.combine( f_eprot_m )
        eprot_featvec_ch = eprot_featvec( eprot_featvec_input ).collect()
        
        // build features from proteomehd
        proteomehd_featvec_input = examples.combine( f_proteomehd_m )
        proteomehd_featvec_ch = proteomehd_featvec( proteomehd_featvec_input )
                                    .collect()

        // build features from mitchell2023
        mitchell2023_featvec_input = examples.combine( f_mitchell2023_m )
        mitchell2023_featvec_ch = mitchell2023_featvec( mitchell2023_featvec_input )
                                    .collect()
        
        // build features from gtex
        gtex_featvec_input = examples.combine( f_gtex_m )
        gtex_featvec_ch = gtex_featvec( gtex_featvec_input )
                            .collect()

        // build features from ptmdb
        ptmdb_featvec_input = examples.combine( f_ptmdb_m )
        ptmdb_featvec_ch = ptmdb_featvec( ptmdb_featvec_input )
                                .collect()

        // build features from depmap
        dependency_featvec_input = examples.combine( f_dependency_m )
        dependency_featvec_ch = dependency_featvec( dependency_featvec_input )
                                    .collect()

        // build features from orthology information
        orthogroup_featvec_input = examples
                                    .combine( f_orthogroup_m )
                                    .combine( orthogroup_pcs )
        orthogroup_featvec_ch = orthogroup_featvec( orthogroup_featvec_input )
                                    .collect()

        // build features from kinase-substrate predictions
        ivkaphe_featvec_input = examples
                                    .combine( f_ivkaphe_m )
                                    .combine( ivkaphe_kinases )
        ivkaphe_featvec_ch = ivkaphe_featvec( ivkaphe_featvec_input )
                                .collect()

        // build features from humap
        humap_score_input = examples.combine( f_humap_m )
        humap_score_ch = humap_score( humap_score_input )
                            .collect()
        
        // join features together
        neg_dataset = features_tables( eprot_featvec_ch,
                                       proteomehd_featvec_ch,
                                       mitchell2023_featvec_ch,
                                       gtex_featvec_ch,
                                       ptmdb_featvec_ch,
                                       dependency_featvec_ch,
                                       orthogroup_featvec_ch,
                                       ivkaphe_featvec_ch,
                                       humap_score_ch )

        dataset = add_examples( neg_dataset, 
                                X_test,
                                y_test )
                                
        X_test = dataset.X_test
        y_test = dataset.y_test
    
    emit:
        X_test
        y_test

}


workflow BUILD_NETWORK {

    take:
        f_gtex_m
        f_eprot_m
        f_proteomehd_m
        f_mitchell2023_m
        f_dependency_m
        f_orthogroup_m
        orthogroup_pcs
        f_ptmdb_m
        f_ubiquitination_m
        f_humap_m
        f_lopit2025_m
        gene_list
        pos_edges
        model
        metabolic_model
    
    main:
        // model
        //model = Channel.fromPath("${out_dir}/training_signalling_xgb/forest_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pkl")

        // generate chunks of gene pairs, map them to UUIDs
        gene_pairs = gene_pairs( gene_list, 200000 )
                        .flatten()
                        .map{ file -> tuple( file.baseName, file ) }

        // eprot features
        eprot_features_input = gene_pairs.combine( f_eprot_m )
        eprot_features = eprot_features( eprot_features_input )
    
        // proteomehd features
        proteomehd_features_input = gene_pairs.combine( f_proteomehd_m )
        proteomehd_features = proteomehd_features( proteomehd_features_input )
        
        // mitchell2023 features
        mitchell2023_features_input = gene_pairs.combine( f_mitchell2023_m )
        mitchell2023_features = mitchell2023_features( mitchell2023_features_input )

        // gtex features
        gtex_features_input = gene_pairs.combine( f_gtex_m )
        gtex_features = gtex_features( gtex_features_input )
        
        // ptmdb features
        ptmdb_features_input = gene_pairs.combine( f_ptmdb_m )
        ptmdb_features = ptmdb_features( ptmdb_features_input )

        // ubiquitination features
        ubiquitination_features_input = gene_pairs.combine( f_ubiquitination_m )
        ubiquitination_features = ubiquitination_features( ubiquitination_features_input )
        
        // dependency features
        dependency_features_input = gene_pairs.combine( f_dependency_m )
        dependency_features = dep_features( dependency_features_input )
        
        // orthogroup features
        orthogroup_features_input = gene_pairs
                                    .combine( f_orthogroup_m )
                                    .combine( orthogroup_pcs )
        orthogroup_features = orthogroup_features( orthogroup_features_input )
        
        // humap features
        humap_features_input = gene_pairs.combine( f_humap_m )
        humap_features = humap_features( humap_features_input )

        // lopit2025 features
        lopit2025_features_input = gene_pairs.combine( f_lopit2025_m )
        lopit2025_features = lopit2025_features( lopit2025_features_input )
        
        // join all features together, by gene pairs chunk's UUID
        edges_features_input = eprot_features
                                    .join( proteomehd_features )
                                    .join( mitchell2023_features )
                                    .join( gtex_features )
                                    .join( ptmdb_features )
                                    .join( dependency_features )
                                    .join( orthogroup_features )
                                    .join( ubiquitination_features )
                                    .join( humap_features )
                                    .join( lopit2025_features )

        // compute edges of the protverse
        edges_features = compute_edge_features( edges_features_input )
        edges_chunks = compute_edge_weights( edges_features, model )

        // concatenate all edges in a single table
        edges_unfiltered = concat_edges(edges_chunks.collect())
        edges_unfiltered_w_id = edges_unfiltered
                                    .flatMap()
                                    .map{file -> tuple( file.baseName, file )}

        // export edges in parquet
        /*TMPedges_parquet = edges_tsv_2_parquet( edges_unfiltered )TMP*/

        // generate chunks of gene lists
        /*ref_genes = get_ref_genes( gene_list, 64 ).flatMap().map{it -> it}

        // extract adjacency vector for each gene
        gene_adjacency_vector_input = edges_parquet.combine( ref_genes )
        adj_vectors = gene_adjacency_vector( gene_adjacency_vector_input )
                        .flatten().collect()

        // compute correlations of adjacency vectors to build the final graph (TODO)    
        //correlate_adjacencies_input = gene_pairs.combine( adj_vectors )
        corr_of_adj = correlate_adjacencies( gene_pairs, adj_vectors )
                        .flatten().collect()

        hypergraph_edges_unfiltered = concatenate( corr_of_adj )*/

        // filter for weight (this generates "wcsn/edges_${params.wholecellnet_edge_min_threshold}minScore.tsv")
        //edges = filter_edge_score( hypergraph_edges_unfiltered )
        edges_filtered = filter_edge_score( edges_unfiltered )

        // transform edge score (this generates "wcsn/edges_${params.wholecellnet_edge_min_threshold}minScore_transformed.tsv")
        edges = transform_edge_score( edges_filtered )

        // add/overwrite positive example edges to protverse
        // this seems actually to lead to worse results
        //edges = add_edges( edges, pos_edges )

        /*TMPw_net_stats( edges_unfiltered_w_id )TMP*/

        /*edges_reactome_genes = filter_reactome_genes( edges,
                                                      gene_list )*/

        // metabolic model edges
        metabolism_edges_chunks = metabolism_compute_edge_weights( edges_features,
                                                                   metabolic_model )
        metabolism_edges_unfiltered = metabolism_concat_edges( metabolism_edges_chunks.collect() )
        // filter edge score (this generates "wcmn/edges_${params.wholecellnet_edge_min_threshold}minScore.tsv")
        metabolism_edges_filtered = metabolism_filter_edge_score( metabolism_edges_unfiltered )
        // transform edge score (this generates "wcmn/edges_${params.wholecellnet_edge_min_threshold}minScore_transformed.tsv")
        metabolism_edges = metabolism_transform_edge_score( metabolism_edges_filtered )


    emit:
        edges
        //edges_parquet
        //edges_reactome_genes
        metabolism_edges

}


workflow PUBLISH_CONFIG {

    main:
        config_ch = Channel.fromPath("${projectDir}/build_protverse.config")
        val_ch = Channel.of('workflow/build_protverse.config')
        
        publish( config_ch, val_ch )

}


workflow {

    // Get gene interaction and omics data from external sources
    reactome = REACTOME()
    gtex = GTEX()
    eprot = EPROT()
    mitchell2023 = MITCHELL2023()
    dependency = DEPENDENCY()
    ptmdb = PTMDB()
    //tahoe = TAHOE()
    lopit2025 = LOPIT2025()
    orthogroup = ORTHOGROUPS()
    ubiquitination = UBIQUITINATION()
    humap3 = HUMAP3()
    proteomehd = PROTEOMEHD()
    uniprot = UNIPROT()

    // Get gene id translations
    id_dict = MAP_IDS( 
        uniprot.gene_list,
        uniprot.idmapping,
        reactome.gene_list,
        reactome.metabolism_gene_list,
        gtex.gene_list,
        eprot.gene_list,
        orthogroup.gene_list,
        dependency.gene_list,
        ubiquitination.gene_list,
        humap3.gene_list,
        ptmdb.gene_list,
        proteomehd.gene_list,
        mitchell2023.gene_list,
        lopit2025.gene_list,
        reactome.hs_gene_list
    )


    // filter and translate datasets to only include tractable gene set
    col_reactome = Channel.from('signalling')
    col_metabolism = Channel.from('metabolism')
    col_humap3 = Channel.from('humap3')
    gene_lists__dict__col = reactome.mapped_cofun_ids
                                .combine(id_dict.unfiltered)
                                .combine(col_reactome)
    metabolism_gene_lists__dict__col = reactome.metabolism_mapped_cofun_ids
                                            .combine(id_dict.unfiltered)
                                            .combine(col_metabolism)
    f_reactome_lists = FILTER_REACTOME( gene_lists__dict__col,
                                        id_dict.unfiltered,
                                        col_reactome )
    f_metabolism_lists = FILTER_METABOLISM( metabolism_gene_lists__dict__col,
                                            id_dict.unfiltered,
                                            col_metabolism)
    f_gtex_m = FILTER_GTEX( gtex.matrix,
                            id_dict.unfiltered,
                            'gtex' )
    f_eprot_m = FILTER_EPROT( eprot.matrix,
                              id_dict.unfiltered,
                              'eprot' )
    f_proteomehd_m = FILTER_PROTEOMEHD( proteomehd.matrix,
                                        id_dict.unfiltered,
                                        'proteomehd' )
    f_mitchell2023_m = FILTER_MITCHELL2023( mitchell2023.matrix,
                                            id_dict.unfiltered,
                                            'mitchell2023' )
    f_lopit2025_m = FILTER_LOPIT2025( lopit2025.table,
                                id_dict.unfiltered,
                                'lopit2025' )
    f_dependency_m = FILTER_DEPENDENCY( dependency.matrix,
                                        id_dict.unfiltered,
                                        'depmap' )
    f_orthogroup_m = FILTER_ORTHOGROUPS( orthogroup.gene2orthogroup,
                                         id_dict.unfiltered,
                                         'orthogroup' )             
    f_ptmdb_m = FILTER_PTMDB( ptmdb.matrix,
                              id_dict.unfiltered,
                              'ptmdb' )
    f_ubiquitination_m = FILTER_UBIQUITINATION( ubiquitination.table,
                                                id_dict.unfiltered,
                                                'ubiquitination' )
    f_humap3_m = FILTER_HUMAP3( humap3.table,
                                id_dict.unfiltered,
                                col_humap3 )


    bernett2024 = BERNETT2024( uniprot.uniprot2gene_name_and_synonym_dict )


    // signalling interaction: generate training data
    dataset = GENERATE_EXAMPLES(
        f_reactome_lists,
        id_dict.unfiltered,
        f_gtex_m,
        f_eprot_m,
        f_proteomehd_m,
        f_mitchell2023_m,
        f_dependency_m,
        f_orthogroup_m,
        orthogroup.pcs,
        f_ptmdb_m,
        f_ubiquitination_m,
        f_humap3_m,
        f_lopit2025_m
    )


    // generate train, validation, and test datasets for physical PPI task
    dataset_bernett2024 = GENERATE_EXAMPLES_BERNETT2024(
        bernett2024.train_ids,
        bernett2024.valid_ids,
        bernett2024.test_ids,
        id_dict.unfiltered,
        f_gtex_m,
        f_eprot_m,
        f_proteomehd_m,
        f_mitchell2023_m,
        f_dependency_m,
        f_orthogroup_m,
        orthogroup.pcs,
        f_ptmdb_m,
        f_ubiquitination_m,
        f_humap3_m,
        f_lopit2025_m
    )


    // metabolic interaction: generate training data
    dataset_metabolism = GENERATE_EXAMPLES_METABOLISM(
        f_metabolism_lists,
        id_dict.unfiltered,
        f_gtex_m,
        f_eprot_m,
        f_proteomehd_m,
        f_mitchell2023_m,
        f_dependency_m,
        f_orthogroup_m,
        orthogroup.pcs,
        f_ptmdb_m,
        f_ubiquitination_m,
        f_humap3_m,
        f_lopit2025_m
    )


    pos_edges = dataset.pos_examples


    // fit and evaluate model on the physical PPI task
    physical_ppi_model = FIT_BERNETT2024(
        dataset_bernett2024.X_train,
        dataset_bernett2024.X_valid,
        dataset_bernett2024.X_test,
        dataset_bernett2024.y_train,
        dataset_bernett2024.y_valid,
        dataset_bernett2024.y_test
    )


    // fit and evaluate single-source models on the physical PPI task
    FIT_BERNETT2024_SINGLE_SOURCE(
        dataset_bernett2024.X_train,
        dataset_bernett2024.X_valid,
        dataset_bernett2024.X_test,
        dataset_bernett2024.y_train,
        dataset_bernett2024.y_valid,
        dataset_bernett2024.y_test
    )


    // fit and evaluate model on the signalling interaction task
    signalling_model = FIT_SIGNALLING(
        dataset.X_train,
        dataset.X_test,
        dataset.y_train,
        dataset.y_test
    )


    // fit and evaluate single-source models on the signalling interaction task
    FIT_SIGNALLING_SINGLE_SOURCE(
        dataset.X_train,
        dataset.X_test,
        dataset.y_train,
        dataset.y_test
    )


    // fit and evaluate model on the metabolic interaction task
    metabolic_model = FIT_METABOLISM(
        dataset_metabolism.X_train,
        dataset_metabolism.X_test,
        dataset_metabolism.y_train,
        dataset_metabolism.y_test
    )


    // fit and evaluate single-source models on the metabolic interaction task
    FIT_METABOLISM_SINGLE_SOURCE(
        dataset_metabolism.X_train,
        dataset_metabolism.X_test,
        dataset_metabolism.y_train,
        dataset_metabolism.y_test
    )

    // build ProtVerse
    protverse = BUILD_NETWORK(
        f_gtex_m,
        f_eprot_m,
        f_proteomehd_m,
        f_mitchell2023_m,
        f_dependency_m,
        f_orthogroup_m,
        orthogroup.pcs,
        f_ptmdb_m,
        f_ubiquitination_m,
        f_humap3_m,
        f_lopit2025_m,
        id_dict.unfiltered,
        pos_edges,
        signalling_model,
        metabolic_model
    )


    // save nextflow.config 
    PUBLISH_CONFIG()

}
