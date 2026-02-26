#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
from a list of genes, generate a list of non redundant pairwise combinations
.META: gene_pairs.tsv
1   gene A name
2   gene B name
*/
process gene_pairs {

    input:
        path 'input/gene_list.txt'
        val n

    output:
        path 'gene_pairs/*.tsv'

    script:
    """
    mkdir -p gene_pairs

    list2combinations.py \\
        input/gene_list.txt \\
        4 \\
        > gene_pairs.tsv

    split -l ${n} -a 16 -x gene_pairs.tsv gene_pairs/ --additional-suffix .tsv
    """    

}

/*
Compute features from E-PROT data

.META: gene_pairs.tsv
1   gene 1
2   gene 2

.META: eprot_input_vector.tsv
1   gene 1
2   gene 2
3   eprot pearson
4   eprot detection ratio gene 1
5   eprot detection ratio gene 2
6   eprot co-detection ratio
*/
process eprot_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '16G'

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/eprot_table.tsv')

    output:
        tuple val(id),
              file('eprot_input_vector.tsv')

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    eprot_features_parallel.py  \
        input/gene_pairs.tsv \
        input/eprot_table.tsv \
        ${params.starmap_n_proc_buildNet} \
        > eprot_input_vector.tsv
    """ 

}

/*
Compute features from gene proteomehd data

.META: gene_pairs.tsv
1   gene 1
2   gene 2

.META: proteomehd_input_vector.tsv
1   gene 1
2   gene 2
3   proteomehd pearson
4   proteomehd co-detection ratio
*/
process proteomehd_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '16G'

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/proteomehd_table.tsv')

    output:
        tuple val(id),
              file('proteomehd_input_vector.tsv')

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    proteomehd_features_parallel.py  \
        input/gene_pairs.tsv \
        input/proteomehd_table.tsv \
        ${params.starmap_n_proc_buildNet} \
        > proteomehd_input_vector.tsv
    """ 

}

/*
Compute features from gene mitchell2023 data

.META: gene_pairs.tsv
1   gene 1
2   gene 2

.META: proteomehd_input_vector.tsv
1   gene 1
2   gene 2
... mitchell2023 features
*/
process mitchell2023_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '16G'

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/mitchell2023_table.tsv')

    output:
        tuple val(id),
              file('mitchell2023_input_vector.tsv')

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    mitchell2023_features_parallel.py  \
        input/gene_pairs.tsv \
        input/mitchell2023_table.tsv \
        ${params.starmap_n_proc_buildNet} \
        > mitchell2023_input_vector.tsv
    """ 

}

/*
Compute features from GTEx data

.META: gene_pairs.tsv
3   gene 1
4   gene 2

.META: gtex_input_vector.tsv
1   gene 1
2   gene 2
3   GTEx pearson
4   GTEx co-detection ratio
*/
process gtex_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '16G'

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/gene_tpm.tsv')

    output:
        tuple val(id),
              file('gtex_input_vector.tsv')

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    gtex_features_parallel.py  \
        input/gene_pairs.tsv \
        input/gene_tpm.tsv \
        ${params.starmap_n_proc_buildNet} \
        > gtex_input_vector.tsv
    """ 

}

/*
Compute features from PTMDB data

NOTE: for two given proteins a and b, corresponding phospho data are picked
      from the pair of phosphosites p_a and p_b -one for each protein-
      according to the following logic:
      - maximize pearson(p_a, p_b) ** 2
      - if pearson(p_a, p_b)==NaN for any (p_a, p_b), then choose (p_a, p_b)
        as the phosphosites detected more frequently

.META: gene_pairs.tsv
1   gene 1
2   gene 2

.META: ptmdb_input_vector.tsv
1   gene 1
2   gene 2
3   PTMDB pearson max
4   PTMDB pearson min
5   PTMDB co-detection ratio
*/
process ptmdb_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '16G'

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/all_logFC.tsv')

    output:
        tuple val(id),
              file('ptmdb_input_vector.tsv')

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    cat input/all_logFC.tsv | sed 's/\t/_/' > all_logFC.tsv

    ptmdb_features_parallel.py  \
        input/gene_pairs.tsv \
        all_logFC.tsv \
        ${params.starmap_n_proc_buildNet} \
        > ptmdb_input_vector.tsv
    """ 

}


/*
Compute features from ubiquitination data

.META: gene_pairs.tsv
1   gene 1
2   gene 2

.META: ptmdb_input_vector.tsv
1   gene 1
2   gene 2
[...]   [features]
*/
process ubiquitination_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '16G'

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/all_logFC.tsv')

    output:
        tuple val(id),
              file('ubiquitination_input_vector.tsv')

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    cat input/all_logFC.tsv | sed 's/\t/_/' > all_logFC.tsv

    ubiquitination_features_parallel.py  \
        input/gene_pairs.tsv \
        all_logFC.tsv \
        ${params.starmap_n_proc_buildNet} \
        > ubiquitination_input_vector.tsv
    """ 

}

/*
Compute features from gene dependency data

.META: gene_pairs.tsv
1   gene 1
2   gene 2

.META: dependency_input_vector.tsv
1   gene 1
2   gene 2
[...]   [features]
*/
process dep_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '20G'

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/dependency_table.tsv')

    output:
        tuple val(id),
              file('dependency_input_vector.tsv')

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    dep_features_parallel.py  \
        input/gene_pairs.tsv \
        input/dependency_table.tsv \
        ${params.starmap_n_proc_buildNet} \
        > dependency_input_vector.tsv
    """ 

}

/*
Compute features from orthogroup data

.META: gene_pairs.tsv
1   gene 1
2   gene 2

.META: orthogroup_input_vector.tsv
1   gene 1
2   gene 2
3   highest orthogroup pcs for the gene pair
4   bool for same orthogroup
*/
process orthogroup_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '32G'

    debug "false"

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/orthogroup_table.tsv'),
              file('input/orthogroup_pcs.tsv')

    output:
        tuple val(id),
              file('orthogroup_input_vector.tsv')

    script:
    """
    orthogroup_features_parallel.py  \
        input/gene_pairs.tsv \
        input/orthogroup_table.tsv \
        input/orthogroup_pcs.tsv \
        ${params.starmap_n_proc_buildNet} \
        > orthogroup_input_vector.tsv
    """ 

}

/*
Compute features from IV-KAPHE data

.META: gene_pairs.tsv
3   gene 1
4   gene 2

.META: ivkaphe_input_vector.tsv
1   gene 1
2   gene 2
3   ivkaphe max score
4   score exists
*/
process ivkaphe_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '16G'

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/ivkaphe_table.tsv'),
              file('input/kinases.txt')

    output:
        tuple val(id),
              file('ivkaphe_input_vector.tsv')

    script:
    """
    ivkaphe_features_parallel.py  \
        input/gene_pairs.tsv \
        input/ivkaphe_table.tsv \
        input/kinases.txt \
        ${params.starmap_n_proc_buildNet} \
        > ivkaphe_input_vector.tsv
    """ 

}

/*
Compute features from hu.MAP data

.META: gene_pairs.tsv
3   gene 1
4   gene 2

.META: gtex_input_vector.tsv
1   gene 1
2   gene 2
3   humap_score
4   humap_reported
*/
process humap_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '32G'

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/humap.tsv')

    output:
        tuple val(id),
              file('humap_input_vector.tsv')

    script:
    """
    humap_features_parallel.py  \
        input/gene_pairs.tsv \
        input/humap.tsv \
        ${params.starmap_n_proc_buildNet} \
        > humap_input_vector.tsv
    """ 

}


/*
Compute features from cotranslocation data

.META: gene_pairs.tsv
1   gene 1
2   gene 2

.META: lopit2025_input_vector.tsv
1   gene 1
2   gene 2
[...]   [features]
*/
process lopit2025_features {

    cpus "${params.starmap_n_proc_buildNet}"
    memory '20G'

    input:
        tuple val(id),
              file('input/gene_pairs.tsv'),
              file('input/lopit2025_table.tsv')

    output:
        tuple val(id),
              file('lopit2025_input_vector.tsv')
              
    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    lopit2025_features_parallel.py  \
        input/gene_pairs.tsv \
        input/lopit2025_table.tsv \
        ${params.starmap_n_proc_buildNet} \
        > lopit2025_input_vector.tsv
    """ 

}


/*
Merge the partial input vectors from each of the databases

.META: wcsn/edges.tsv.gz
1   index (<Uniprot AC gene 1>_<Uniprot AC gene 2>)
2   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
3   <feature 1>
4   <feature 2>
5   <feature 3>
... ...
*/
process compute_edge_features {

    memory '16G'

    input:
        tuple val(id),
              file('input/eprot.tsv'),
              file('input/proteomehd.tsv'),
              file('input/mitchell2023.tsv'),
              file('input/gtex.tsv'),
              file('input/ptmdb.tsv'),
              file('input/dependency.tsv'),
              file('input/orthogroup.tsv'),
              file('input/ubiquitination.tsv'),
              file('input/humap.tsv'),
              file('input/lopit2025.tsv')

    output:
        path "features.tsv"

    script:
    """
    mkdir -p edges_features

    cat \
        <(cat input/gtex.tsv | sed -n '1p' | cut -f3- \
            | awk '{print "index\\t"\$0}') \
        <(cat input/gtex.tsv \
            | sed '1d' | grep -v "label" \
            | awk 'NF{print \$1"_"\$2"\\t"\$0}' | cut -f1,4-) \
        > edges_features/gtex.tsv

    cat \
        <(cat input/eprot.tsv | sed -n '1p' | cut -f3- \
            | awk '{print "index\\t"\$0}') \
        <(cat input/eprot.tsv \
            | sed '1d' | grep -v "label" \
            | awk 'NF{print \$1"_"\$2"\\t"\$0}' | cut -f1,4-) \
        > edges_features/eprot.tsv

    cat \
        <(cat input/proteomehd.tsv | sed -n '1p' | cut -f3- \
            | awk '{print "index\\t"\$0}') \
        <(cat input/proteomehd.tsv \
            | sed '1d' | grep -v "label" \
            | awk 'NF{print \$1"_"\$2"\\t"\$0}' | cut -f1,4-) \
        > edges_features/proteomehd.tsv

    cat \
        <(cat input/mitchell2023.tsv | sed -n '1p' | cut -f3- \
            | awk '{print "index\\t"\$0}') \
        <(cat input/mitchell2023.tsv \
            | sed '1d' | grep -v "label" \
            | awk 'NF{print \$1"_"\$2"\\t"\$0}' | cut -f1,4-) \
        > edges_features/mitchell2023.tsv
    
    cat \
        <(cat input/ptmdb.tsv | sed -n '1p' | cut -f3- \
            | awk '{print "index\\t"\$0}') \
        <(cat input/ptmdb.tsv \
            | sed '1d' | grep -v "label" \
            | awk 'NF{print \$1"_"\$2"\\t"\$0}' | cut -f1,4-) \
        > edges_features/ptmdb.tsv

    cat \
        <(cat input/ubiquitination.tsv | sed -n '1p' | cut -f3- \
            | awk '{print "index\\t"\$0}') \
        <(cat input/ubiquitination.tsv \
            | sed '1d' | grep -v "label" \
            | awk 'NF{print \$1"_"\$2"\\t"\$0}' | cut -f1,4-) \
        > edges_features/ubiquitination.tsv

    cat \
        <(cat input/orthogroup.tsv | sed -n '1p' | cut -f3- \
            | awk '{print "index\\t"\$0}') \
        <(cat input/orthogroup.tsv \
            | sed '1d' | grep -v "label" \
            | awk 'NF{print \$1"_"\$2"\\t"\$0}' | cut -f1,4-) \
        > edges_features/orthogroup.tsv

    cat \
        <(cat input/dependency.tsv | sed -n '1p' | cut -f3- \
            | awk '{print "index\\t"\$0}') \
        <(cat input/dependency.tsv \
            | sed '1d' | grep -v "label" \
            | awk 'NF{print \$1"_"\$2"\\t"\$0}' | cut -f1,4-) \
        > edges_features/dependency.tsv

    cat \
        <(cat input/humap.tsv | sed -n '1p' | cut -f3- \
            | awk '{print "index\\t"\$0}') \
        <(cat input/humap.tsv \
            | sed '1d' | grep -v "label" \
            | awk 'NF{print \$1"_"\$2"\\t"\$0}' | cut -f1,4-) \
        > edges_features/humap.tsv

    cat \
        <(cat input/lopit2025.tsv | sed -n '1p' | cut -f3- \
            | awk '{print "index\\t"\$0}') \
        <(cat input/lopit2025.tsv \
            | sed '1d' | grep -v "label" \
            | awk 'NF{print \$1"_"\$2"\\t"\$0}' | cut -f1,4-) \
        > edges_features/lopit2025.tsv

    edges_features.py \
        edges_features/eprot.tsv \
        edges_features/proteomehd.tsv \
        edges_features/mitchell2023.tsv \
        edges_features/gtex.tsv \
        edges_features/ptmdb.tsv \
        edges_features/ubiquitination.tsv \
        edges_features/dependency.tsv \
        edges_features/orthogroup.tsv \
        edges_features/humap.tsv \
        edges_features/lopit2025.tsv \
        > features.tsv
    """

}

/*
given a list of gene pairs with features, and a model, compute the
predicted signalling functional relatedness score. 

.META:
1   gene_a
2   gene_b
3   edge score
*/
process compute_edge_weights {

    memory '32G'

    input:
        path 'input/features.tsv'
        path 'input/model.pkl'

    output:
        path 'edges.tsv'

    script:
    """
    computeEdges.py \
        input/features.tsv \
        input/model.pkl \
    | sed 's/_/\\t/' \
    > edges.tsv
    """

}


/*
concatenate all edges of the wcsn (whole-cell signalling network)
*/
process concat_edges {

    publishDir "${out_dir}",
                pattern: "wcsn/edges.tsv",
                mode: 'copy'

    input:
        path 'input/*.tsv'

    output:
        path "wcsn/edges.tsv"

    script:
    """
    mkdir -p wcsn

    cat input/*.tsv | awk 'NF' | sort \
        > wcsn/edges.tsv
    """

}


/*
concatenate all edges of the wcmn (whole-cell metabolic network)
*/
process metabolism_concat_edges {

    publishDir "${out_dir}",
                pattern: "wcmn/edges.tsv",
                mode: 'copy'

    input:
        path 'input/*.tsv'

    output:
        path "wcmn/edges.tsv"

    script:
    """
    mkdir -p wcmn

    cat input/*.tsv | awk 'NF' | sort \
        > wcmn/edges.tsv
    """

}




/*
generate parquet format of tsv edge list
*/
process edges_tsv_2_parquet {

    memory "${params.adj_dask_memory}"
    cpus "${params.adj_dask_n_proc}"

    input:
        path 'input/file.tsv'

    output:
        path "edges.parquet.tar.gz"

    script:
    """
    edges_tsv_2_parquet.py

    tar -zcvf edges.parquet.tar.gz edges.parquet
    """

}


/*
get list of reference genes, split nito chucks, UUID-named
*/
process get_ref_genes {

    input:
        path 'input/gene_list.txt'
        val n

    output:
        path "ref_genes/*txt"

    script:
    """
    mkdir -p ref_genes

    cat input/gene_list.txt \
        | grep -v "gene_name"  | cut -f1 | awk 'NF' \
        | sort | uniq \
        > genes.txt

    split -l ${n} -a 16 -x genes.txt ref_genes/ --additional-suffix .txt
    """

}


/*
get the edge weights for a given gene, for a list of genes
*/
process gene_adjacency_vector {

    memory "${params.adj_dask_memory}"
    cpus "${params.adj_dask_n_proc}"

    input:
        tuple file('input/edges.parquet.tar.gz'),
              file('input/genes.txt')

    output:
        path "adj_vectors/*.tsv"

    script:
    """
    mkdir -p adj_vectors

    tar -zxvf input/edges.parquet.tar.gz

    adjacency_vector.py \
        edges.parquet \
        input/genes.txt \
        adj_vectors    
    """

}


/*
compute Pearson correlation coefficient between pairs of genes;
the vectors are given by the weights in the adjacency lists of those genes
*/
process correlate_adjacencies {

    memory "32G"
    cpus "${params.corr_adj_dask_n_proc}"

    input:
        tuple val(id),
              file('input/gene_pairs.tsv')
        path 'input/adj/*.tsv'

    output:
        path "wcsn/adj_corr_edges/*.tsv"

    script:
    """
    mkdir -p wcsn/adj_corr_edges adj

    #target_file=\$(find input/adj -type l -exec readlink {} \\; -quit)
    #adj_dir=\$(dirname \${target_file})

    for file in \$(ls input/adj); do \
        target_file=\$(readlink input/adj/"\$file"); \
        target_basename=\$(basename "\$target_file"); \
        ln -s "\$target_file" "adj/\$target_basename"; \
    done

    correlate_adjacencies.py \
        input/gene_pairs.tsv \
        adj \
        > wcsn/adj_corr_edges/${id}.tsv
    """

}


/*
filter wcsn to only include edges having
weight > ${params.wholecellnet_edge_min_threshold}
and rescale them by substracting ${params.wholecellnet_edge_min_threshold} 
and dividing by 1 - ${params.wholecellnet_edge_min_threshold}
*/
process filter_edge_score {

    publishDir "${out_dir}",
                pattern: "wcsn/edges_${params.wholecellnet_edge_min_threshold}minScore.tsv",
                mode: 'copy'
    
    input:
        path 'input/edges.tsv'
    
    output:
        path "wcsn/edges_${params.wholecellnet_edge_min_threshold}minScore.tsv"
    
    script:
    """
    mkdir -p wcsn
    
    cat input/edges.tsv \\
    | awk '\$3>${params.wholecellnet_edge_min_threshold}{S=(\$3-${params.wholecellnet_edge_min_threshold})/(1-${params.wholecellnet_edge_min_threshold}); print \$1"\\t"\$2"\\t"S}' \\
        > wcsn/edges_${params.wholecellnet_edge_min_threshold}minScore.tsv
    """

}


/*
filter wcmn to only include edges having
weight > ${params.wholecellnet_edge_min_threshold}
and rescale them by substracting ${params.wholecellnet_edge_min_threshold} 
and dividing by 1 - ${params.wholecellnet_edge_min_threshold}
*/
process metabolism_filter_edge_score {

    publishDir "${out_dir}",
                pattern: "wcmn/edges_${params.wholecellnet_edge_min_threshold}minScore.tsv",
                mode: 'copy'
    
    input:
        path 'input/edges.tsv'
    
    output:
        path "wcmn/edges_${params.wholecellnet_edge_min_threshold}minScore.tsv"
    
    script:
    """
    mkdir -p wcmn
    
    cat input/edges.tsv \\
    | awk '\$3>${params.wholecellnet_edge_min_threshold}{S=(\$3-${params.wholecellnet_edge_min_threshold})/(1-${params.wholecellnet_edge_min_threshold}); print \$1"\\t"\$2"\\t"S}' \\
        > wcmn/edges_${params.wholecellnet_edge_min_threshold}minScore.tsv
    """

}


/*
transform edge scores by applying the specified function
*/
process transform_edge_score {

    publishDir "${out_dir}",
                pattern: "wcsn/*.tsv",
                mode: 'copy'
    
    input:
        path 'input/edges.tsv'
    
    output:
        path "wcsn/edges_${params.wholecellnet_edge_min_threshold}minScore_transformed.tsv"
    
    script:
    """
    mkdir -p wcsn
    
    cat input/edges.tsv \\
        | awk '{s=exp(2*log(\$3)); print \$1"\\t"\$2"\\t"s}' \\
        > wcsn/edges_${params.wholecellnet_edge_min_threshold}minScore_transformed.tsv
    """

}


/*
transform edge scores by applying the specified function
*/
process metabolism_transform_edge_score {

    publishDir "${out_dir}",
                pattern: "wcmn/*.tsv",
                mode: 'copy'
    
    input:
        path 'input/edges.tsv'
    
    output:
        path "wcmn/edges_${params.wholecellnet_edge_min_threshold}minScore_transformed.tsv"
    
    script:
    """
    mkdir -p wcmn
    
    cat input/edges.tsv \\
    | awk '{s=exp(2*log(\$3)); print \$1"\\t"\$2"\\t"s}' \\
        > wcmn/edges_${params.wholecellnet_edge_min_threshold}minScore_transformed.tsv
    """

}


/*
add edges specified in separate lists
only add non-redundant edges, overwriting existing identical edges
weights are set to 1 for all the newly added edges
*/
process add_edges {

    publishDir "${out_dir}",
                pattern: "wcsn/edges_refined_${params.wholecellnet_edge_min_threshold}minScore.tsv",
                mode: 'copy'

    input:
        path 'input/wcsn.tsv'
        path 'input/edges.tsv'

    output:
        path "wcsn/edges_refined_${params.wholecellnet_edge_min_threshold}minScore.tsv"

    script:
    """
    mkdir -p wcsn

    cat input/edges.tsv | cut -f3,4 > to_add.tsv
    cat input/edges.tsv | cut -f3,4 | awk '{print \$2"\\t"\$1}' > to_add_rev.tsv
    cat to_add.tsv to_add_rev.tsv > to_remove.tsv
    cat input/wcsn.tsv | grep -v -f <(cat to_remove.tsv) > wcsn.tsv
    cat to_add.tsv | awk '{print \$0"\\t1"}' > new_edges.tsv
    cat wcsn.tsv new_edges.tsv | sort \
        > wcsn/edges_refined_${params.wholecellnet_edge_min_threshold}minScore.tsv
    """

}

/*
make some plots showing the properties of a weighted graph
*/
process w_net_stats {

    memory "32G"

    publishDir "${out_dir}",
                pattern: "net_stats/${id}/edge_weight.pdf",
                mode: 'copy'

    input:  
        tuple val(id),
              file('input/net.tsv')

    output:
        path "net_stats/${id}/edge_weight.pdf"

    script:
    """
    mkdir -p net_stats/${id}

    w_net_stats.py \
        input/net.tsv \
        net_stats/${id}/edge_weight.pdf 
    """

}


/*
filter edges by only picking the ones where
both genes are annotated in Reactome
*/
process filter_reactome_genes {

    publishDir "${out_dir}",
                pattern: "wcsn/*.tsv",
                mode: 'copy'
    
    input:
        path 'input/edges.tsv'
        path 'input/gene_dict.tsv'

    output:
        path "wcsn/edges_reactomeGenes_${params.wholecellnet_edge_min_threshold}minScore.tsv"
    
    script:
    """
    mkdir -p wcsn

    cat input/gene_dict.tsv \\
        | cut -f1,2 | sed '1d' | awk '\$2!="NA"{print \$1}' \
        > gene_list.txt

    awk 'NR==FNR {genes[\$1]; next} (\$1 in genes) && (\$2 in genes)' \\
        gene_list.txt input/edges.tsv \\
        > wcsn/edges_reactomeGenes_${params.wholecellnet_edge_min_threshold}minScore.tsv
    """

}