#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Download STRING full protein network data with scored links between homo
sapiens proteins (combined score)
*/
process dl_stringdb_net {

    publishDir "${out_dir}", pattern: "databases/stringdb/*", mode: 'copy'

    output:
        path 'databases/stringdb/9606.protein.links.v11.5.txt.gz'

    script:
    """
    mkdir -p databases/stringdb/
    wget -O databases/stringdb/9606.protein.links.v11.5.txt.gz \
        ${params.url_string_hs_links}
    """
}

/*
Download STRING full protein network data with scored links between homo
sapiens proteins (detailed scores)
*/
process dl_stringdb_net_detailed {

    publishDir "${out_dir}", pattern: "databases/stringdb/*", mode: 'copy'

    output:
        path 'databases/stringdb/9606.protein.detailed.links.v11.5.txt.gz'

    script:
    """
    mkdir -p databases/stringdb/
    wget -O databases/stringdb/9606.protein.detailed.links.v11.5.txt.gz \
        ${params.url_string_hs_links_detailed}
    """
}

/*
Filter STRING full protein network data with scored links between homo
sapiens proteins, to only inlcude links with combined score greater than
${params.string_combined_score_treshold}
Also remove the organism identifier prepended to the ENSP id

.META:
1   protein 1 ENSP id
2   protein 2 ENSP id
3   combined score
*/
process parse_stringdb {

    publishDir "${out_dir}", pattern: "databases/stringdb/*", mode: 'copy'

    input:
        path 'input/stringdb.gz'

    output:
        path "databases/stringdb/links_scoreGt${params.string_combined_score_treshold}.tsv"

    script:
    """
    mkdir -p databases/stringdb/
    zcat input/stringdb.gz | sed '1d' \
        | awk '\$3>${params.string_combined_score_treshold}' \
        | sed 's/....\\.//g' | tr " " "\t" \
        > databases/stringdb/links_scoreGt${params.string_combined_score_treshold}.tsv
    """
}

/*
Filter STRING full protein network data with scored links between homo
sapiens proteins, to only inlcude links with combined score greater than
${params.string_combined_score_treshold}
Also remove the organism identifier prepended to the ENSP id

.META:
1   protein 1 ENSP id
2   protein 2 ENSP id
3   stringdb_neighborhood
4   stringdb_fusion
5   stringdb_cooccurence
6   stringdb_coexpression
7   stringdb_experimental
8   stringdb_database
9   stringdb_textmining
10  stringdb_combined_score
*/
process parse_stringdb_detailed {

    publishDir "${out_dir}", pattern: "databases/stringdb/*", mode: 'copy'

    input:
        path 'input/stringdb_detailed.gz'

    output:
        path "databases/stringdb/links_detailed.tsv"

    script:
    """
    mkdir -p databases/stringdb/
    zcat input/stringdb_detailed.gz | sed '1d' \
        | sed 's/....\\.//g' | tr " " "\t" \
        > databases/stringdb/links_detailed.tsv
    """
}


/*
Get list of ENSP id in occurring in the parsed network
*/
process get_stringdb_genes {

    input:
        path 'input/stringdb.tsv'

    output:
        path 'genes.txt'

    script:
    """
    cat input/stringdb.tsv | tr '\t' '\n' | sort | uniq \
        > genes.txt
    """
}

/*
Build the partial input vector from STRINGDB data (combined score)

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0)
3   gene 1
4   gene 2

.META: stringdb_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0)
3   gene 1
4   gene 2
5   stringdb combined score
*/
process stringdb_score {

    input:
        tuple val(gene_pairs),
              file('input/stringdb_table.tsv')

    output:
        path 'stringdb_input_vector.tsv'

    script:
    """
    echo -e "${gene_pairs}" | awk 'NF' > gene_pairs.tsv

    stringdb_input_vector.py  \
        gene_pairs.tsv \
        input/stringdb_table.tsv \
        > stringdb_input_vector.tsv
    """ 

}

/*
Build the partial input vector from STRINGDB data (detailed scores)

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0)
3   gene 1
4   gene 2

.META: stringdb_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0)
3   gene 1
4   gene 2
5   stringdb_neighborhood
6   stringdb_fusion
7   stringdb_cooccurence
8   stringdb_coexpression
9   stringdb_experimental
10  stringdb_database
11  stringdb_textmining
*/
process stringdb_score_detailed {

    memory '16G'
    
    input:
        tuple val(gene_pairs),
              file('input/stringdb_table_detailed.tsv')

    output:
        path 'stringdb_input_vector.tsv'

    script:
    """
    echo -e "${gene_pairs}" | awk 'NF' > gene_pairs.tsv

    cat input/stringdb_table_detailed.tsv | grep -v "ENSP" \
        > stringdb_table_detailed.tsv

    stringdb_input_vector_detailed.py  \
        gene_pairs.tsv \
        stringdb_table_detailed.tsv \
        > stringdb_input_vector.tsv
    """ 

}