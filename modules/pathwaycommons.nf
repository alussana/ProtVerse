#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
download intact interaction data in text format from
PathwayCommons <https://www.pathwaycommons.org>
*/
process dl_pc13_intact {

    publishDir "${out_dir}",
                pattern: "databases/PathwayCommons12/intact.hgnc.sif.gz",
                mode: 'copy'

    output:
        path 'databases/PathwayCommons12/intact.hgnc.sif.gz'

    script:
    """
    mkdir -p databases/PathwayCommons12

    wget -O databases/PathwayCommons12/intact.hgnc.sif.gz \
        ${params.url_pc13_intact}
    """
}

/*
download biogrid interaction data in text format from
PathwayCommons <https://www.pathwaycommons.org>
*/
process dl_pc13_biogrid {

    publishDir "${out_dir}",
                pattern: "databases/PathwayCommons12/biogrid.hgnc.sif.gz",
                mode: 'copy'

    output:
        path 'databases/PathwayCommons12/biogrid.hgnc.sif.gz'

    script:
    """
    mkdir -p databases/PathwayCommons12

    wget -O databases/PathwayCommons12/biogrid.hgnc.sif.gz \
        ${params.url_pc13_biogrid}
    """
}

/*
obtain unweighted edges from intact interaction data

.META:
1   interactor A
2   interactor B
*/
process unweighted_intact {

    publishDir "${out_dir}",
                pattern: "databases/PathwayCommons12/intact_edges.tsv",
                mode: 'copy'

    input:
        path 'input/pc13_intact.sif.gz'

    output:
        path 'databases/PathwayCommons12/intact_edges.tsv'

    script:
    """
    mkdir -p databases/PathwayCommons12

    zcat input/pc13_intact.sif.gz \
        | cut -f1,3 | sed '1d' | awk 'NF' \
        > databases/PathwayCommons12/intact_edges.tsv
    """
}

/*
obtain unweighted edges from biogrid interaction data

.META:
1   interactor A
2   interactor B
*/
process unweighted_biogrid {

    publishDir "${out_dir}",
                pattern: "databases/PathwayCommons12/biogrid_edges.tsv",
                mode: 'copy'

    input:
        path 'input/pc13_biogrid.sif.gz'

    output:
        path 'databases/PathwayCommons12/biogrid_edges.tsv'

    script:
    """
    mkdir -p databases/PathwayCommons12

    zcat input/pc13_biogrid.sif.gz \
        | cut -f1,3 | sed '1d' | awk 'NF' \
        > databases/PathwayCommons12/biogrid_edges.tsv
    """
}

/*
download reactome edges in sif format for hgnc nomenclature

.META:
1   source gene name    A2M
2   interaction type    controls-state-change-of
3   target gene name    RAC1

353636 interactions are distributed in 12 types:

```
  52989 catalysis-precedes
  13072 chemical-affects
   5173 consumption-controlled-by
   3757 controls-expression-of
   3415 controls-phosphorylation-of
   5343 controls-production-of
 116870 controls-state-change-of
   5325 controls-transport-of
   1927 controls-transport-of-chemical
 141797 in-complex-with
    599 reacts-with
   3369 used-to-produce
*/
process dl_pc12_reactome {

    publishDir "${out_dir}",
            pattern: 'databases/pathwaycommons/PathwayCommons12.reactome.hgnc.sif.gz',
            mode: 'copy'

    output:
        path 'databases/pathwaycommons/PathwayCommons12.reactome.hgnc.sif.gz'

    shell:
    """
    mkdir -p databases/pathwaycommons

    wget -O databases/pathwaycommons/PathwayCommons12.reactome.hgnc.sif.gz \
        '${params.url_pc12_reactome_hgnc_sif}'
    """

}

/*
filter for controls-phosphorylation-of and controls-state-change-of
interaction types

.META: databases/pathwaycommons/signalling_specific.hgnc.sif.gz
1   source gene name    A2M
2   interaction type    controls-state-change-of
3   target gene name    RAC1

.META: databases/pathwaycommons/signalling_gene_occurrences.tsv
1   gene name
2   numbe of total occurrences in databases/pathwaycommons/signalling_specific.hgnc.sif.gz
*/
process filter_sig_interactions {

    publishDir "${out_dir}",
            pattern: 'databases/pathwaycommons/signalling_specific.hgnc.sif.gz',
            mode: 'copy'

    publishDir "${out_dir}",
            pattern: 'databases/pathwaycommons/signalling_gene_occurrences.tsv',
            mode: 'copy'

    input:
        path 'input/PathwayCommons12.reactome.hgnc.sif.gz'

    output:
        path 'databases/pathwaycommons/signalling_specific.hgnc.sif.gz', emit: sif
        path 'databases/pathwaycommons/signalling_gene_occurrences.tsv'

    shell:
    """
    mkdir -p databases/pathwaycommons

    zcat input/PathwayCommons12.reactome.hgnc.sif.gz \
        | grep -f \
            <(echo -e "controls-phosphorylation-of\\ncontrols-state-change-of") \
        | gzip \
        > databases/pathwaycommons/signalling_specific.hgnc.sif.gz

    zcat databases/pathwaycommons/signalling_specific.hgnc.sif.gz \
        | cut -f1,3 | tr '\\t' '\\n' | sort | uniq -c | tr -s " " \
        | sed 's/^ //' | awk '{print \$2"\\t"\$1}' \
        > databases/pathwaycommons/signalling_gene_occurrences.tsv
    """

}

/*
get reactome signalling-specific, unweighted genes
*/
process get_reactome_sig_edges {

    publishDir "${out_dir}",
            pattern: 'databases/pathwaycommons/reactome/net.tsv',
            mode: 'copy'

    input:
        path 'input/edges.sif.gz'

    output:
        path 'databases/pathwaycommons/reactome/net.tsv'

    shell:
    """
    mkdir -p databases/pathwaycommons/reactome

    zcat input/edges.sif.gz \
        | cut -f1,3 | sort | uniq \
        > databases/pathwaycommons/reactome/net.tsv
    """

}

/*
get reactome sif genes
*/
process get_reactome_genes {

    publishDir "${out_dir}",
            pattern: 'databases/pathwaycommons/reactome/genes.txt',
            mode: 'copy'

    input:
        path 'input/edges.sif.gz'

    output:
        path 'databases/pathwaycommons/reactome/genes.txt'

    shell:
    """
    mkdir -p databases/pathwaycommons/reactome

    zcat input/edges.sif.gz \
        | cut -f1,3 | tr '\\t' '\\n' | sort | uniq \
        > databases/pathwaycommons/reactome/genes.txt
    """

}


/*
get reactome signalling-specific, unweighted genes
*/
process make_translate_reactome_edges {

    publishDir "${out_dir}",
            pattern: 'databases/pathwaycommons/reactome/*.tsv.gz',
            mode: 'copy'

    input:
        path 'input/edges.sif.gz'
        path 'input/dict.tsv'

    output:
        path 'databases/pathwaycommons/reactome/edges_translated.tsv.gz'

    shell:
    """
    mkdir -p databases/pathwaycommons/reactome

    zcat input/edges.sif.gz \\
        | cut -f1,3 | sort -u \\
        | grep -v -f <(echo -e ":") \\
        > databases/pathwaycommons/reactome/edges.tsv

    tr_fast.py \\
        input/dict.tsv \\
        databases/pathwaycommons/reactome/edges.tsv \\
        1 \\
        3 \\
        1,2 \\
        1 \\
        > databases/pathwaycommons/reactome/edges_translated_wo_header.tsv 

    cat \\
        <(echo -e "source\ttarget") \\
        databases/pathwaycommons/reactome/edges_translated_wo_header.tsv \\
        | gzip > databases/pathwaycommons/reactome/edges_translated.tsv.gz 

    cat \\
        <(echo -e "source\ttarget") \\
        databases/pathwaycommons/reactome/edges.tsv \\
        | gzip > databases/pathwaycommons/reactome/edges_translated.tsv.gz
    """

}
