#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Download Reactome reactions IDs and metadata
- Identifier mapping files, Lowest level pathway diagram / Subset of the pathway
- Pathway Information, Pathways hierarchy relationship
- Pathway Information, Complete List of Pathways
*/
process dl_reactome {

    publishDir "${out_dir}", pattern: "databases/reactome/*", mode: 'copy'

    output:
        path 'databases/reactome/ReactomePathways.txt',
            emit: pathway_ids
        path 'databases/reactome/ReactomePathwaysRelation.txt',
            emit: hierarchy
        path 'databases/reactome/Ensembl2Reactome.txt',
            emit: uniprot_ids

    script:
    """
    mkdir -p databases/reactome/
    wget -O databases/reactome/ReactomePathways.txt ${params.url_ReactomePathways}
    wget -O databases/reactome/ReactomePathwaysRelation.txt ${params.url_ReactomePathwaysRelation}
    wget -O databases/reactome/Ensembl2Reactome.txt ${params.url_Ensembl2Reactome}
    """
}

/*
download Reactome_2022 from the EnrichR libraries
<https://maayanlab.cloud/Enrichr/#libraries>
*/
process dl_reactome2022_enrichr {

    publishDir "${out_dir}",
               pattern: "databases/reactome/*tsv",
               mode: 'copy'

    output:
        tuple val('Reactome_2022'),
              file('databases/reactome/enrichr_Reactome_2022.tsv')

    script:
    """
    mkdir -p databases/reactome
    wget \
        -O databases/reactome/enrichr_Reactome_2022.tsv \
        '${params.url_reactome2022_enrichr}'
    """

}

/*
plot modules size distribution etc
*/
process modules_stats {

    publishDir "${out_dir}",
               pattern: "modules/${net}/*.pdf",
               mode: 'copy'

    input:
        tuple val(net),
              file('input/modules.tsv')

    output:
        path "modules/${net}/logSizeDistrib.pdf"
        path "modules/${net}/overlapRatioDistrib.pdf"
        path "modules/${net}/pleiotropyDistrib.pdf"
        path "modules/${net}/sizeDistrib.pdf"

    script:
    """
    mkdir -p modules/${net}

    cat input/modules.tsv | cut -f2- > modules.tsv

    modules_stats.py \
        modules.tsv \
        modules/${net}/logSizeDistrib.pdf \
        modules/${net}/overlapRatioDistrib.pdf \
        modules/${net}/pleiotropyDistrib.pdf \
        modules/${net}/sizeDistrib.pdf
    """

}

/*
== NOT USED ==

Download Reactome pairwise human interactions in tabular format
.META: https://reactome.org/download-data?id=62
*/
process dl_reactome_tab {

    publishDir "${out_dir}", pattern: "databases/reactome/*", mode: 'copy'

    output:
        path 'databases/reactome/reactome.homo_sapiens.interactions.psi-mitab.txt'

    script:
    """
    mkdir -p databases/reactome/
    wget -O databases/reactome/reactome.homo_sapiens.interactions.psi-mitab.txt \
        ${params.url_ReactomeInteractions}
    """
}

/*
Get a list of pathway IDs for all the reactions rooted at signal transduction
(i.e. node ID ${params.ReactomeSignalTransduction_id}) and related to
Homo sapiens
*/
process sig_react_ids {

    publishDir "${out_dir}", 
        pattern: "databases/reactome/hs_sig_react_ids.txt",
        mode: 'copy'

    input:
        path 'input/ReactomePathways.txt'
        path 'input/ReactomePathwaysRelation.txt'
        path 'input/Ensembl2Reactome.txt'

    output:
        path 'databases/reactome/hs_sig_react_ids.txt'

    script:
    """
    mkdir -p databases/reactome/
    
    visitReactomeGraph.py \
        input/ReactomePathwaysRelation.txt \
        ${params.ReactomeSignalTransduction_id} \
        | sort | uniq \
        > sig_react_ids.txt
    
    cat input/Ensembl2Reactome.txt \
        | awk -F "\t" '\$6~/Homo sapiens/{print \$2}' \
        | sort | uniq \
        > hs_path_ids.txt
    
    cat hs_path_ids.txt | grep -f sig_react_ids.txt \
        > databases/reactome/hs_sig_react_ids.txt
    """

}


/*
Get a list of pathway IDs for all the reactions rooted at Metabolism
(i.e. node ID ${params.ReactomeMetabolism_id}) and related to
Homo sapiens
*/
process metabolism_react_ids {

    publishDir "${out_dir}", 
        pattern: "databases/reactome/metabolism_hs_react_ids.txt",
        mode: 'copy'

    input:
        path 'input/ReactomePathways.txt'
        path 'input/ReactomePathwaysRelation.txt'
        path 'input/Ensembl2Reactome.txt'

    output:
        path 'databases/reactome/metabolism_hs_react_ids.txt'

    script:
    """
    mkdir -p databases/reactome/
    
    visitReactomeGraph.py \
        input/ReactomePathwaysRelation.txt \
        ${params.ReactomeMetabolism_id} \
        | sort | uniq \
        > metabolism_react_ids.txt
    
    cat input/Ensembl2Reactome.txt \
        | awk -F "\t" '\$6~/Homo sapiens/{print \$2}' \
        | sort | uniq \
        > metabolism_hs_path_ids.txt

    cat metabolism_hs_path_ids.txt | grep -f metabolism_react_ids.txt \
        > databases/reactome/metabolism_hs_react_ids.txt
    """

}


/*
Get a list of pathway IDs for all the reactions rooted at 
(i.e. node ID ${params.ReactomeImmuneSystem_id}) and related to
Homo sapiens, which also contain the word "ignal"
*/
process immune_system_react_ids {

    publishDir "${out_dir}", 
        pattern: "databases/reactome/hs_imm_react_ids.txt",
        mode: 'copy'

    input:
        path 'input/ReactomePathways.txt'
        path 'input/ReactomePathwaysRelation.txt'
        path 'input/Ensembl2Reactome.txt'

    output:
        path 'databases/reactome/hs_imm_react_ids.txt'

    script:
    """
    mkdir -p databases/reactome/
    
    visitReactomeGraph.py \
        input/ReactomePathwaysRelation.txt \
        ${params.ReactomeImmuneSystem_id} \
        | sort | uniq \
        > sig_react_ids.txt
    
    cat input/Ensembl2Reactome.txt \
        | awk -F "\\t" '\$6~/Homo sapiens/ && \$4~/ignal/{print \$2}' \
        | sort | uniq \
        > hs_path_ids.txt
    
    cat hs_path_ids.txt | grep -f sig_react_ids.txt \
        > databases/reactome/hs_imm_react_ids.txt
    """

}


/*
Generate a list of Ensembl ids for each reaction
*/
process ensembl_react_lists {

    input:
        path 'input/hs_sig_react_ids.txt'
        path 'input/Ensembl2Reactome.txt'

    output:
        path 'Ensembl_ids/*txt'

    script:
    """
    mkdir -p Ensembl_ids

    cat input/Ensembl2Reactome.txt | grep -f input/hs_sig_react_ids.txt \
        | awk -F '\t' '\$6~/Homo sapiens/{print \$1"\t"\$2}' \
        > Ensembl_react_relation.tsv

    for react_id in \$(cat input/hs_sig_react_ids.txt); do \
        cat Ensembl_react_relation.tsv | grep "\${react_id}" \
        | cut -f1 | sort | uniq > Ensembl_ids/\${react_id}.txt; \
    done
    """
}


/*
Generate a list of Ensembl ids for each metabolic reaction
*/
process metabolism_ensembl_react_lists {

    input:
        path 'input/metabolism_hs_react_ids.txt'
        path 'input/Ensembl2Reactome.txt'

    output:
        path 'Ensembl_ids/*txt'

    script:
    """
    mkdir -p Ensembl_ids

    cat input/Ensembl2Reactome.txt | grep -f input/metabolism_hs_react_ids.txt \
        | awk -F '\t' '\$6~/Homo sapiens/{print \$1"\t"\$2}' \
        > Ensembl_react_relation.tsv

    for react_id in \$(cat input/metabolism_hs_react_ids.txt); do \
        cat Ensembl_react_relation.tsv | grep "\${react_id}" \
        | cut -f1 | sort | uniq > Ensembl_ids/\${react_id}.txt; \
    done
    """
}


/*
=== NOT USED ===

Make dictionary of relevant reactome genes
Before querying the UniProt resource for translation, the ids are truncated
like: P11362-19 --> P11362

NOTE: there can be multiple ENSG IDs mapping to the same UniProt ID
      (see Q9NY61 on Ensembl as an example: 
      ENSG00000275700 and ENSG00000276072 - both corresponding to AATF)
      Apparently duplicated ENSG occurr when the gene was mapped on a non
      canonical chr, and then mapped on the correct location and integrated
      among the CCDS set (how to get only the CCDS ENSG ID?)

.META:
1   uniprot ID
2   ENSEMBL ID
*/
process uniprot2ensg {

    publishDir "${out_dir}",
                pattern: "databases/reactome/uniprot2ensg.tsv",
                mode: 'copy'

    input:
        path 'input/*.txt'

    output:
        path 'databases/reactome/uniprot2ensg.tsv'

    script:
    """
    mkdir -p databases/reactome

    cat input/*.txt | sed 's/-.*\$//g' | sort | uniq \
        | uniprot2ensg.py | sed '1d' \
        > databases/reactome/uniprot2ensg.tsv
    """

}

/*
Get the list of gene names (UniProt AC) found in the Signal Transduction
subgraph
*/
process get_signalling_genes {

    input:
        path 'input/*.txt'

    output:
        path 'genes.txt'

    script:
    """
    cat input/*.txt | sed 's/-.*\$//g' | sort | uniq \
        > genes.txt
    """

}


/*
Get the list of gene names (UniProt AC) found in the Metabolism
subgraph
*/
process get_metabolism_genes {

    input:
        path 'input/*.txt'

    output:
        path 'genes.txt'

    script:
    """
    cat input/*.txt | sed 's/-.*\$//g' | sort | uniq \
        > genes.txt
    """

}


/*
List positive training examples from all the gene pairs that
participate in the given Reactome ID

.META:
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0)
3   gene 1
4   gene 2
*/
process pos_examples {

    input:
        tuple val(react_id),
              path('input/list.txt')

    output:
        path "${react_id}.tsv"

    script:
    """
    positive_examples.py \
            input/list.txt \
            ${react_id} \
        | sort | uniq \
        > ${react_id}.tsv
    """

}

/*
remove redundancies in the positive edges examples
*/
process pos_ex_unique_edges {

    input:
        path 'input/*'

    output:
        env(n), emit: n
        path 'non_redundant_pos_ex.tsv', emit: edges

    script:
    """
    cat input/* | awk 'NF' | cut -f3,4 | sort | uniq | awk 'NF{print \$1"_"\$2}' > for
    cat input/* | awk 'NF' | cut -f3,4 | sort | uniq | awk 'NF{print \$2"_"\$1}' > rev
    cat for | grep -f <(cat rev) | tr '_' '\\t' > dup_by_rev
    cat input/* | awk 'NF' | cut -f3,4 | sort | uniq -d > dup_by_for
    cat dup_by_for dup_by_rev | sort | uniq > dup
    cat input/* | awk 'NF' | grep -v -f <(cat dup) > pos_ex_uniq_only
    cat dup | awk '{print "1\\tmultiple\\t"\$1"\\t"\$2}' > pos_ex_dup
    cat pos_ex_uniq_only pos_ex_dup | sort | uniq > non_redundant_pos_ex.tsv
    n=\$(cat non_redundant_pos_ex.tsv | wc -l | cut -d " " -f1)
    """

}

/*
For a given Reactome ID, pick negative training examples, i.e. gene pairs
that do not participate in the same reaction, where one gene is annotated to
the given Reactome ID, and the other is not

Examaples are restricted to Homo sapiens only
TODO: consider whether to restrict to Reactome signalling only
      (i.e. add grep -f input/hs_sig_react_ids.txt)

params.neg_examples_proportion_per_ReactomeID defines the amount of negative
examples for a given Reactome ID as the factor that has to be multiplied with
the amount of positive examples generate from the same ReactomeID
    e.g. params.neg_examples_proportion_per_ReactomeID = 1 implies the amount
         of negative examples being equal to the amount of positive examples
         params.neg_examples_proportion_per_ReactomeID = 2 implies the amount
         of negative examples being twice the amount of positive examples

NOTE: a negative example in one pathway might correspond to a positive example
      in another pathway, due to proteins partecipating in multiple pathways
      Those negative examples are dropped later in modules/training.nf

.META:
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0)
3   Uniprot ID of gene 1
4   Uniprot ID of gene 2
*/
process neg_examples {

    input:
        val N
        tuple val(react_id),
              path('input/list.txt')
        path 'input/uniprot2reactome.tsv'

    output:
        path 'examples.tsv'

    script:
    """
    cat input/uniprot2reactome.tsv \
        | awk -F '\t' '\$6~/Homo sapiens/' \
        | grep -v -f input/list.txt \
        | grep -v "${react_id}" \
        | cut -f1 | sort | uniq \
        > not_list.txt
    negative_examples.py \
        ${N} \
        input/list.txt \
        not_list.txt \
        > examples.tsv
    """

}

/*
generate the "clusters" file aggregating individual reactome pathways
to be used in gene set enrichment tests

.META:
1   reactome id
2   tab-separated gene names
*/
process list_terms {

    publishDir "${out_dir}",
                pattern: "databases/reactome/hs_sig_react_sets.tsv",
                mode: 'copy'

    input:
        path 'input/*.txt'

    output:
        tuple val('reactome'),
              file('databases/reactome/hs_sig_react_sets.tsv')

    script:
    """
    mkdir -p databases/reactome

    for file in \$(ls input/); do \
        name=\$(readlink -f "input/\$file"); \
        name="\${name##*/}"; name="\${name%.*}"; \
        genes=\$(cat input/\$file | tr '\\n' '\\t'); \
        echo -e "\${name}\\t\${genes}"; \
    done > databases/reactome/hs_sig_react_sets.tsv
    """

}


/*
generate the "clusters" file aggregating individual reactome pathways
to be used in gene set enrichment tests

.META:
1   reactome id
2   tab-separated gene names
*/
process metabolism_list_terms {

    publishDir "${out_dir}",
                pattern: "databases/reactome/metabolism_hs_react_sets.tsv",
                mode: 'copy'

    input:
        path 'input/*.txt'

    output:
        tuple val('reactome'),
              file('databases/reactome/metabolism_hs_react_sets.tsv')

    script:
    """
    mkdir -p databases/reactome

    for file in \$(ls input/); do \
        name=\$(readlink -f "input/\$file"); \
        name="\${name##*/}"; name="\${name%.*}"; \
        genes=\$(cat input/\$file | tr '\\n' '\\t'); \
        echo -e "\${name}\\t\${genes}"; \
    done > databases/reactome/metabolism_hs_react_sets.tsv
    """

}


/*
translate the "clusters" file obtained from process list_terms()
using a dictionary

translation goes from ENSG id to gene symbol

.META:
1   reactome id
2   tab-separated gene names (gene symbols)
*/
process translate_tab {

    publishDir "${out_dir}",
                pattern: "databases/reactome/hs_sig_react_sets_translated.tsv",
                mode: 'copy'

    input:
        tuple val('reactome'),
              file('input/tab.tsv')
        path 'input/dictionary.txt'

    output:
        tuple val('reactome'),
              file('databases/reactome/hs_sig_react_sets_translated.tsv')

    script:
    """
    mkdir -p databases/reactome
    
    translateTab.py \
        input/dictionary.txt \
        input/tab.tsv \
        reactome \
        gene_name \
        all \
        0 \
        > tab_translated.tsv

    paste \
        <(cat input/tab.tsv | cut -f1) \
        <(cat tab_translated.tsv) \
        > databases/reactome/hs_sig_react_sets_translated.tsv
    """

}

/*
plot density of set size for leaf nodes in 
Signal Transduction Reactome subgraph
*/
process sig_react_sets_stats {

    publishDir "${out_dir}",
                pattern: "databases/reactome/hs_sig_react_sets_stats.pdf",
                mode: 'copy'

    input:
        path 'input/hs_sig_react_sets_translated.tsv'

    output:
        path 'databases/reactome/hs_sig_react_sets_stats.pdf'

    script:
    """
    mkdir -p databases/reactome

    sets_stats.py \
        input/hs_sig_react_sets_translated.tsv \
        databases/reactome/hs_sig_react_sets_stats.pdf
    """    
}

/*
translate the tab separated genes in new line-separated sets for the
reactome in homo sapiens using a dictionary

translation goes from ENSG id to gene symbol

if a set does not contain any genes after the translation because the
translation does not exist, it will be discarded

.META:
1   reactome id
2   tab-separated gene names (gene symbols)
*/
process translate_tab_all {

    publishDir "${out_dir}",
                pattern: "databases/reactome/${sets_name}_translated.tsv",
                mode: 'copy'

    input:
        tuple val(sets_name),
              file('input/tab.tsv')
        path 'input/dictionary.txt'

    output:
        tuple val(sets_name),
              file("databases/reactome/${sets_name}_translated.tsv")

    script:
    """
    mkdir -p databases/reactome
    
    translateTab.py \
        input/dictionary.txt \
        input/tab.tsv \
        reactome \
        gene_name \
        all \
        0 \
        > tab_translated.tsv

    paste \
        <(cat input/tab.tsv | cut -f1) \
        <(cat tab_translated.tsv) \
        | awk 'NF>1' > databases/reactome/${sets_name}_translated.tsv
    """

}

/*
all gene sets for human Reactome

.META:
1   reactome id
2   tab-separated gene names
*/
process hs_react_sets {

    publishDir "${out_dir}",
                pattern: "databases/reactome/hs_react_sets.tsv",
                mode: 'copy'

    input:
        path 'input/ReactomePathways.txt'
        path 'input/Ensembl2Reactome.txt'

    output:
        tuple val('hs_reactome'),
              file('databases/reactome/hs_react_sets.tsv')

    script:
    """
    mkdir -p databases/reactome

    cat input/Ensembl2Reactome.txt \
        | grep "Homo sapiens" \
        > databases/reactome/hs_Ensembl2Reactome.txt

    cat input/ReactomePathways.txt \
        | grep "Homo sapiens" | cut -f1 | awk 'NF' \
        > databases/reactome/hs_react_ids.txt

    for react_id in \$(cat databases/reactome/hs_react_ids.txt); do \
        genes=\$(cat databases/reactome/hs_Ensembl2Reactome.txt \
            | grep \${react_id} | cut -f1 | tr '\\n' '\\t'); \
        echo -e "\${react_id}\t\${genes}"; \
    done \
    | awk 'NF>1' > databases/reactome/hs_react_sets.tsv
    """
}


/*
gene list from genes (ENSG) that are annotated in the human reactome
*/
process hs_gene_list {

    publishDir "${out_dir}",
                pattern: "databases/reactome/hs_react_genes.txt",
                mode: 'copy'

    input:
        tuple val(gene_set_name),
              file('input/hs_react_sets.tsv')
    
    output:
        path 'databases/reactome/hs_react_genes.txt'

    script:
    """
    mkdir -p databases/reactome

    cat input/hs_react_sets.tsv \
        | cut -f2- | tr '\\t' '\\n' | awk 'NF' \
        | sort | uniq > databases/reactome/hs_react_genes.txt
    """

}

/*
*/
process signal_transduction_reactome_leaves {

    publishDir "${out_dir}",
                pattern: "databases/reactome/signal_transduction_reactome_leaves.tsv",
                mode: 'copy'

    input:
        path 'input/*.tsv'
    
    output:
        tuple val('reactome_signalling'),
              file('databases/reactome/signal_transduction_reactome_leaves.tsv')

    script:
    """
    mkdir -p databases/reactome
    mkdir -p terms

    for file in \$(ls input/); do \
        linkpath=\$(readlink -f input/"\$file"); \
        filebasename=\$(basename \$linkpath | sed 's/.tsv//'); \
        cat \
            <(echo "\$filebasename") \
            <(cat \$linkpath) \
            | tr '\\n' '\\t' \
            | sed -r 's/\\t\$/\\n/' \
        > terms/\${filebasename}.tsv; \
    done

    cat terms/*.tsv \
        > databases/reactome/signal_transduction_reactome_leaves.tsv

    """

}