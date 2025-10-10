#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Use the UniProt API to retrieve domain information for the proteins of the
human reference proteome

Note: not all entries have at least one gene name

.META:
1   entry                           id
2   length                          length
3   topological domain              ft_topo_dom
4   transmembrane                   ft_transmem
5   InterProt cross-reference       xref_interpro       colon-separated
6   gene names                      gene_names          space-separated, primary name is first, synonyms follow
*/
process dl_ref_proteome {

    publishDir "${out_dir}",
            pattern: 'databases/uniprot/hs_proteome.tsv.gz',
            mode: 'copy'

    output:
        path 'databases/uniprot/hs_proteome.tsv.gz'

    shell:
    """
    mkdir -p databases/uniprot

    wget -O databases/uniprot/hs_proteome.tsv.gz '${params.uniprot_query}'
    """

}

/*
.META:
1. UniProtKB-AC 
2. ID_type 
3. ID
*/
process dl_human_ids {

    publishDir "${out_dir}",
                pattern: 'databases/uniprot/HUMAN_9606_idmapping.dat.gz',
                mode: 'copy'

    output:
        path 'databases/uniprot/HUMAN_9606_idmapping.dat.gz'

    shell:
    """
    mkdir -p databases/uniprot

    wget -P databases/uniprot '${params.url_human_proteome}'
    """

}

/*
filter the UniProt id mapping information to only incude mapping to:

Gene_Name
GeneID
Ensembl
Gene_Synonym

NOTE: Ensembl gets truncated e.g. ENSG00000281151.2 --> ENSG00000281151
*/
process filter_idmapping {

    publishDir "${out_dir}",
                pattern: 'databases/uniprot/idmapping.tsv',
                mode: 'copy'

    input:
        path 'input/HUMAN_9606_idmapping.dat.gz'

    output:
        path 'databases/uniprot/idmapping.tsv'

    script:
    """
    mkdir -p databases/uniprot

    echo -e "Gene_Name\\nGeneID\\nEnsembl\\nGene_Synonym" \
        > idtypes.txt

    zcat input/HUMAN_9606_idmapping.dat.gz \
        | grep -w -f idtypes.txt \
        | sed -r 's/(Ensembl\\tENSG[0-9]+)\\.[0-9]+\$/\\1/' \
        | sed -r 's/(UniProtKB-ID\\t[0-9A-Z]+)_HUMAN\$/\\1/' \
        > databases/uniprot/idmapping.tsv
    """

}

/*
define list of reference primary gene names for the human reference proteome
*/
process get_proteome_genes {

    publishDir "${out_dir}",
                pattern: 'databases/uniprot/proteome_ref_genes.txt',
                mode: 'copy'

    input:
        path 'input/hs_proteome.tsv.gz'

    output:
        path 'databases/uniprot/proteome_ref_genes.txt'

    script:
    """
    mkdir -p databases/uniprot
    
    zcat input/hs_proteome.tsv.gz \
        | cut -f6 | sed '1d' | awk 'NF' | cut -d ' ' -f1 \
        | sort | uniq \
        > databases/uniprot/proteome_ref_genes.txt
    """
    
}

/*
assign to each of the reference proteome gene names
a corresponding translation in each of the datasets

then modify reference gene names to remove underscores
e.g. HERVK_113 --> HERVK113
*/
process translate_ids {
    
    cpus "${params.starmap_n_proc}"

    publishDir "${out_dir}",
                pattern: 'databases/uniprot/id_translations.tsv',
                mode: 'copy'

    input:
        path 'input/ref_genes.txt'
        path 'input/idmapping.tsv'
        path 'input/signalling_genes.txt'
        path 'input/metabolism_genes.txt'
        path 'input/gtex_genes.txt'
        path 'input/eprot_genes.txt'
        path 'input/orthogroups_genes.txt'
        path 'input/depmap_genes.txt'
        path 'input/ubiquitination_genes.txt'
        path 'input/humap3_genes.txt'
        path 'input/ptmdb_genes.txt'
        path 'input/proteomehd_genes.txt'
        path 'input/mitchell2023_genes.txt'
        path 'input/reactome_genes.txt'

    output:
        path 'databases/uniprot/id_translations.tsv'

    script:
    """
    mkdir -p databases/uniprot

    create_gene_dict.py \
        input/ref_genes.txt \
        input/idmapping.tsv \
        input/signalling_genes.txt \
        input/gtex_genes.txt \
        input/eprot_genes.txt \
        input/orthogroups_genes.txt \
        input/depmap_genes.txt \
        input/ubiquitination_genes.txt \
        input/humap3_genes.txt \
        input/ptmdb_genes.txt \
        input/proteomehd_genes.txt \
        input/mitchell2023_genes.txt \
        input/metabolism_genes.txt \
        input/reactome_genes.txt \
        ${params.starmap_n_proc} \
        > databases/uniprot/id_translations_w_underscores.tsv

    cat databases/uniprot/id_translations_w_underscores.tsv \
        | sed '1d' | cut -f1 | awk 'NF' \
        | sed 's/_//' \
        > databases/uniprot/id_translations_index.txt

    cat databases/uniprot/id_translations_w_underscores.tsv \
        | sed -n '1p' \
        > databases/uniprot/id_translations_header.tsv

    cat databases/uniprot/id_translations_w_underscores.tsv \
        | sed '1d' | cut -f2- | awk 'NF' \
        > databases/uniprot/id_translations_values.tsv

    cat \
        <(cat databases/uniprot/id_translations_header.tsv) \
        <(cat \
            <(paste databases/uniprot/id_translations_index.txt \
                    databases/uniprot/id_translations_values.tsv \
            ) \
        ) \
        > databases/uniprot/id_translations.tsv
    """

}

/*
the dictionary of tractable genes is here defined including only
genes that can be mapped to all the data sources
*/
process filter_id_dict {

    publishDir "${out_dir}",
                pattern: 'filtered_data/id_dict.tsv',
                mode: 'copy'

    input:
        path 'input/id_dict.tsv'

    output:
        path 'filtered_data/id_dict.tsv'

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    from pathlib import Path
    
    Path('filtered_data').mkdir(parents=True, exist_ok=True)

    gene_dict_file = 'input/id_dict.tsv'

    gene_dict = pd.read_csv(gene_dict_file, sep='\t', index_col=0, dtype=str)

    neglected_colnames = ['reactome', 'humap']

    cols = [i for i in gene_dict.columns if i not in neglected_colnames]
    
    tractable_genes = list(
        gene_dict.loc[:, cols].dropna(axis=0).index
    )
    
    tractable_gene_dict = gene_dict.loc[tractable_genes, :]
    
    tractable_gene_dict.to_csv(
        'filtered_data/id_dict.tsv', sep='\\t',
        header=True, index=True, na_rep='NA'
    )
    """

}

/*
== NOT USED ==

3-fields tab-delimited file

provides mapping in this order:

gene name  --->  UniProt AC ---> refseq ID

all gene names are mapped to the corresponding UniProt ACs,
and the UniProt ACs are then mapped to the corresponding refseq gene IDs

.META:
1. gene name 
2. UniProt AC referred to gene name
3. refseq ID referred to UniProt AC
*/
process name2uniprot2refseq {

    conda 'envs/idmapping_pandas.yml'

    publishDir "${out_dir}",
                pattern: 'uniprot/name2uniprot2refseq.tsv',
                mode: 'copy'

    input:
        path 'input/mapping.tsv.gz'

    output:
        path 'uniprot/name2uniprot2refseq.tsv'

    script:
    """
    mkdir -p uniprot

    name2ac2refseq.py input/mapping.tsv.gz \
        > uniprot/name2uniprot2refseq.tsv
    """

}

/*
== NOT USED ==

get sequence length from human reference proteome's proteins

only for entries having at least one gene name

.META:
1   gene names              space-separated, primary name is first, synonyms follow
2   length
*/
process length {

    publishDir "${out_dir}",
            pattern: 'databases/uniprot/length.tsv',
            mode: 'copy'

    input:
        path 'input/hs_proteome.tsv.gz'

    output:
        path 'databases/uniprot/length.tsv'

    shell:
    """
    mkdir -p databases/uniprot

    zcat input/hs_proteome.tsv.gz \
        | awk -F "\t" '\$6!=""{print \$6"\t"\$2}' \
        > databases/uniprot/length.tsv
    """

}

/*
== NOT USED ==
== NOT COMPLETE == 

get topological domain information from human reference proteome's proteins

only for entries having at least one gene name

.META:
1   gene names              space-separated, primary name is first, synonyms follow
2   topological domains     colon-separated
*/
process topo_domain {

    publishDir "${out_dir}",
            pattern: 'databases/uniprot/topo_domain.tsv.gz',
            mode: 'copy'

    input:
        path 'input/hs_proteome.tsv.gz'

    output:
        path 'databases/uniprot/topo_domain.tsv.gz'

    shell:
    """
    mkdir -p databases/uniprot

    zcat input/hs_proteome.tsv.gz \
        | awk -F "\\t" '\$6!=""{print \$6"\\t"\$3}' \
        | sed -r 's/[; ]*TOPO_DOM [\\.0-9]+; \\/note="([a-zA-Z ]*)"; \\/evidence="[a-zA-Z:0-9_\\|, ]+"/\\1;/g' \
        | sed 's/;\$//' \
        | gzip > databases/uniprot/topo_domain.tsv.gz
    """

}

/*
== NOT USED ==
== NOT COMPLETE ==

get transmembrane domain information from human reference proteome's proteins

only for entries having at least one gene name

.META:
1   gene names              space-separated, primary name is first, synonyms follow
2   transmembrane           colon-separated
*/
process transmem {

    publishDir "${out_dir}",
            pattern: 'databases/uniprot/transmem.tsv.gz',
            mode: 'copy'

    input:
        path 'input/hs_proteome.tsv.gz'

    output:
        path 'databases/uniprot/transmem.tsv.gz'

    shell:
    """
    mkdir -p databases/uniprot

    zcat input/hs_proteome.tsv.gz \
        | awk -F "\\t" '\$6!=""{print \$6"\\t"\$4}' \
        | sed -r 's/[; ]*TRANSMEM [0-9]+\\.\\.[0-9]+; \\/note="([a-zA-Z ]*)"; \\/evidence="[a-zA-Z:0-9|]+"/\\1;/g' \
        | sed 's/;\$//' \
        | gzip > databases/uniprot/transmem.tsv.gz
    """

}

/*
== NOT USED ==

TODO: add process to filter for interpro domains only (not families)

get interpro domain annotation from human reference proteome's proteins

only for entries having at least one gene name

.META:
1   gene names              space-separated, primary name is first, synonyms follow
2   interpro                colon-separated
*/
process interpro {

    publishDir "${out_dir}",
            pattern: 'databases/uniprot/interpro.tsv.gz',
            mode: 'copy'

    input:
        path 'input/hs_proteome.tsv.gz'

    output:
        path 'databases/uniprot/interpro.tsv.gz'

    shell:
    """
    mkdir -p databases/uniprot

    zcat input/hs_proteome.tsv.gz \
        | awk -F "\t" '\$6!=""{print \$6"\t"\$5}' \
        | sed 's/;\$//' \
        | gzip > databases/uniprot/interpro.tsv.gz
    """

}


/*
filter the UniProt id mapping information to only incude mapping to:

Gene_Name
Gene_Synonym
*/
process uniprot2gene_name_and_synonym {

    publishDir "${out_dir}",
                pattern: 'databases/uniprot/uniprot2gene_synonym.tsv',
                mode: 'copy'

    input:
        path 'input/HUMAN_9606_idmapping.dat.gz'

    output:
        path 'databases/uniprot/uniprot2gene_synonym.tsv'

    script:
    """
    mkdir -p databases/uniprot

    echo -e "Gene_Name\\nGene_Synonym" \
        > idtypes.txt

    zcat input/HUMAN_9606_idmapping.dat.gz \
        | grep -w -f idtypes.txt \
        | sed -r 's/(UniProtKB-ID\\t[0-9A-Z]+)_HUMAN\$/\\1/' \
        > databases/uniprot/uniprot2gene_synonym.tsv
    """
}


/*
translate all the words in input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the third column
the second column of input/dict.tsv specifies whether the translation is to a
"Gene_Name" or a "Gene_Synonym"
UniProt AC ids are preferentially translated to Gene_Name; if a gene name is
not found in the dictionary, then we look for a Gene_Synonym 
keep untranslated rows
*/
process translate_ac_to_gene_name {

    input:
        tuple val(id), file('input/file.tsv')
        path 'input/dict.tsv'

    output:
        tuple val(id), file('translated_file.tsv')
        
    script:
    """
    translator_ac_2_gene_name_or_synonym.py \
        input/dict.tsv \
        input/file.tsv \
        1 \
        3 \
        all \
        1 \
        2 \
        > translated_file.tsv
    """
}


/*
Visualize coverage of reference proteome genes in the different datasets
*/
process viz_ref_genes_coverage {

    publishDir "${out_dir}",
                pattern: 'databases/uniprot/*pdf',
                mode: 'copy'

    input:
        path "input/ref_genes_id_mapping.tsv"

    output:
        path "databases/uniprot/*.pdf"
        
    script:
    """
    mkdir -p databases/uniprot
    
    viz_ref_genes_coverage.py \
        input/ref_genes_id_mapping.tsv \
        databases/uniprot/missing_cumsum.pdf \
        databases/uniprot/missing_bar.pdf \
        databases/uniprot/missing_matrix.pdf \
        databases/uniprot/missing_heatmap.pdf \
        databases/uniprot/missing_dendrogram.pdf
    """
}