#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Download relative protein abundance from an Expression Atlas proteomics
baseline experiment, given an E-PROT identifier (e.g. "29")
Note: E-PROT-46 not found (https://www.ebi.ac.uk/gxa/experiments/E-PROT-46/Results)
*/
process download_eprot {

    publishDir "${out_dir}",
                pattern: "databases/eprot/experiments/e-prot-${ID}.tsv",
                mode: 'copy'

    input:
        val ID

    output:
        path "databases/eprot/experiments/e-prot-${ID}.tsv"

    script:
    """
    mkdir -p databases/eprot/experiments
    URL=\$(echo '${params.url_eprot_template}' | sed 's/XX/${ID}/')
    wget -O databases/eprot/experiments/e-prot-${ID}.tsv "\${URL}"
    """

}

/*
Generate protein abundance tables from E-PROT experiments

.META:
1   ENSG ID
2   <e-prot-id>__<sample_1> abundance
3   <e-prot-id>__<sample_2> abundance
... ...
*/
process parse_eprot {

    input:
        tuple val(id), file('input/exp.tsv')

    output:
        path "databases/eprot/parsed/${id}.tsv"

    script:
    """
    mkdir -p databases/eprot/parsed

    old_cols=\$(cat input/exp.tsv | sed -n '5p' | cut -f3- \
               | sed 's/ /_/g' | tr '\t' '\n')

    new_cols=\$(for col in \${old_cols}; do \
                    echo "${id}__\${col}"; \
                done | tr '\n' '\t' | sed -e 's/\t\$//')
    cat \
        <(echo -e "Gene_ID\t\${new_cols}") \
        <(cat input/exp.tsv | sed '1,5d' | cut -f1,3- ) \
        > databases/eprot/parsed/${id}.tsv
    """

}

/*
Merge the E-PROT experiments from PARSE_EPROT() into a single table

13329 ENSG (rows)
137 samples (columns)

Note: here the invalid gene id "55872" found in E-PROT-29 is replaced with
      "ENSG00000168078", according to
      https://github.com/biothings/mygene.info/issues/94
*/
process merge_eprot_exp {

    publishDir "${out_dir}",
                pattern: "databases/eprot/all_eprot.tsv",
                mode: 'copy'

    input:
        path 'input/e-prot-*.tsv'

    output:
        path 'databases/eprot/all_eprot.tsv', emit: table

    script:
    """
    mkdir -p databases/eprot
    
    for eprot_file in \$(ls input/); do \
        exp_id=\$(basename \${eprot_file} .tsv); \
        echo "\${exp_id}" >> databases/eprot/exp_ids.txt; \
        sed -i 's/^55872/ENSG00000168078/' input/\${eprot_file}; \
    done

    mergeEPROTexp.py \
        databases/eprot/exp_ids.txt \
        input \
        databases/eprot/reported_ensg_aboundance_count.txt \
        > databases/eprot/all_eprot.tsv
    """

}

/*
Get the list of gene names in the merged EPROT matrix to be translated
*/
process get_eprot_genes {

    publishDir "${out_dir}",
                pattern: "databases/eprot/genes.txt",
                mode: 'copy'

    input:
        path 'input/all_eprot.tsv'

    output:
        path 'databases/eprot/genes.txt'

    script:
    """
    mkdir -p databases/eprot

    cat input/all_eprot.tsv | cut -f1 | sed '1d' | awk 'NF' \
        | sort | uniq \
        > databases/eprot/genes.txt
    """

}


/*
Build the partial input vector from E-PROT data

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2

.META: eprot_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions) 
3   gene 1
4   gene 2
... eprot features
*/
process eprot_featvec {

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/eprot_table.tsv')

    output:
        path 'eprot_input_vector.tsv'

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi
    
    eprot_input_vector.py  \
        input/gene_pairs.tsv \
        input/eprot_table.tsv \
        > eprot_input_vector.tsv
    """ 

}