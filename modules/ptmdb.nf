/*
Extract individual phosphoproteomics experiments results from an exported
portion of ptmdb. Each individual experiment is a text file

ptmdb is an SQL database of post-translation modifications
<https://github.com/evocellnet/ptmdb>

The data used here have been exported from ptmdb by Girolamo Giudice
<ggiudice@ebi.ac.uk>

NOTE: files contain many duplicated lines

.META: *txt
1   FS (Functional Score)
2   uniprot (ID)
3   lfc (LogFC)
4   pos (phosphosite position)
*/
process open_ptmdb {

    publishDir "${out_dir}", 
                pattern: "databases/ptmdb/experiments/*.txt",
                mode: 'copy'

    input:
        path 'ptmdb.zip'

    output:
        path "databases/ptmdb/experiments/*.txt"

    script:
    """
    mkdir -p databases/ptmdb/experiments
    unzip ptmdb.zip
    for id in \$(ls all_db2); do \
        for exp in \$(ls all_db2/\${id}); do \
            DEST=\$(echo "databases/ptmdb/experiments/\${id}__\${exp}" | tr -d '()'); \
            cp all_db2/\${id}/\${exp} \${DEST}; \
         done; \
    done
    """

}

/*
Generate logFC tables from ptmdb experiments

NOTE: some phosphosites might have been duplicated during parsing from PTMDB
      (https://github.com/evocellnet/ptmdb). This may be at least partially
      due to some of them coming from different data sources, which may or may
      not have an associated funtional score

ISSUE: even when not considering phosphosites with no functional score ("NA"),
       there will be some duplicated phosphosites having different lfc!

.META:
1   uniprot_pos
2   <exp_name>_logFC
*/
process parse_ptmdb {

    publishDir "${out_dir}", 
                pattern: "databases/ptmdb/parsed/${id}.tsv",
                mode: 'copy'

    input:
        tuple val(id), file('input/exp.tsv')

    output:
        path "databases/ptmdb/parsed/${id}.tsv"

    script:
    """
    mkdir -p databases/ptmdb/parsed
    cat \
        <(echo -e "uniprot_pos\t${id}_logFC") \
        <(cat input/exp.tsv \
            | sed '1d' \
            | awk '{print \$2"_"\$4"\t"\$3}' \
            | sort | uniq) \
        > databases/ptmdb/parsed/${id}.tsv
    """

}

/*
Merge the phosphoproteomics experiment from PARSE_PTMDB() in a single table

This includes data for 72635 phosphosites belonging to 7784 unique proteins 

ISSUE: because of what commented in PARSE_PTMDB, here mergePTMDBexp.py keeps
       only the first occurrence of a rows that have duplicated indeces in
       order to merge the experiments in a single dataframe

.META:
1   <uniprot id>
2   <position>
3   logFC in experiment 1
4   logFC in experiment 2
... ...
*/
process merge_ptmdb_exp {

    memory '8G'

    publishDir "${out_dir}",
                pattern: 'databases/ptmdb/all_logFC.tsv',
                mode: 'copy'

    input:
        path 'input/*.tsv'

    output:
        path 'databases/ptmdb/all_logFC.tsv'

    script:
    """
    mkdir -p databases/ptmdb

    for exp_file in \$(ls input/); do \
        exp_id=\$(basename \${exp_file} .tsv); \
        echo "\${exp_id}" >> databases/ptmdb/exp_ids.txt; \
    done

    mergePTMDBexp.py \
        databases/ptmdb/exp_ids.txt \
        input \
        | sed 's/_/\t/' > databases/ptmdb/all_logFC.tsv
    """

}

/*
Get the list of gene names in the ptmdb to be translated
*/
process get_ptmdb_genes {

    publishDir "${out_dir}",
                pattern: 'databases/ptmdb/genes.txt',
                mode: 'copy'

    input:
        path 'input/ptmdb_matrix.tsv'

    output:
        path 'databases/ptmdb/genes.txt'

    script:
    """
    mkdir -p databases/ptmdb

    cat input/ptmdb_matrix.tsv \
        | cut -f1 | sed '1d' | cut -f1 -d "_" | awk 'NF' \
        | sort | uniq \
        > databases/ptmdb/genes.txt
    """

}

/*
Build the partial input vector from PTMDB data

NOTE: for two given proteins a and b, corresponding phospho data are picked
      from the pair of phosphosites p_a and p_b -one for each protein-
      according to the following logic:
      - maximize pearson(p_a, p_b) ** 2
      - if pearson(p_a, p_b)==NaN for any (p_a, p_b), then choose (p_a, p_b)
        as the phosphosites detected more frequently

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2

.META: ptmdb_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2
... ptmdb features
*/
process ptmdb_featvec {

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/all_logFC.tsv')

    output:
        path 'ptmdb_input_vector.tsv'

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi
    
    cat input/all_logFC.tsv | sed 's/\\t/_/' > all_logFC.tsv

    ptmdb_input_vector.py  \
        input/gene_pairs.tsv \
        all_logFC.tsv \
        > ptmdb_input_vector.tsv
    """ 

}