#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
Publish the ubiquitination data table to the output directory

.META: input/source.csv
 1  Protein                                     sp|A0A024RBG1|NUD4B_HUMAN
 2  Protein_pos                                 134
 3  PTM_residue                                 K
 4  Decoy_mod                                   0
 5  PTM_FLR_category                            Gold
 6  Sum_of_PSM_counts(5%FLR)                    61.0
 7  PXD006201_peptide_mod_pos                   VLQC[Carbamidomethyl]HK[Ubiq]PVHAEYLEK-6
 8  PXD006201_FLR                               0.0049405020033906
 9  PXD006201_BinomialScore                     0.9967026
10  PXD006201_peptidoform_PSMcount(5%FLR)       ...
11  PXD019854_peptide_mod_pos                   ...
12  PXD019854_FLR                               ...
13  PXD019854_BinomialScore
14  PXD019854_peptidoform_PSMcount(5%FLR)
15  PXD022367_peptide_mod_pos
16  PXD022367_FLR
17  PXD022367_BinomialScore
18  PXD022367_peptidoform_PSMcount(5%FLR)
19  PXD023218_peptide_mod_pos
20  PXD023218_FLR
21  PXD023218_BinomialScore
22  PXD023218_peptidoform_PSMcount(5%FLR)
23  PXD023889_peptide_mod_pos
24  PXD023889_FLR
25  PXD023889_BinomialScore
26  PXD023889_peptidoform_PSMcount(5%FLR)
27  PXD025890_peptide_mod_pos
28  PXD025890_FLR
29  PXD025890_BinomialScore
30  PXD025890_peptidoform_PSMcount(5%FLR)
31  PXD027328_peptide_mod_pos
32  PXD027328_FLR
33  PXD027328_BinomialScore
34  PXD027328_peptidoform_PSMcount(5%FLR)
35  PXD037009_peptide_mod_pos
36  PXD037009_FLR
37  PXD037009_BinomialScore
38  PXD037009_peptidoform_PSMcount(5%FLR)
39  PXD019692_peptide_mod_pos
40  PXD019692_FLR
41  PXD019692_BinomialScore
42  PXD019692_peptidoform_PSMcount(5%FLR)
43  PXD020909_peptide_mod_pos
44  PXD020909_FLR
45  PXD020909_BinomialScore
46  PXD020909_peptidoform_PSMcount(5%FLR)
47  PXD003936_peptide_mod_pos
48  PXD003936_FLR
49  PXD003936_BinomialScore
50  PXD003936_peptidoform_PSMcount(5%FLR)
*/
process get_ubiquitination {

    publishDir "${out_dir}",
                pattern: "databases/ubiquitination/source.tsv",
                mode: 'copy'

    input:
        path 'input/source.csv'

    output:
        path 'databases/ubiquitination/source.tsv'

    script:
    """
    mkdir -p databases/ubiquitination

    cat input/source.csv \\
        | sed 's/,/\\t/g' \\
        | awk 'NF' \\
        > databases/ubiquitination/source.tsv
    """
}


/*
Filter ubiquitination sites with FDR<0.01 ("Silver" and "Gold" categories)

.META: databases/ubiquitination/source.tsv
[...]
*/
process parse_ubiquitination {

    publishDir "${out_dir}",
                pattern: "databases/ubiquitination/table.tsv",
                mode: 'copy'

    input:
        path 'input/table.tsv'

    output:
        path 'databases/ubiquitination/table.tsv'

    script:
    """
    mkdir -p databases/ubiquitination

    cat \\
        <(echo -e "UniProtAC\\tgene_name\\t\$(cat input/table.tsv | sed -n '1p' | cut -f2-)") \\
        <(cat input/table.tsv | sed '1d' | sed 's/|/\\t/' | cut -f2- | sed 's/|/\\t/;s/_HUMAN//') \\
        > table.tsv

    parse_ubiquitination.py \
        table.tsv \
        > databases/ubiquitination/table.tsv
    """
}


/* 
list unique ids found in the ubiquitination table
*/
process get_ubiquitination_genes {

    publishDir "${out_dir}",
                pattern: "databases/ubiquitination/genes.txt",
                mode: 'copy'

    input:
        path 'input/table.tsv'

    output:
        path 'databases/ubiquitination/genes.txt'

    script:
    """
    mkdir -p databases/ubiquitination

    cat input/table.tsv | cut -f1 | sed '1d' | awk 'NF' \
        | sort | uniq \
        > databases/ubiquitination/genes.txt
    """
}


/*
Build the partial input vector from ubiquitination data

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

.META: ubiquitination_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2
... ubiquitination features
*/
process ubiquitnation_featvec {

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/all_logFC.tsv')

    output:
        path 'ubiquitination_input_vector.tsv'

    script:
    """
    cat input/all_logFC.tsv | sed 's/\\t/_/' > all_logFC.tsv

    ubiquitination_input_vector.py  \
        input/gene_pairs.tsv \
        all_logFC.tsv \
        > ubiquitination_input_vector.tsv
    """ 

}


/*
Build the partial input vector from ubiquitination data

NOTE: there is only one ubiquitination site per protein, the features are
      computed considering the measurements to be protein-specific rather
      than site-specific

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2

.META: ubiquitination_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2
... ubiquitination features
*/
process ubiquitination_featvec {

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/table.tsv')

    output:
        path 'ubiquitination_input_vector.tsv'

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi

    cat input/table.tsv | sed 's/_/\\t/' | cut -f1,3- > table.tsv

    ubiquitination_input_vector.py  \
        input/gene_pairs.tsv \
        table.tsv \
        > ubiquitination_input_vector.tsv
    """ 

}
