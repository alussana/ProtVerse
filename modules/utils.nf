#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
publish a file giving it an arbitrary path
*/
process publish {

    publishDir "${out_dir}",
                pattern: "${outputFile}",
                mode: 'copy'

    input:
        path 'input/file'
        val outputFile

    output:
        path "${outputFile}"

    script:
    """
    mkdir -p \$(dirname ${outputFile})
    
    cp -r input/file ${outputFile}
    """

}

/*
split a single input file into chunks of n rows,
each named with a unique identifier

typically, each chunk in the output channel is mapped to its id with
.flatten().map{ file -> tuple( file.baseName, file ) }
*/
process split {

    input:
        path 'input/file'
        val n

    output:
        path "output/*"

    script:
    """
    mkdir -p output
    
    split -l ${n} -a 16 -x input/file output/
    """

}

/*
concatenate input files removing empty lines
*/
process concatenate {

    input:
        path 'input/*'
    
    output:
        path 'output/file'
    
    script:
    """
    mkdir -p output

    cat input/* | awk 'NF' > output/file
    """

}

/*
Query the Uniprot database with a list of gene IDs in a given
${from_nomenclature} to get their translation in the given ${to_nomenclature}

See https://www.uniprot.org/help/api_idmapping for available gene
nomenclatures

.META: dict.txt
1   gene id in ${from_nomenclature}
2   gene id in ${to_nomenclature}
*/
process map_ids {

    input:
        path 'input/gene_list.txt'
        val from_nomenclature
        val to_nomenclature

    output:
        path 'dict.txt'

    script:
    """
    cat input/gene_list.txt \
        | idmapping.py \
            ${from_nomenclature} \
            ${to_nomenclature} \
            9606 \
        | sed '1d' | awk 'NF' \
        > dict.txt
    """

}

/*
Translate all the words specified in the first tab-separated column of
input/dictionary.txt with the corresponding word found in the second
column
*/
process translate {

    publishDir "${out_dir}",
                pattern: "filtered_data/${id}.tsv",
                mode: 'copy'

    input:
        tuple val(id),
              file('input/target.txt')
        path 'input/dictionary.txt'
        val dict_col

    output:
        path "filtered_data/${id}.tsv"

    script:
    """
    mkdir -p filtered_data/

    col=\$( \
        cat input/dictionary.txt | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${dict_col}"{print NR}' \
    )

    cat input/dictionary.txt | sed '1d' \
        | awk -v col="\$col" '{print \$"\\t"\$1}' \
        > dictionary.txt
    
    translator.awk \
        dictionary.txt \
        input/target.txt \
        > filtered_data/${id}.tsv
    """

}

/*
Translate all the words in the first field of input/target.txt 
specified in the first tab-separated column of input/dictionary.txt
with the corresponding word found in the second column
do not discard untranslated rows
*/
process translatepy {

    publishDir "${out_dir}",
                pattern: "filtered_data/${id}.tsv",
                mode: 'copy'

    input:
        tuple val(id),
              file('input/target.txt')
        path 'input/dictionary.txt'
        val dict_col

    output:
        path "filtered_data/${id}.tsv"

    script:
    """
    mkdir -p filtered_data/

    col=\$( \
        cat input/dictionary.txt | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${dict_col}"{print NR}' \
    )

    cat input/dictionary.txt | sed '1d' \
        | awk -v col="\$col" '{print \$col"\\t"\$1}' \
        | awk 'NF==2' | sort | uniq | grep -v -w "NA" > dictionary.txt
    
    translator.py \
        dictionary.txt \
        input/target.txt \
        1 \
        2 \
        1 \
        1 \
        > filtered_data/${id}.tsv
    """

}

/*
translate the gene names in the first column of a matrix
if multiple translations for a word are possible, all of them are reported
in the translated matrix by duplicating the row as many times as the
available translations
*/
process translate_expand_matrix {

    publishDir "${out_dir}",
                pattern: "filtered_data/${id}.tsv",
                mode: 'copy'

    input:
        tuple val(id),
              file('input/target.txt')
        path 'input/dictionary.txt'
        val dict_col

    output:
        path "filtered_data/${id}.tsv"

    script:
    """
    mkdir -p filtered_data/

    col=\$( \
        cat input/dictionary.txt | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${dict_col}"{print NR}' \
    )

    cat input/dictionary.txt | sed '1d' \
        | awk -v col="\$col" '{print \$col"\\t"\$1}' \
        | awk 'NF==2' | sort | uniq | grep -v -w "NA" > dictionary.txt
    
    translate_expand_matrix.py \
        dictionary.txt \
        input/target.txt \
        1 \
        2 \
        > filtered_data/${id}.tsv
    """
}

/*
translate the gene names in the first two column of a table
if multiple translations for a gene are possible, all of them are reported
in the translated table by duplicating the row as many times as the
possible combinations of translations
*/
process translate_expand_pairs {

    publishDir "${out_dir}",
                pattern: "filtered_data/${id}.tsv",
                mode: 'copy'

    input:
        tuple val(id),
              file('input/target.txt'),
              file('input/dictionary.txt'),
              val(dict_col)

    output:
        path "filtered_data/${id}*.tsv"

    script:
    """
    mkdir -p filtered_data/

    col=\$( \
        cat input/dictionary.txt | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${dict_col}"{print NR}' \
    )

    cat input/dictionary.txt | sed '1d' \
        | awk -v col="\$col" '{print \$col"\\t"\$1}' \
        | awk 'NF==2' | sort | uniq | grep -v -w "NA" > dictionary.txt

    uuid=\$(cat /proc/sys/kernel/random/uuid)
    
    translate_expand_pairs.py \
        dictionary.txt \
        input/target.txt \
        1 \
        2 \
        0 \
        > filtered_data/${id}_\${uuid}.tsv
    """
}

/*
Translate all the words in the first field of input/target.txt 
specified in the first tab-separated column of input/dictionary.txt
with the corresponding word found in the second column
do not discard untranslated rows
Do not publishDir on output
*/
process translatepy_no_pub {

    input:
        tuple val(id),
              file('input/target.txt'),
              file('input/dictionary.txt'),
              val(dict_col)

    output:
        path "filtered_data/${id}.tsv"

    script:
    """
    mkdir -p filtered_data/

    col=\$( \
        cat input/dictionary.txt | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${dict_col}"{print NR}' \
    )

    cat input/dictionary.txt | sed '1d' \
        | awk -v col="\$col" '{print \$col"\\t"\$1}' \
        | awk 'NF==2' | sort | uniq | grep -v -w "NA" > dictionary.txt
    
    translator.py \
        dictionary.txt \
        input/target.txt \
        1 \
        2 \
        1 \
        1 \
        > filtered_data/${id}.tsv
    """

}

/*
translate the gene names in a list
if multiple translations for a gene are possible, all of them are reported
in the translated list by duplicating the row as many times as the
possible translations
*/
process translate_expand_list {

    input:
        tuple val(id),
              file('input/target.txt'),
              file('input/dictionary.txt'),
              val(dict_col)

    output:
        path "filtered_data/${id}.tsv"

    script:
    """
    mkdir -p filtered_data/

    col=\$( \
        cat input/dictionary.txt | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${dict_col}"{print NR}' \
    )

    cat input/dictionary.txt | sed '1d' \
        | awk -v col="\$col" '{print \$col"\\t"\$1}' \
        | awk 'NF==2' | sort | uniq | grep -v -w "NA" > dictionary.txt
    
    translate_expand_list.py \
        dictionary.txt \
        input/target.txt \
        1 \
        2 \
        > filtered_data/${id}.tsv
    """

}

/*
Translate all the words specified in the first tab-separated column of
input/dictionary.txt with the corresponding word found in the second
column
Do not publishDir on output
*/
process translate_no_pub {

    input:
        tuple val(id),
              file('input/target.txt'),
              file('input/dictionary.txt'),
              val(dict_col)

    output:
        path "filtered_data/${id}.tsv"

    script:
    """
    mkdir -p filtered_data/

    col=\$( \
        cat input/dictionary.txt | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${dict_col}"{print NR}' \
    )

    cat input/dictionary.txt | sed '1d' \
        | awk -v col="\$col" '{print \$col"\\t"\$1}' \
        | awk 'NF==2' | grep -v -w "NA" > dictionary.txt
    
    translator.awk \
        dictionary.txt \
        input/target.txt \
        > filtered_data/${id}.tsv
    """

}

/*
Translate all the words specified in the first tab-separated column of
input/dictionary.txt with the corresponding word found in the second
column
Do not publishDir on output
*/
process translate_no_pub_rand_id {

    input:
        tuple val(id),
              file('input/target.txt'),
              file('input/dictionary.txt'),
              val(dict_col)

    output:
        path "filtered_data/${id}_*.tsv"

    script:
    """
    mkdir -p filtered_data/

    col=\$( \
        cat input/dictionary.txt | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${dict_col}"{print NR}' \
    )

    cat input/dictionary.txt | sed '1d' \
        | awk -v col="\$col" '{print \$col"\t"\$1}' \
        | awk 'NF==2' | sort | uniq | grep -v -w "NA" > dictionary.txt

    uuid=\$(cat /proc/sys/kernel/random/uuid)
    
    translator.awk \
        dictionary.txt \
        input/target.txt \
        | grep -v -w <(cat dictionary.txt | cut -f1 | sort | uniq) \
        > filtered_data/${id}_\${uuid}.tsv
    """

}

/*
Translate all the words in the first two fields of input/target.txt 
specified in the first tab-separated column of input/dictionary.txt
with the corresponding word found in the second column
Do not keep lines where both the first two fields had a translation
Do not publishDir on output
*/
process translatepy_no_pub_rand_id {

    input:
        tuple val(id),
              file('input/target.txt'),
              file('input/dictionary.txt'),
              val(dict_col)

    output:
        path "filtered_data/${id}_*.tsv"

    script:
    """
    mkdir -p filtered_data/

    col=\$( \
        cat input/dictionary.txt | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${dict_col}"{print NR}' \
    )

    cat input/dictionary.txt | sed '1d' \
        | awk -v col="\$col" '{print \$col"\t"\$1}' \
        | awk 'NF==2' | sort | uniq | grep -v -w "NA" > dictionary.txt

    uuid=\$(cat /proc/sys/kernel/random/uuid)
    
    translator.py \
        dictionary.txt \
        input/target.txt \
        1 \
        2 \
        1,2 \
        0 \
        > filtered_data/${id}_\${uuid}.tsv
    """

}

/*
Filter the rows of a data matrix (input/matrix.tsv) based on the indeces,
keeping only those found in a given gene id dictionary (input/gene_dict.tsv)
at the specified column col
*/
process filter_data_matrix_rows {

    input:
        path 'input/matrix.tsv'
        path 'input/gene_dict.tsv'
        val col
    
    output:
        path 'filtered_matrix.tsv'

    script:
    """
    col_n=\$( \
        cat input/gene_dict.tsv | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${col}"{print NR}' \
    )

    cat input/matrix.tsv \
        | grep -w -f <( cat input/gene_dict.tsv | sed '1d' \
                        | cut -f\${col_n} | grep -w -v "NA" \
                        | sort | uniq ) \
        > filtered_matrix.tsv
    """

}

/*
Filter the rows of a data matrix (input/matrix.tsv) based on the indeces,
keeping only those found in a given gene id dictionary (input/gene_dict.tsv)
at the specified column col
For id-mapped input

NOTE: for reasons not understood processes that return no gene ids exit with
      a non-0 status, even though reproducing the command in the nf work dir
      does not. For this reason a 0 exit status is forced in this script

*/
process filter_data_matrix_rows_w_id {

    input:
        tuple val(id),
              file('input/matrix.tsv'),
              file('input/gene_dict.tsv'),
              val(col)
    
    output:
        tuple val(id),
              file("${id}.tsv")              

    script:
    """
    col_n=\$( \
        cat input/gene_dict.tsv | sed -n '1p' \
        | tr '\\t' '\\n' | awk '\$0=="${col}"{print NR}' \
    )

    cat input/matrix.tsv \
        | grep -w -f <( cat input/gene_dict.tsv | sed '1d' \
                        | cut -f\${col_n} | grep -w -v "NA" \
                        | sort | uniq ) \
        > ${id}.tsv \
    #|| true
    """

}

/*
Concatenate input text files removing empty lines
*/
process concat {

    input:
        path 'input/*'

    output:
        path 'concat.txt'

    script:
    """
    cat input/* | awk 'NF' > concat.txt
    """

}

/*
concatenate input text files removing empty lines, then split in smaller files
of n lines each
*/
process cat_and_split {

    input:
        path 'input/*'
        val n

    output:
        path 'cat_and_split/*.tsv'

    script:
    """
    mkdir -p cat_and_split

    cat input/* | awk 'NF' > cat.txt
    
    split -l ${n} -a 16 -x cat.txt cat_and_split/ --additional-suffix .tsv
    """

}

/*
Concatenate input text files removing empty lines
and publishDir using id
Add header as first line of the output
*/
process concat_w_id {

    debug false

    publishDir "${out_dir}",
                pattern: "filtered_data/${id}",
                mode: 'copy'

    input:
        val id
        path 'input/*'
        val header

    output:
        path "filtered_data/${id}"

    script:
    """
    mkdir -p filtered_data/
    cat \
        <(echo -e "${header}") \
        <(cat input/*) \
        | awk 'NF' \
        > filtered_data/${id}
    """

}


/*
3-fields tab-delimited file

provides mapping in this order:

IDa --->  UniProt AC ---> IDb

all ids of nomenclature IDa are mapped to the corresponding UniProt ACs,
and the UniProt ACs are then mapped to the corresponding ids of nomenclature
IDb. For example, all ENSP ids are mapped to their UniProt AC. Then, each of
the UniProt AC identifiers is mapped to HGNC ids.

.META:
1. IDa
2. UniProt AC referred to IDa
3. IDb referred to UniProt AC
*/
process IDa2uniprot2IDb {

    publishDir "${out_dir}",
                pattern: "databases/uniprot/${IDa}2uniprot2${IDb}.tsv",
                mode: 'copy'

    input:
        path 'input/mapping.tsv'
        val IDa
        val IDb

    output:
        path "databases/uniprot/${IDa}2uniprot2${IDb}.tsv"

    script:
    """
    mkdir -p databases/uniprot

    grep -w -f <(echo -e "${IDa}\\n${IDb}") \\
        < input/mapping.tsv \\
        | gzip > mapping.tsv.gz

    IDa2uniprot2IDb.py \\
        mapping.tsv.gz \\
        ${IDa} \\
        ${IDb} \\
        > databases/uniprot/${IDa}2uniprot2${IDb}.tsv
    """

}


/*
3-fields tab-delimited file

provides mapping in this order:

IDa --->  UniProt AC ---> IDb

all ids of nomenclature IDa are mapped to the corresponding UniProt ACs,
and the UniProt ACs are then mapped to the corresponding ids of nomenclature
IDb. For example, all ENSP ids are mapped to their UniProt AC. Then, each of
the UniProt AC identifiers is mapped to HGNC ids.

.META:
1. IDa
2. UniProt AC referred to IDa
3. IDb referred to UniProt AC
*/
process IDa2uniprot2IDb_vec {

    publishDir "${out_dir}",
                pattern: "databases/uniprot/${IDa}2uniprot2${IDb}.tsv",
                mode: 'copy'

    input:
        path 'input/mapping.tsv'
        val IDa
        val IDb

    output:
        path "databases/uniprot/${IDa}2uniprot2${IDb}.tsv"

    script:
    """
    mkdir -p databases/uniprot

    grep -w -f <(echo -e "${IDa}\\n${IDb}") \\
        < input/mapping.tsv \\
        | gzip > mapping.tsv.gz

    IDa2uniprot2IDb_vec.py \\
        mapping.tsv.gz \\
        ${IDa} \\
        ${IDb} \\
        > databases/uniprot/${IDa}2uniprot2${IDb}.tsv
    """

}
