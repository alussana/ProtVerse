#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
download KEGG_2021_Human from EnrichR gene sets repository
<https://maayanlab.cloud/Enrichr/#libraries>
*/
process dl_kegg_2021_human {

    publishDir "${out_dir}",
               pattern: "databases/kegg/KEGG_2021_Human.tsv",
               mode: 'copy'

    output:
        tuple val('KEGG_2021_Human'),
              file('databases/kegg/KEGG_2021_Human.tsv')

    script:
    """
    mkdir -p databases/kegg
    
    wget -O databases/kegg/KEGG_2021_Human.tsv \
        '${params.url_kegg2021_enrichr}'
    """

}