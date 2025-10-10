#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
download gene expression across tissues from human protein atlas
<https://www.proteinatlas.org/search/has_protein_data_in:Tissue>
*/
process dl_hpa_tissue_rna {

    publishDir "${out_dir}",
                pattern: "databases/hpa/hpa_tissueRNA.tsv",
                mode: 'copy'

    output:
        path 'databases/hpa/hpa_tissueRNA.tsv'

    script:
    """
    mkdir -p databases/hpa

    wget -O databases/hpa/hpa_tissueRNA.tsv '${params.url_hpa_tissueRNA}'
    """

}