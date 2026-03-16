#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Get gene dependencies fold changes computed in
<https://doi.org/10.1038/s41467-021-21898-7>

908 cell lines (columns) (model_id from <https://cellmodelpassports.sanger.ac.uk/>)
17486 genes (rows) (HUGO gene symbols)
*/
process download_essentiality_matrices {

    publishDir "${out_dir}",
                pattern: "databases/depMap/*txt",
                mode: 'copy'

    output:
        path 'databases/depMap/CERES_FC.txt', 
            emit: CERES_FC
        path 'databases/depMap/CRISPRcleanR_FC.txt',
            emit: CRISPRcleanR_FC

    script:
    """
    mkdir -p databases/depMap/
    wget --user-agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36" \\
         --referer="https://example.com/" \\
         --content-on-error --content-disposition -O databases/depMap/26705222.zip \\
         "${params.url_depMap}"

    unzip databases/depMap/26705222.zip

    mv integrated_Sanger_Broad_essentiality_matrices_20201201/CERES_FC.txt \\
        databases/depMap/CERES_FC.txt

    mv integrated_Sanger_Broad_essentiality_matrices_20201201/CRISPRcleanR_FC.txt \\
        databases/depMap/CRISPRcleanR_FC.txt

    rm -fr integrated_Sanger_Broad_essentiality_matrices_20201201
    """

}

/*
Get the list of gene names in the essentiality matrices to be translated
*/
process get_depmap_genes {

    publishDir "${out_dir}",
                pattern: "databases/depMap/genes.txt",
                mode: 'copy'

    input:
        path 'input/gene_dependencies.txt'

    output:
        path 'databases/depMap/genes.txt'

    script:
    """
    mkdir -p databases/depMap
    
    cat input/gene_dependencies.txt | cut -f1 | awk 'NF' \
        | sort | uniq \
        > databases/depMap/genes.txt
    """

}

/*
Build the partial input vector from gene dependency data

.META: gene_pairs.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2

.META: dependency_input_vector.tsv
1   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
2   Reactome ID (or NA if field 1 == 0) ("multiple" if the two genes are annotated in multiple reactions)
3   gene 1
4   gene 2
... depmap features
*/
process dependency_featvec {

    input:
        tuple file('input/gene_pairs.tsv'),
              file('input/dependency_table.tsv')

    output:
        path 'dependency_input_vector.tsv'

    script:
    """
    if [ -z "\${PYTHONPATH:-}" ]; then \\
        export PYTHONPATH="${projectDir}/src"; \\
    else \\
        export PYTHONPATH="${projectDir}/src:\$PYTHONPATH"; \\
    fi
    
    dependency_input_vector.py  \
        input/gene_pairs.tsv \
        input/dependency_table.tsv \
        > dependency_input_vector.tsv
    """ 

}


/*
get available files from the DepMap API endpoint
<https://depmap.org/portal/api/download/files>

.META:
     1  release
     2  release_date
     3  filename
     4  url
     5  md5_hash
*/
process download_depmap_file_list {

    publishDir "${out_dir}",
                pattern: "databases/depMap/*.csv",
                mode: 'copy'

    output:
        path 'databases/depMap/files.csv' 

    script:
    """
    mkdir -p databases/depMap/
    
    curl 'https://depmap.org/portal/api/download/files' \\
        > databases/depMap/files.csv
    """

}


/*
download mutations deemed to be likely loss of function (LoF) for the cell
lines given the API endopoint file list information

data release DepMap Public 25Q3

- ModelID: "ACH-*" Model IDs
- IsDefaultEntryForModel: whether or not each sequencing is selected to represent the model in model-level datasets
- IsDefaultEntryForMC: whether or not each sequencing is selected to represent the model condition in model condition-level datasets
- StrippedCellLineName
- DepMapCode
- Lineage
- ModelConditionID
- SourceModelCondition
- CellFormat
- GrowthMedia
- GrowthPattern
- ProfileID
- SequencingPlatform
- SequencingID: unique identifier of this table
- DataType
- Stranded
- SharedToDbGaP
- SequencingDate
*/
process download_depmap_lof_mutations {

    publishDir "${out_dir}",
                pattern: "databases/depMap/*.csv",
                mode: 'copy'

    input:
        path 'input/files.csv'

    output:
        path 'databases/depMap/lof_mutations.csv'

    script:
    """
    mkdir -p databases/depMap/
    
    url=\$(cat input/files.csv \\
            | awk -F "," '\$1=="DepMap Public 25Q3"' \\
            | awk -F "," '\$3=="OmicsSomaticMutationsMatrixDamaging.csv"' \\
            | awk -F "," '{print \$4}')

    curl "\${url}" > databases/depMap/lof_mutations.csv
    """

}


/*
download mutations deemed to be likely cancer drivers
given the API endopoint file list information

data release DepMap Public 25Q3

*/
process download_depmap_hotspot_mutations {

    publishDir "${out_dir}",
                pattern: "databases/depMap/*.csv",
                mode: 'copy'

    input:
        path 'input/files.csv'

    output:
        path 'databases/depMap/hotspot_mutations.csv'

    script:
    """
    mkdir -p databases/depMap/
    
    url=\$(cat input/files.csv \\
            | awk -F "," '\$1=="DepMap Public 25Q3"' \\
            | awk -F "," '\$3=="OmicsSomaticMutationsMatrixHotspot.csv"' \\
            | awk -F "," '{print \$4}')

    curl "\${url}" > databases/depMap/hotspot_mutations.csv
    """
}


/*
download gene dependency probability estimates for all models in the integrated
gene effect, given the API endopoint file list information

data release DepMap Public 25Q3
*/
process download_depmap_gene_dependency {

    publishDir "${out_dir}",
                pattern: "databases/depMap/*.csv",
                mode: 'copy'

    input:
        path 'input/files.csv'

    output:
        path 'databases/depMap/gene_dependency.csv'

    script:
    """
    mkdir -p databases/depMap/
    
    url=\$(cat input/files.csv \\
            | awk -F "," '\$1=="DepMap Public 25Q3"' \\
            | awk -F "," '\$3=="CRISPRGeneDependency.csv"' \\
            | awk -F "," '{print \$4}')

    curl "\${url}" > databases/depMap/gene_dependency.csv
    """

}


/*
download metadata describing all cancer models/cell lines which are referenced
by a dataset contained within the DepMap portal, given the API endopoint file
list information

data release DepMap Public 25Q3

.META:
     1  ModelID
     2  PatientID
     3  CellLineName
     4  StrippedCellLineName
     5  DepmapModelType
     6  OncotreeLineage
     7  OncotreePrimaryDisease
     8  OncotreeSubtype
     9  OncotreeCode
    10  PatientSubtypeFeatures
    11  RRID
    12  Age
    13  AgeCategory
    14  Sex
    15  PatientRace
    16  PrimaryOrMetastasis
    17  SampleCollectionSite
    18  SourceType
    19  SourceDetail
    20  CatalogNumber
    21  ModelType
    22  TissueOrigin
    23  ModelDerivationMaterial
    24  ModelTreatment
    25  PatientTreatmentStatus
    26  PatientTreatmentType
    27  PatientTreatmentDetails
    28  Stage
    29  StagingSystem
    30  PatientTumorGrade
    31  PatientTreatmentResponse
    32  GrowthPattern
    33  OnboardedMedia
    34  FormulationID
    35  SerumFreeMedia
    36  PlateCoating
    37  EngineeredModel
    38  EngineeredModelDetails
    39  CulturedResistanceDrug
    40  PublicComments
    41  CCLEName
    42  HCMIID
    43  PediatricModelType
    44  ModelAvailableInDbgap
    45  ModelSubtypeFeatures
    46  WTSIMasterCellID
    47  SangerModelID
    48  COSMICID
    49  ModelIDAlias
*/
process download_depmap_models_info {

    publishDir "${out_dir}",
                pattern: "databases/depMap/*.csv",
                mode: 'copy'

    input:
        path 'input/files.csv'

    output:
        path 'databases/depMap/models.csv'

    script:
    """
    mkdir -p databases/depMap/
    
    url=\$(cat input/files.csv \\
            | awk -F "," '\$1=="DepMap Public 25Q3"' \\
            | awk -F "," '\$3=="Model.csv"' \\
            | awk -F "," '{print \$4}')

    curl "\${url}" > databases/depMap/models.csv
    """

}


/*
only select default lof entries for models
*/
process make_seed_list {

    input:
        path 'input/lof_mutations.csv'
        val n

    output:
        path "rwr_lof/seeds/*.txt"

    script:
    """
    mkdir -p rwr_lof/seeds

    mkdir -p seeds_tmp

    cat input/lof_mutations.csv \\
        | sed -n '1p' | cut -d ',' -f 7- | tr ',' '\\n' | awk 'NF' \\
        > genes.txt

    cat input/lof_mutations.csv \\
        | awk -F ',' '\$5=="Yes"' \\
        | sed '1d' | cut -d ',' -f 3,7- | tr ',' '\\t' | awk 'NF' \\
        > mutations.tsv

    while IFS=\$'\\t' read -r id mutations; do \\
        echo \${mutations} | tr ' ' '\\n' > seeds_tmp/\${id}.txt; \\
        paste genes.txt seeds_tmp/\${id}.txt > seeds_tmp//\${id}.tsv; \\
        rm seeds_tmp/\${id}.txt; \\
        awk -F "\\t" '\$2>=${n}{print \$1}' seeds_tmp//\${id}.tsv > rwr_lof/seeds/\${id}.txt; \\
        rm seeds_tmp//\${id}.tsv
        done < mutations.tsv
    """

}


/*
only select default lof entries for models
*/
process make_seed_list_eq {

    input:
        path 'input/lof_mutations.csv'
        val n

    output:
        path "rwr_lof/seeds/*.txt"

    script:
    """
    mkdir -p rwr_lof/seeds

    mkdir -p seeds_tmp

    cat input/lof_mutations.csv \\
        | sed -n '1p' | cut -d ',' -f 7- | tr ',' '\\n' | awk 'NF' \\
        > genes.txt

    cat input/lof_mutations.csv \\
        | awk -F ',' '\$5=="Yes"' \\
        | sed '1d' | cut -d ',' -f 3,7- | tr ',' '\\t' | awk 'NF' \\
        > mutations.tsv

    while IFS=\$'\\t' read -r id mutations; do \\
        echo \${mutations} | tr ' ' '\\n' > seeds_tmp/\${id}.txt; \\
        paste genes.txt seeds_tmp/\${id}.txt > seeds_tmp//\${id}.tsv; \\
        rm seeds_tmp/\${id}.txt; \\
        awk -F "\\t" '\$2==${n}{print \$1}' seeds_tmp//\${id}.tsv > rwr_lof/seeds/\${id}.txt; \\
        rm seeds_tmp//\${id}.tsv
        done < mutations.tsv
    """

}


/*
edit the header of input/lof_mutations.csv so that gene identifier,
reported as:
"<Gene_Synonym or Gene_Name> (<GeneID>)"
become:
"<Gene_Name>" 
*/
process translate_depmap_lof_table {

    publishDir "${out_dir}",
                pattern: "databases/depMap/*.csv",
                mode: 'copy'

    input:
        path 'input/lof_mutations.csv'
        path 'input/dict.tsv'

    output:
        path 'databases/depMap/lof_mutations_translated.csv'

    script:
    """
    mkdir -p databases/depMap/
    
    cat input/lof_mutations.csv \\
        | sed -n '1p' | cut -d ',' -f7- | tr "," "\\n" | awk 'NF' \\
        | tr -s " " "\\t" | tr -d '()' \\
        > target_name_id.tsv

    tr.py \\
        input/dict.tsv \\
        target_name_id.tsv \\
        1 \\
        3 \\
        1 \\
        1 \\
        > target_name_id_translated.tsv

    paste -d "," \\
        <(cat input/lof_mutations.csv  | sed -n '1p' | cut -d ',' -f-6) \\
        <(cat target_name_id_translated.tsv | cut -f1 | awk 'NF' | tr "\\n" "," ) \\
        | awk 'NF' > header_translated.csv

    cat \\
        header_translated.csv \\
        <(sed '1d' input/lof_mutations.csv) \\
        > databases/depMap/lof_mutations_translated.csv
    """

}


/*
edit the header of input/gene_dependency.csv so that gene identifiers,
reported as:
"<Gene_Synonym or Gene_Name> (<GeneID>)"
become:
"<Gene_Name>" 
*/
process translate_depmap_dependency_info {

    publishDir "${out_dir}",
                pattern: "databases/depMap/*.csv",
                mode: 'copy'

    input:
        path 'input/gene_dependency.csv'
        path 'input/dict.tsv'

    output:
        path 'databases/depMap/gene_dependency_translated.csv'

    script:
    """
    mkdir -p databases/depMap/
    
    cat input/gene_dependency.csv \\
        | sed -n '1p' | tr "," "\\n" | awk 'NF' \\
        | tr -s " " "\\t" | tr -d '()' \\
        > target_name_id.tsv

    tr.py \\
        input/dict.tsv \\
        target_name_id.tsv \\
        1 \\
        3 \\
        1 \\
        1 \\
        > target_name_id_translated.tsv

    cat target_name_id_translated.tsv \\
        | cut -f1 | awk 'NF' | tr "\\n" "," | sed 's/,\$//' \\
         > header_translated.csv

    cat \\
        <(awk '{print ","\$0}' header_translated.csv) \\
        <(sed '1d' input/gene_dependency.csv) \\
        > databases/depMap/gene_dependency_translated.csv
    """

}
