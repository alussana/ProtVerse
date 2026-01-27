#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Get PPI data leakage-free dataset from

Judith Bernett, David B Blumenthal, Markus List
"Cracking the black box of deep sequence-based protein–protein interaction
prediction"
Briefings in Bioinformatics
<https://doi.org/10.1093/bib/bbae076>
*/
process dl_bernett_2024 {

    publishDir "${out_dir}",
                pattern: 'databases/bernett2024/*',
                mode: 'copy'

    output:
        path "databases/bernett2024/21591618.zip"

    script:
    """
    mkdir -p databases/bernett2024
    
    wget --user-agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36" \\
         --referer="https://example.com/" \\
         --content-on-error --content-disposition \\
         -O databases/bernett2024/21591618.zip \\
         ${params.url_bernett_2024}
    """

}


/*
From 
<https://figshare.com/articles/dataset/PPI_prediction_from_sequence_gold_standard_dataset/21591618>

* Big dataset: 
  163,192 training points (Intra-1),
  59,260 validation points (Intra-0),
  52,048 test points (Intra-2)) 
  + corresponding protein sequences from Swissprot
* No direct data leakage: proteins from training are not contained in
  validation or test, proteins from validation are not in training or test,
  proteins from test are not in validation or training
* Minimized sequence similarity between training, validation, test because
  whole human proteome was split with KaHIP such that sequence similarities
  are minimized w.r.t. length-normalized bitscores
* Redundancy-reduction with CD-HIT: inside of the datasets, no proteins with
  >40% pairwise sequence similarity

translate all UniProt AC specified in the first tab-separated column of
input/dict.tsv with the corresponding gene names found in the third column

the second column of input/dict.tsv specifies whether the translation is to a
"Gene_Name" or a "Gene_Synonym"

UniProt AC ids are preferentially translated to Gene_Name; if a gene name is
not found in the dictionary, then we look for a Gene_Synonym 

do not keep untranslated rows

final data point counts are:
  52035 test_ids.tsv
  52048 test_ids_untr.tsv
 163134 train_ids.tsv
 163192 train_ids_untr.tsv
  59246 valid_ids.tsv
  59260 valid_ids_untr.tsv
*/
process parse_and_translate_bernett_2024 {

    publishDir "${out_dir}",
                pattern: 'databases/bernett2024/*',
                mode: 'copy'

    input:
        path "input/21591618.zip"
        path "input/dict.tsv"
    
    output:
        path "databases/bernett2024/train_ids_untr.tsv", emit: train_ids_untr
        path "databases/bernett2024/valid_ids_untr.tsv", emit: valid_ids_untr
        path "databases/bernett2024/test_ids_untr.tsv", emit: test_ids_untr
        path "databases/bernett2024/train_ids.tsv", emit: train_ids
        path "databases/bernett2024/valid_ids.tsv", emit: valid_ids
        path "databases/bernett2024/test_ids.tsv", emit: test_ids

    script:
    """
    mkdir -p databases/bernett2024

    unzip input/21591618.zip

    cat \
        <(cat Intra1_pos_rr.txt | awk '{print \$1"\\t"\$2"\\t"1}') \
        <(cat Intra1_neg_rr.txt | awk '{print \$1"\\t"\$2"\\t"0}') \
        > databases/bernett2024/train_ids_untr.tsv

    cat \
        <(cat Intra0_pos_rr.txt | awk '{print \$1"\\t"\$2"\\t"1}') \
        <(cat Intra0_neg_rr.txt | awk '{print \$1"\\t"\$2"\\t"0}') \
        > databases/bernett2024/valid_ids_untr.tsv

    cat \
        <(cat Intra2_pos_rr.txt | awk '{print \$1"\\t"\$2"\\t"1}') \
        <(cat Intra2_neg_rr.txt | awk '{print \$1"\\t"\$2"\\t"0}') \
        > databases/bernett2024/test_ids_untr.tsv

    translator_ac_2_gene_name_or_synonym.py \
        input/dict.tsv \
        databases/bernett2024/train_ids_untr.tsv \
        1 \
        3 \
        1,2 \
        0 \
        2 \
        | sort | uniq > databases/bernett2024/train_ids.tsv

    translator_ac_2_gene_name_or_synonym.py \
        input/dict.tsv \
        databases/bernett2024/valid_ids_untr.tsv \
        1 \
        3 \
        1,2 \
        0 \
        2 \
        | sort | uniq > databases/bernett2024/valid_ids.tsv
    
    translator_ac_2_gene_name_or_synonym.py \
        input/dict.tsv \
        databases/bernett2024/test_ids_untr.tsv \
        1 \
        3 \
        1,2 \
        0 \
        2 \
        | sort | uniq > databases/bernett2024/test_ids.tsv
    """
}


/*
concatenate input text files removing empty lines, then split in smaller files
of n lines each

before than, columns of input files are rearranged and a dummy column added
to match the existing format in the input of *_featvec() processes 
*/
process bernett2024_conform_cat_and_split {

    input:
        path "input/train_ids.tsv"
        path "input/valid_ids.tsv"
        path "input/test_ids.tsv"
        val n

    output:
        path 'cat_and_split/*.tsv'

    script:
    """
    cat input/train_ids.tsv \
        | awk '{print \$3"\\tx\\t"\$1"\\t"\$2}' \
        > train_ids.tsv

    cat input/valid_ids.tsv \
        | awk '{print \$3"\\tx\\t"\$1"\\t"\$2}' \
        > valid_ids.tsv

    cat input/test_ids.tsv \
        | awk '{print \$3"\\tx\\t"\$1"\\t"\$2}' \
        > test_ids.tsv

    mkdir -p cat_and_split

    cat *.tsv | awk 'NF' > cat.txt
    
    split -l ${n} -a 16 -x cat.txt cat_and_split/ --additional-suffix .tsv
    """

}


/*
Merge the partial input vectors from each of the databases

Exclude SREK1_ZRANB2: for unknown reasons it has 4 additional and empty columns

.META: training_bernett2024/*tsv.gz
1   index (<Gene Name 1>_<Gene Name 2>)
2   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
3   <feature 1>
4   <feature 2>
5   <feature 3>
... ...
*/
process bernett2024_features_tables {

    publishDir "${out_dir}",
                pattern: 'training_bernett2024/examples.tsv.gz',
                mode: 'copy'

    input:
        path 'input/eprot_*.tsv'
        path 'input/proteomehd_*.tsv'
        path 'input/mitchell2023_*.tsv'
        path 'input/gtex_*.tsv'
        path 'input/ptmdb_*.tsv'
        path 'input/ubiquitination_*.tsv'
        path 'input/dependency_*.tsv'
        path 'input/orthogroup_*.tsv'
        path 'input/humap3_*.tsv'
        path 'input/lopit2025_*.tsv'

    output:
        path "training_bernett2024/examples.tsv.gz"

    script:
    """
    mkdir -p training_bernett2024
    
    cat input/eprot_*.tsv | gzip > eprot.tsv.gz
    cat input/proteomehd_*.tsv | gzip > proteomehd.tsv.gz
    cat input/mitchell2023_*.tsv | gzip > mitchell2023.tsv.gz
    cat input/gtex_*.tsv | gzip > gtex.tsv.gz
    cat input/ptmdb_*.tsv | gzip > ptmdb.tsv.gz
    cat input/ubiquitination_*.tsv | gzip > ubiquitination.tsv.gz
    cat input/dependency_*.tsv | gzip > dependency.tsv.gz
    cat input/orthogroup_*.tsv | gzip > orthogroup.tsv.gz
    cat input/humap3_*.tsv | gzip > humap3.tsv.gz
    cat input/lopit2025_*.tsv | gzip > lopit2025.tsv.gz

    cat \
        <(zcat gtex.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat gtex.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training_bernett2024/gtex.tsv.gz

    cat \
        <(zcat eprot.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat eprot.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training_bernett2024/eprot.tsv.gz

    cat \
        <(zcat proteomehd.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat proteomehd.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training_bernett2024/proteomehd.tsv.gz

    cat \
        <(zcat mitchell2023.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat mitchell2023.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training_bernett2024/mitchell2023.tsv.gz

    cat \
        <(zcat ptmdb.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat ptmdb.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training_bernett2024/ptmdb.tsv.gz

    cat \
        <(zcat ubiquitination.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat ubiquitination.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training_bernett2024/ubiquitination.tsv.gz

    cat \
        <(zcat orthogroup.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat orthogroup.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training_bernett2024/orthogroup.tsv.gz

    cat \
        <(zcat dependency.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat dependency.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training_bernett2024/dependency.tsv.gz

    cat \
        <(zcat humap3.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat humap3.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training_bernett2024/humap3.tsv.gz

    cat \
        <(zcat lopit2025.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat lopit2025.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training_bernett2024/lopit2025.tsv.gz

    bernett2024_join_features.py \\
        training_bernett2024/eprot.tsv.gz \\
        training_bernett2024/proteomehd.tsv.gz \\
        training_bernett2024/mitchell2023.tsv.gz \\
        training_bernett2024/gtex.tsv.gz \\
        training_bernett2024/ptmdb.tsv.gz \\
        training_bernett2024/ubiquitination.tsv.gz \\
        training_bernett2024/dependency.tsv.gz \\
        training_bernett2024/orthogroup.tsv.gz \\
        training_bernett2024/humap3.tsv.gz \\
        training_bernett2024/lopit2025.tsv.gz \\
        | grep -v "SREK1_ZRANB2" \\
        | gzip > training_bernett2024/examples.tsv.gz
    """

}


/*
Plot PC1,PC2 and PC3,PC4 and PC5,PC6 for features in class 1 and 0 separately
*/
process bernett2024_features_pca {

    publishDir "${out_dir}",
            pattern: 'training_bernett2024/*.pdf',
            mode: 'copy'

    input:
        path 'input/examples.tsv.gz'

    output:
        path 'training_bernett2024/*.pdf'

    script:
    """
    mkdir -p training_bernett2024

    features_pca.py \
        input/examples.tsv.gz \
        training_bernett2024/pca12_omics_features.pdf \
        training_bernett2024/pca34_omics_features.pdf \
        training_bernett2024/pca56_omics_features.pdf \
        training_bernett2024/pca12_stringdb_features.pdf \
        training_bernett2024/pca34_stringdb_features.pdf \
        training_bernett2024/pca56_stringdb_features.pdf
    """

}


/*
[...]
*/
process bernett2024_make_data_splits {

    publishDir "${out_dir}",
            pattern: "training_bernett2024/*.pkl",
            mode: 'copy'

    input:
        path "input/examples.tsv.gz"
        path "input/train_ids.tsv"
        path "input/valid_ids.tsv"
        path "input/test_ids.tsv"

    output:
        path "training_bernett2024/X_train.pkl", emit: X_train
        path "training_bernett2024/X_valid.pkl", emit: X_valid
        path "training_bernett2024/X_test.pkl", emit: X_test
        path "training_bernett2024/y_train.pkl", emit: y_train
        path "training_bernett2024/y_valid.pkl", emit: y_valid
        path "training_bernett2024/y_test.pkl", emit: y_test

    script:
    """
    mkdir -p training_bernett2024

    cat input/train_ids.tsv | awk '{print \$1"_"\$2}' > train_ids.txt
    cat input/valid_ids.tsv | awk '{print \$1"_"\$2}' > valid_ids.txt
    cat input/test_ids.tsv | awk '{print \$1"_"\$2}' > test_ids.txt
    
    bernett2024_make_data_splits.py \
        input/examples.tsv.gz \
        train_ids.txt \
        valid_ids.txt \
        test_ids.txt \
        training_bernett2024/X_train.pkl \
        training_bernett2024/X_valid.pkl \
        training_bernett2024/X_test.pkl \
        training_bernett2024/y_train.pkl \
        training_bernett2024/y_valid.pkl \
        training_bernett2024/y_test.pkl
    """
}


/*
create violin plots of feature distributions splitted by example's class
*/
process bernett2024_features_distrib {

    publishDir "${out_dir}",
                pattern: 'training_bernett2024/*.pdf',
                mode: 'copy'

    input:
        path 'input/examples.tsv.gz'

    output:
        path "training_bernett2024/*.pdf"

    script:
    """
    mkdir -p training_bernett2024

    features_distrib.py \\
        input/examples.tsv.gz \\
        training_bernett2024/features_distrib_h_to_55.pdf \\
        training_bernett2024/features_distrib_h_from_55.pdf \\
        training_bernett2024/features_distrib_v.pdf
    """

}


/*
Train and evaluate a Random Forest classifier 
*/
process bernett2024_train_valid_rf {

    publishDir "${out_dir}",
                pattern: 'training_bernett2024/*.pdf',
                mode: 'copy'

    memory '64G'
    cpus "${params.bernett2024_rf_n_jobs}"

    publishDir "${out_dir}",
                pattern: 'training_bernett2024/*.pdf',
                mode: 'copy'

    publishDir "${out_dir}",
                pattern: 'training_bernett2024/*.pkl',
                mode: 'copy'

    input:
        path 'input/examples_X_train.pkl'
        path 'input/examples_X_valid.pkl'
        path 'input/examples_y_train.pkl'
        path 'input/examples_y_valid.pkl'

    output:
        path "training_bernett2024/*.pkl", emit: model
        path "training_bernett2024/*.pdf"

    script:
    """
    mkdir -p training_bernett2024
    bernett2024_train_valid_rf.py \
        ${params.bernett2024_rf_n_jobs} \
        ${params.bernett2024_rf_max_depth} \
        ${params.bernett2024_rf_n_estimators} \
        ${params.bernett2024_rf_min_samples_split} \
        input/examples_X_train.pkl \
        input/examples_X_valid.pkl \
        input/examples_y_train.pkl \
        input/examples_y_valid.pkl \
        training_bernett2024/roc_train_maxDepth${params.bernett2024_rf_max_depth}_nEst${params.bernett2024_rf_n_estimators}_minSplit${params.bernett2024_rf_min_samples_split}.pdf \
        training_bernett2024/roc_valid_maxDepth${params.bernett2024_rf_max_depth}_nEst${params.bernett2024_rf_n_estimators}_minSplit${params.bernett2024_rf_min_samples_split}.pdf \
        training_bernett2024/prc_train_maxDepth${params.bernett2024_rf_max_depth}_nEst${params.bernett2024_rf_n_estimators}_minSplit${params.bernett2024_rf_min_samples_split}.pdf \
        training_bernett2024/prc_valid_maxDepth${params.bernett2024_rf_max_depth}_nEst${params.bernett2024_rf_n_estimators}_minSplit${params.bernett2024_rf_min_samples_split}.pdf \
        training_bernett2024/forest_maxDepth${params.bernett2024_rf_max_depth}_nEst${params.bernett2024_rf_n_estimators}_minSplit${params.bernett2024_rf_min_samples_split}.pkl \
        training_bernett2024/mdi_importance_maxDepth${params.bernett2024_rf_max_depth}_nEst${params.bernett2024_rf_n_estimators}_minSplit${params.bernett2024_rf_min_samples_split}.pdf \
        training_bernett2024/perm_importance_train_maxDepth${params.bernett2024_rf_max_depth}_nEst${params.bernett2024_rf_n_estimators}_minSplit${params.bernett2024_rf_min_samples_split}.pdf \
        training_bernett2024/perm_importance_valid_maxDepth${params.bernett2024_rf_max_depth}_nEst${params.bernett2024_rf_n_estimators}_minSplit${params.bernett2024_rf_min_samples_split}.pdf
    """
}


/*
Train RF models using only indivdual features and plot the ROC and PRC curves
*/
process bernett2024_train_valid_rf_reduced {

    publishDir "${out_dir}",
                pattern: 'training_bernett2024/*.pdf',
                mode: 'copy'

    cpus "${params.bernett2024_rf_n_jobs}"
    memory '64G'

    input:
        path 'input/examples_X_train.pkl'
        path 'input/examples_X_valid.pkl'
        path 'input/examples_y_train.pkl'
        path 'input/examples_y_valid.pkl'

    output:
        path "training_bernett2024/*.pdf"

    script:
    """
    mkdir -p training_bernett2024

    bernett2024_train_valid_rf_reduced.py \
        ${params.bernett2024_rf_n_jobs} \
        ${params.bernett2024_rf_reduced_max_depth} \
        ${params.bernett2024_rf_reduced_n_estimators} \
        ${params.bernett2024_rf_min_samples_split} \
        input/examples_X_train.pkl \
        input/examples_X_valid.pkl \
        input/examples_y_train.pkl \
        input/examples_y_valid.pkl \
        training_bernett2024/roc_valid_reduced_maxDepth${params.bernett2024_rf_reduced_max_depth}_nEst${params.bernett2024_rf_reduced_n_estimators}_minSplit${params.bernett2024_rf_min_samples_split}.pdf \
        training_bernett2024/roc_valid_restricted_maxDepth${params.bernett2024_rf_reduced_max_depth}_nEst${params.bernett2024_rf_reduced_n_estimators}_minSplit${params.bernett2024_rf_min_samples_split}.pdf \
    """

}


/*
Train and evaluate a XGBoost classifier 
*/
process bernett2024_train_valid_xgb {

    publishDir "${out_dir}",
                pattern: 'training_bernett2024_xgb/*.pdf',
                mode: 'copy'

    memory '64G'
    cpus "${params.bernett2024_xgb_n_jobs}"

    publishDir "${out_dir}",
                pattern: 'training_bernett2024_xgb/*.pdf',
                mode: 'copy'

    publishDir "${out_dir}",
                pattern: 'training_bernett2024_xgb/*.pkl',
                mode: 'copy'

    input:
        path 'input/examples_X_train.pkl'
        path 'input/examples_X_valid.pkl'
        path 'input/examples_y_train.pkl'
        path 'input/examples_y_valid.pkl'

    output:
        path "training_bernett2024_xgb/*.pkl", emit: model
        path "training_bernett2024_xgb/*.pdf"

    script:
    """
    mkdir -p training_bernett2024_xgb

    bernett2024_train_valid_xgb.py \
        ${params.bernett2024_xgb_n_jobs} \
        ${params.bernett2024_xgb_max_depth} \
        ${params.bernett2024_xgb_n_estimators} \
        ${params.bernett2024_xgb_lr} \
        input/examples_X_train.pkl \
        input/examples_X_valid.pkl \
        input/examples_y_train.pkl \
        input/examples_y_valid.pkl \
        training_bernett2024_xgb/roc_train_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}.pdf \
        training_bernett2024_xgb/roc_valid_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}.pdf \
        training_bernett2024_xgb/prc_train_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}.pdf \
        training_bernett2024_xgb/prc_valid_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}.pdf \
        training_bernett2024_xgb/forest_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}.pkl \
        training_bernett2024_xgb/mdi_importance_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}.pdf \
        training_bernett2024_xgb/perm_importance_train_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}.pdf \
        training_bernett2024_xgb/perm_importance_valid_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}.pdf
    """
}


/*
Train XGB models using only indivdual sources and plot the ROC and PRC curves
*/
process bernett2024_xgb_single_source {

    publishDir "${out_dir}",
                pattern: 'training_bernett2024_xgb/*.pdf',
                mode: 'copy'

    cpus "${params.bernett2024_xgb_n_jobs}"
    memory '64G'

    input:
        path 'input/examples_X_train.pkl'
        path 'input/examples_X_valid.pkl'
        path 'input/examples_y_train.pkl'
        path 'input/examples_y_valid.pkl'

    output:
        path "training_bernett2024_xgb/*.pdf"

    script:
    """
    mkdir -p training_bernett2024_xgb

    bernett2024_xgb_single_source.py \
        ${params.bernett2024_xgb_n_jobs} \
        ${params.bernett2024_xgb_max_depth} \
        ${params.bernett2024_xgb_n_estimators} \
        ${params.bernett2024_xgb_lr} \
        input/examples_X_train.pkl \
        input/examples_X_valid.pkl \
        input/examples_y_train.pkl \
        input/examples_y_valid.pkl \
        training_bernett2024_xgb/roc_valid_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}_single_source.pdf \
        training_bernett2024_xgb/roc_valid_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}_single_feature.pdf \
        training_bernett2024_xgb/prc_valid_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}_single_source.pdf \
        training_bernett2024_xgb/prc_valid_maxDepth${params.bernett2024_xgb_max_depth}_nEst${params.bernett2024_xgb_n_estimators}_lr${params.bernett2024_xgb_lr}_single_feature.pdf
    """

}