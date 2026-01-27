#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Generate N negative examples picking N random pairs of genes from the list of
"tractable" genes

.META:
1   bool (0)
2   Reactome ID (NA)
3   Gene 1
4   Gene 2
*/
process neg_examples_bck {

    cpus "${params.starmap_n_proc_buildNet}"

    input:
        path 'input/gene_dict.tsv'
        val N

    output:
        path 'neg.tsv'

    script:
    """
    cat input/gene_dict.tsv | sed '1d' | cut -f1 | awk 'NF' > genes.txt
    
    #N=\$(echo "${N} * 10" | bc)
    N=${N}

    for i in \$(seq 1 \$N); do \
        sample_negative_example.sh genes.txt & \
    done \
    | grep -v -w "Done" > neg.tsv
    """

}


/*
Generate N negative examples picking N random pairs of genes from the list of
genes annotated in Reactome

.META:
1   bool (0)
2   Reactome ID (NA)
3   Gene 1
4   Gene 2
*/
process neg_examples {

    cpus "${params.starmap_n_proc_buildNet}"

    input:
        path 'input/gene_dict.tsv'
        val N

    output:
        path 'neg.tsv'

    script:
    """
    cat input/gene_dict.tsv | sed '1d' | awk '\$NF!="NA"{print \$1}' > genes.txt
    
    #N=\$(echo "${N} * 10" | bc)
    N=${N}

    sample_gene_pairs.sh -i genes.txt -n \$N -o neg_pairs.tsv --seed 42

    cat neg_pairs.tsv | awk '{print "0\tNA\t"\$0}' > neg.tsv
    """

}


/*
Merge the partial input vectors from each of the databases

.META: training/*tsv.gz
1   index (<Gene Name 1>_<Gene Name 2>)
2   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
3   <feature 1>
4   <feature 2>
5   <feature 3>
... ...
*/
process features_tables {

    publishDir "${out_dir}",
                pattern: 'training/examples.tsv.gz',
                mode: 'copy'

    input:
        path 'input/eprot_*.tsv'
        path 'input/proteomehd_*.tsv'
        path 'input/mitchell2023_*.tsv'
        path 'input/gtex_*.tsv'
        path 'input/ptmdb_*.tsv'
        path 'input/dependency_*.tsv'
        path 'input/ubiquitination_*.tsv'
        path 'input/orthogroup_*.tsv'
        path 'input/humap_*.tsv'
        path 'input/lopit2025_*.tsv'
        //path 'input/ivkaphe_*.tsv'

    output:
        path "training/examples.tsv.gz"

    script:
    """
    mkdir -p training
    
    cat input/eprot_*.tsv | gzip > eprot.tsv.gz
    cat input/proteomehd_*.tsv | gzip > proteomehd.tsv.gz
    cat input/mitchell2023_*.tsv | gzip > mitchell2023.tsv.gz
    cat input/gtex_*.tsv | gzip > gtex.tsv.gz
    cat input/ptmdb_*.tsv | gzip > ptmdb.tsv.gz
    cat input/dependency_*.tsv | gzip > dependency.tsv.gz
    cat input/ubiquitination_*.tsv | gzip > ubiquitination.tsv.gz
    cat input/orthogroup_*.tsv | gzip > orthogroup.tsv.gz
    cat input/humap_*.tsv | gzip > humap.tsv.gz
    cat input/lopit2025_*.tsv | gzip > lopit2025.tsv.gz
    #cat input/ivkaphe_*.tsv | gzip > ivkaphe.tsv.gz

    cat \
        <(zcat gtex.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat gtex.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/gtex.tsv.gz

    cat \
        <(zcat eprot.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat eprot.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/eprot.tsv.gz

    cat \
        <(zcat proteomehd.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat proteomehd.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/proteomehd.tsv.gz

    cat \
        <(zcat mitchell2023.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat mitchell2023.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/mitchell2023.tsv.gz

    cat \
        <(zcat ptmdb.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat ptmdb.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/ptmdb.tsv.gz

    cat \
        <(zcat ubiquitination.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat ubiquitination.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/ubiquitination.tsv.gz

    cat \
        <(zcat orthogroup.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat orthogroup.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/orthogroup.tsv.gz

    cat \
        <(zcat dependency.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat dependency.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/dependency.tsv.gz

    cat \
        <(zcat lopit2025.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat lopit2025.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/lopit2025.tsv.gz

    #cat \
    #    <(zcat ivkaphe.tsv.gz | sed -n '1p' | cut -f1,3- \
    #        | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
    #    <(zcat ivkaphe.tsv.gz \
    #        | sed '1d' | grep -v "label" | cut -f1,3- \
    #        | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
    #        | sort | uniq) \
    #    | gzip > training/ivkaphe.tsv.gz

    cat \
        <(zcat humap.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat humap.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/humap.tsv.gz

    join_features.py \
        training/eprot.tsv.gz \
        training/proteomehd.tsv.gz \
        training/mitchell2023.tsv.gz \
        training/gtex.tsv.gz \
        training/ptmdb.tsv.gz \
        training/ubiquitination.tsv.gz \
        training/dependency.tsv.gz \
        training/orthogroup.tsv.gz \
        training/humap.tsv.gz \
        training/lopit2025.tsv.gz \
        | gzip > training/examples.tsv.gz
    """

}


/*
Merge the partial input vectors from each of the databases

.META: training/*tsv.gz
1   index (<Gene Name 1>_<Gene Name 2>)
2   label (boolean: "1" if gene 1 and gene 2 participate in the same reaction)
3   <feature 1>
4   <feature 2>
5   <feature 3>
... ...
*/
process features_tables_metabolism {

    publishDir "${out_dir}",
                pattern: 'training/examples_metabolism.tsv.gz',
                mode: 'copy'

    input:
        path 'input/eprot_*.tsv'
        path 'input/proteomehd_*.tsv'
        path 'input/mitchell2023_*.tsv'
        path 'input/gtex_*.tsv'
        path 'input/ptmdb_*.tsv'
        path 'input/dependency_*.tsv'
        path 'input/ubiquitination_*.tsv'
        path 'input/orthogroup_*.tsv'
        path 'input/humap_*.tsv'
        path 'input/lopit2025_*.tsv'

    output:
        path "training/examples_metabolism.tsv.gz"

    script:
    """
    mkdir -p training
    
    cat input/eprot_*.tsv | gzip > eprot.tsv.gz
    cat input/proteomehd_*.tsv | gzip > proteomehd.tsv.gz
    cat input/mitchell2023_*.tsv | gzip > mitchell2023.tsv.gz
    cat input/gtex_*.tsv | gzip > gtex.tsv.gz
    cat input/ptmdb_*.tsv | gzip > ptmdb.tsv.gz
    cat input/dependency_*.tsv | gzip > dependency.tsv.gz
    cat input/ubiquitination_*.tsv | gzip > ubiquitination.tsv.gz
    cat input/orthogroup_*.tsv | gzip > orthogroup.tsv.gz
    cat input/humap_*.tsv | gzip > humap.tsv.gz
    cat input/lopit2025_*.tsv | gzip > lopit2025.tsv.gz

    cat \
        <(zcat gtex.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat gtex.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/gtex.tsv.gz

    cat \
        <(zcat eprot.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat eprot.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/eprot.tsv.gz

    cat \
        <(zcat proteomehd.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat proteomehd.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/proteomehd.tsv.gz

    cat \
        <(zcat mitchell2023.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat mitchell2023.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/mitchell2023.tsv.gz

    cat \
        <(zcat ptmdb.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat ptmdb.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/ptmdb.tsv.gz

    cat \
        <(zcat ubiquitination.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat ubiquitination.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/ubiquitination.tsv.gz

    cat \
        <(zcat orthogroup.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat orthogroup.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/orthogroup.tsv.gz

    cat \
        <(zcat dependency.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat dependency.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/dependency.tsv.gz

    cat \
        <(zcat humap.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat humap.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/humap.tsv.gz

    cat \
        <(zcat lopit2025.tsv.gz | sed -n '1p' | cut -f1,3- \
            | awk '{print "index\t"\$0}' | cut -f1,2,5-) \
        <(zcat lopit2025.tsv.gz \
            | sed '1d' | grep -v "label" | cut -f1,3- \
            | awk 'NF{print \$2"_"\$3"\t"\$0}' | cut -f1,2,5- \
            | sort | uniq) \
        | gzip > training/lopit2025.tsv.gz

    join_features.py \
        training/eprot.tsv.gz \
        training/proteomehd.tsv.gz \
        training/mitchell2023.tsv.gz \
        training/gtex.tsv.gz \
        training/ptmdb.tsv.gz \
        training/ubiquitination.tsv.gz \
        training/dependency.tsv.gz \
        training/orthogroup.tsv.gz \
        training/humap.tsv.gz \
        training/lopit2025.tsv.gz \
        | gzip > training/examples_metabolism.tsv.gz
    """

}


/*
Plot PC1,PC2 and PC3,PC4 and PC5,PC6 for features in class 1 and 0 separately
*/
process features_pca {

    publishDir "${out_dir}",
            pattern: 'training/*.pdf',
            mode: 'copy'

    input:
        path 'input/examples.tsv.gz'

    output:
        path 'training/*.pdf'

    script:
    """
    mkdir -p training
    features_pca.py \
        input/examples.tsv.gz \
        training/pca12_omics_features.pdf \
        training/pca34_omics_features.pdf \
        training/pca56_omics_features.pdf \
        training/pca12_stringdb_features.pdf \
        training/pca34_stringdb_features.pdf \
        training/pca56_stringdb_features.pdf
    """

}


/*
Plot PC1,PC2 and PC3,PC4 and PC5,PC6 for features in class 1 and 0 separately
*/
process features_pca_metabolism {

    publishDir "${out_dir}",
            pattern: 'training_metabolism/*.pdf',
            mode: 'copy'

    input:
        path 'input/examples.tsv.gz'

    output:
        path 'training_metabolism/*.pdf'

    script:
    """
    mkdir -p training_metabolism
    features_pca.py \
        input/examples.tsv.gz \
        training_metabolism/pca12_omics_features.pdf \
        training_metabolism/pca34_omics_features.pdf \
        training_metabolism/pca56_omics_features.pdf \
        training_metabolism/pca12_stringdb_features.pdf \
        training_metabolism/pca34_stringdb_features.pdf \
        training_metabolism/pca56_stringdb_features.pdf
    """

}


/*
Split examples in 25% test set and 75% training set, maintaining label
frequencies
*/
process split_dataset {

    publishDir "${out_dir}",
            pattern: 'training/*.pkl',
            mode: 'copy'

    input:
        path 'input/examples.tsv.gz'

    output:
        path "training/X_train.pkl", emit: X_train
        path "training/X_test.pkl", emit: X_test
        path "training/y_train.pkl", emit: y_train
        path "training/y_test.pkl", emit: y_test

    script:
    """
    mkdir -p training
    split_folds.py \
        input/examples.tsv.gz \
        training/X_train.pkl \
        training/X_test.pkl \
        training/y_train.pkl \
        training/y_test.pkl
    """

}


/*
Split examples in 25% test set and 75% training set, maintaining label
frequencies
*/
process split_dataset_metabolism {

    publishDir "${out_dir}",
            pattern: 'training/*.pkl',
            mode: 'copy'

    input:
        path 'input/examples.tsv.gz'

    output:
        path "training/X_train_metabolism.pkl", emit: X_train
        path "training/X_test_metabolism.pkl", emit: X_test
        path "training/y_train_metabolism.pkl", emit: y_train
        path "training/y_test_metabolism.pkl", emit: y_test

    script:
    """
    mkdir -p training
    split_folds.py \
        input/examples.tsv.gz \
        training/X_train_metabolism.pkl \
        training/X_test_metabolism.pkl \
        training/y_train_metabolism.pkl \
        training/y_test_metabolism.pkl
    """

}


/*
inspect missingness of data sources in training and test sets,
comparing positive and negative labels; generate histograms
*/
process inspect_examples {

    publishDir "${out_dir}",
                pattern: 'training/*.pdf',
                mode: 'copy'

    input:
        path 'input/examples_X_train.pkl'
        path 'input/examples_X_test.pkl'
        path 'input/examples_y_train.pkl'
        path 'input/examples_y_test.pkl'

    output:
        path "training/inspect_switches_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}.pdf"

    script:
    """
    mkdir -p training
    inspect_switches.py \
        input/examples_X_train.pkl \
        input/examples_X_test.pkl \
        input/examples_y_train.pkl \
        input/examples_y_test.pkl \
        training/inspect_switches_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}.pdf
    """

}


/*
Compute SHAP values on test set and plot beeswarm violin plot
*/
process shap {

    publishDir "${out_dir}",
                pattern: 'training/*',
                mode: 'copy'

    memory '32G'

    input:
        path 'input/examples_X_test.pkl'
        path 'input/examples_y_test.pkl'
        path 'input/model.pkl'

    output:
        path "training/*.png"
        path "training/*.pkl"

    script:
    """
    mkdir -p training
    shap_values.py \
        input/examples_X_test.pkl \
        input/examples_y_test.pkl \
        input/model.pkl \
        training/shap_values.pkl \
        training/shap_test_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}.png \
    """

}

/*
calibrate the rf based on the balanced testing set
*/
process calibrate {

    publishDir "${out_dir}",
                pattern: 'training/*',
                mode: 'copy'

    input:
        path 'input/X_test.pkl'
        path 'input/y_test.pkl'
        path 'input/model.pkl'

    output:
        path "training/forest_calibrated_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}.pkl", emit: model
        path 'training/*.pdf'

    script:
    """
    mkdir -p training
    
    calibrate_RF.py \
        input/X_test.pkl \
        input/y_test.pkl \
        input/model.pkl \
        training/calibration_before_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}.pdf \
        training/calibration_after_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}.pdf \
        training/forest_calibrated_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}.pkl
    """
}

/*
create violin plots of feature distributions splitted by example's class
*/
process features_distrib {

    memory "32G"

    publishDir "${out_dir}",
                pattern: 'training/*.pdf',
                mode: 'copy'

    input:
        path 'input/examples.tsv.gz'

    output:
        path "training/*.pdf"

    script:
    """
    mkdir -p training

    features_distrib.py \\
        input/examples.tsv.gz \\
        training/features_distrib_h_1of2.pdf \\
        training/features_distrib_h_2of2.pdf \\
        training/features_distrib_v.pdf
    """

}


/*
create violin plots of feature distributions splitted by example's class
*/
process features_distrib_metabolism {

    memory "32G"

    publishDir "${out_dir}",
                pattern: 'training/*.pdf',
                mode: 'copy'

    input:
        path 'input/examples.tsv.gz'

    output:
        path "training/*.pdf"

    script:
    """
    mkdir -p training

    features_distrib.py \\
        input/examples.tsv.gz \\
        training/features_distrib_metabolism_h_1of2.pdf \\
        training/features_distrib_metabolism_h_2of2.pdf \\
        training/features_distrib_metabolism_v.pdf
    """

}


/*
Train and evaluate a XGBoost classifier 
*/
process signalling_train_valid_xgb_bck {

    publishDir "${out_dir}",
                pattern: 'training_signalling_xgb/*.pdf',
                mode: 'copy'

    memory '64G'
    cpus "${params.signalling_xgb_n_jobs}"

    publishDir "${out_dir}",
                pattern: 'training_signalling_xgb/*.pdf',
                mode: 'copy'

    publishDir "${out_dir}",
                pattern: 'training_signalling_xgb/*.pkl',
                mode: 'copy'

    input:
        path 'input/examples_X_train.pkl'
        path 'input/examples_X_valid.pkl'
        path 'input/examples_y_train.pkl'
        path 'input/examples_y_valid.pkl'

    output:
        path "training_signalling_xgb/*.pkl", emit: model
        path "training_signalling_xgb/*.pdf"

    script:
    """
    mkdir -p training_signalling_xgb

    signalling_train_valid_xgb.py \
        ${params.signalling_xgb_n_jobs} \
        ${params.signalling_xgb_max_depth} \
        ${params.signalling_xgb_n_estimators} \
        ${params.signalling_xgb_lr} \
        input/examples_X_train.pkl \
        input/examples_X_valid.pkl \
        input/examples_y_train.pkl \
        input/examples_y_valid.pkl \
        training_signalling_xgb/roc_train_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/roc_valid_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/prc_train_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/prc_valid_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/forest_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pkl \
        training_signalling_xgb/mdi_importance_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/perm_importance_train_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/perm_importance_valid_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf

    exit 0
    """
}


/*
Train and evaluate a XGBoost classifier 
*/
process signalling_train_valid_xgb {

    publishDir "${out_dir}",
                pattern: 'training_signalling_xgb/*.pdf',
                mode: 'copy'

    memory '64G'
    cpus "${params.signalling_xgb_n_jobs}"

    publishDir "${out_dir}",
                pattern: 'training_signalling_xgb/*.pdf',
                mode: 'copy'

    publishDir "${out_dir}",
                pattern: 'training_signalling_xgb/*.pkl',
                mode: 'copy'

    input:
        path 'input/examples_X_train.pkl'
        path 'input/examples_X_valid.pkl'
        path 'input/examples_y_train.pkl'
        path 'input/examples_y_valid.pkl'

    output:
        path "training_signalling_xgb/*.pkl", emit: model
        path "training_signalling_xgb/*.pdf"

    script:
    """
    mkdir -p training_signalling_xgb

    signalling_train_valid_xgb.py \
        ${params.signalling_xgb_n_jobs} \
        ${params.signalling_xgb_max_depth} \
        ${params.signalling_xgb_n_estimators} \
        ${params.signalling_xgb_lr} \
        input/examples_X_train.pkl \
        input/examples_X_valid.pkl \
        input/examples_y_train.pkl \
        input/examples_y_valid.pkl \
        training_signalling_xgb/roc_train_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/roc_valid_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/prc_train_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/prc_valid_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/forest_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pkl \
        training_signalling_xgb/mdi_importance_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/perm_importance_train_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf \
        training_signalling_xgb/perm_importance_valid_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}.pdf
    """
}


/*
Train XGB models using only indivdual sources and plot the ROC and PRC curves
*/
process signalling_xgb_single_source {

    publishDir "${out_dir}",
                pattern: 'training_signalling_xgb/*.pdf',
                mode: 'copy'

    cpus "${params.signalling_xgb_n_jobs}"
    memory '64G'

    input:
        path 'input/examples_X_train.pkl'
        path 'input/examples_X_valid.pkl'
        path 'input/examples_y_train.pkl'
        path 'input/examples_y_valid.pkl'

    output:
        path "training_signalling_xgb/*.pdf"

    script:
    """
    mkdir -p training_signalling_xgb

    signalling_xgb_single_source.py \
        ${params.signalling_xgb_n_jobs} \
        ${params.signalling_xgb_max_depth} \
        ${params.signalling_xgb_n_estimators} \
        ${params.signalling_xgb_lr} \
        input/examples_X_train.pkl \
        input/examples_X_valid.pkl \
        input/examples_y_train.pkl \
        input/examples_y_valid.pkl \
        training_signalling_xgb/roc_valid_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}_single_source.pdf \
        training_signalling_xgb/roc_valid_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}_single_feature.pdf \
        training_signalling_xgb/prc_valid_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}_single_source.pdf \
        training_signalling_xgb/prc_valid_maxDepth${params.signalling_xgb_max_depth}_nEst${params.signalling_xgb_n_estimators}_lr${params.signalling_xgb_lr}_single_feature.pdf
    """

}


/*
Train and evaluate a XGBoost classifier 
*/
process metabolism_train_valid_xgb {

    publishDir "${out_dir}",
                pattern: 'training_metabolism_xgb/*.pdf',
                mode: 'copy'

    memory '64G'
    cpus "${params.metabolism_xgb_n_jobs}"

    publishDir "${out_dir}",
                pattern: 'training_metabolism_xgb/*.pdf',
                mode: 'copy'

    publishDir "${out_dir}",
                pattern: 'training_metabolism_xgb/*.pkl',
                mode: 'copy'

    input:
        path 'input/examples_X_train.pkl'
        path 'input/examples_X_valid.pkl'
        path 'input/examples_y_train.pkl'
        path 'input/examples_y_valid.pkl'

    output:
        path "training_metabolism_xgb/*.pkl", emit: model
        path "training_metabolism_xgb/*.pdf"

    script:
    """
    mkdir -p training_metabolism_xgb

    metabolism_train_valid_xgb.py \
        ${params.metabolism_xgb_n_jobs} \
        ${params.metabolism_xgb_max_depth} \
        ${params.metabolism_xgb_n_estimators} \
        ${params.metabolism_xgb_lr} \
        input/examples_X_train.pkl \
        input/examples_X_valid.pkl \
        input/examples_y_train.pkl \
        input/examples_y_valid.pkl \
        training_metabolism_xgb/roc_train_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}.pdf \
        training_metabolism_xgb/roc_valid_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}.pdf \
        training_metabolism_xgb/prc_train_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}.pdf \
        training_metabolism_xgb/prc_valid_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}.pdf \
        training_metabolism_xgb/forest_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}.pkl \
        training_metabolism_xgb/mdi_importance_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}.pdf \
        training_metabolism_xgb/perm_importance_train_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}.pdf \
        training_metabolism_xgb/perm_importance_valid_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}.pdf
    """
}


/*
Train XGB models using only indivdual sources and plot the ROC and PRC curves
*/
process metabolism_xgb_single_source {

    publishDir "${out_dir}",
                pattern: 'training_metabolism_xgb/*.pdf',
                mode: 'copy'

    cpus "${params.metabolism_xgb_n_jobs}"
    memory '64G'

    input:
        path 'input/examples_X_train.pkl'
        path 'input/examples_X_valid.pkl'
        path 'input/examples_y_train.pkl'
        path 'input/examples_y_valid.pkl'

    output:
        path "training_metabolism_xgb/*.pdf"

    script:
    """
    mkdir -p training_metabolism_xgb

    metabolism_xgb_single_source.py \
        ${params.metabolism_xgb_n_jobs} \
        ${params.metabolism_xgb_max_depth} \
        ${params.metabolism_xgb_n_estimators} \
        ${params.metabolism_xgb_lr} \
        input/examples_X_train.pkl \
        input/examples_X_valid.pkl \
        input/examples_y_train.pkl \
        input/examples_y_valid.pkl \
        training_metabolism_xgb/roc_valid_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}_single_source.pdf \
        training_metabolism_xgb/roc_valid_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}_single_feature.pdf \
        training_metabolism_xgb/prc_valid_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}_single_source.pdf \
        training_metabolism_xgb/prc_valid_maxDepth${params.metabolism_xgb_max_depth}_nEst${params.metabolism_xgb_n_estimators}_lr${params.metabolism_xgb_lr}_single_feature.pdf
    """

}