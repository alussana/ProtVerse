#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
Generate a dataset to test the model by adding more negative examples to the
test dataset, after discarding the already existing negative label examples.
*/
process add_examples {

    publishDir "${out_dir}",
                pattern: 'training/*.pkl',
                mode: 'copy'

    input:
        path 'input/dataset.tsv.gz'
        path 'input/X_test.pkl'
        path 'input/y_test.pkl'

    output:
        path "training/X_test_nNeg${params.additional_N_neg_examples}.pkl", emit: X_test
        path "training/y_test_nNeg${params.additional_N_neg_examples}.pkl", emit: y_test

    script:
        """
        mkdir -p training
        add_examples_to_test_set.py \
            input/dataset.tsv.gz \
            input/X_test.pkl \
            input/y_test.pkl \
            training/X_test_nNeg${params.additional_N_neg_examples}.pkl \
            training/y_test_nNeg${params.additional_N_neg_examples}.pkl

        echo "0"
        """

}

/*
Test RF performance on unbalanced dataset
ROC, PR curve; MCC
Use omics features
*/
process rf_omics_predict {

    cpus "${params.rf_n_jobs}"
    memory '32G'

    publishDir "${out_dir}",
                pattern: 'training/*.pdf',
                mode: 'copy'

    input:
        path 'input/X_test.pkl'
        path 'input/y_test.pkl'
        path 'input/model.pkl'

    output:
        path "training/roc_test_omics_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}_nNeg${params.additional_N_neg_examples}.pdf"
        path "training/prc_test_omics_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}_nNeg${params.additional_N_neg_examples}.pdf"

    script:
        """
        mkdir -p training
        test_RF_omics.py \
            input/X_test.pkl \
            input/y_test.pkl \
            input/model.pkl \
            training/roc_test_omics_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}_nNeg${params.additional_N_neg_examples}.pdf \
            training/prc_test_omics_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}_nNeg${params.additional_N_neg_examples}.pdf
        """
}

/*
Test RF performance on unbalanced dataset
ROC, PR curve; MCC
Use stringdb features
*/
process rf_stringdb_predict {

    cpus "${params.rf_n_jobs}"
    memory '32G'

    publishDir "${out_dir}",
                pattern: 'training/*.pdf',
                mode: 'copy'

    input:
        path 'input/X_test.pkl'
        path 'input/y_test.pkl'
        path 'input/model.pkl'

    output:
        path "training/roc_test_stringdb_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}_nNeg${params.additional_N_neg_examples}.pdf"
        path "training/prc_test_stringdb_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}_nNeg${params.additional_N_neg_examples}.pdf"

    script:
        """
        mkdir - p training
        test_RF_stringdb.py \
            input/X_test.pkl \
            input/y_test.pkl \
            input/model.pkl \
            training/roc_test_stringdb_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}_nNeg${params.additional_N_neg_examples}.pdf \
            training/prc_test_stringdb_maxDepth${params.rf_max_depth}_nEst${params.rf_n_estimators}_minSplit${params.rf_min_samples_split}_nNeg${params.additional_N_neg_examples}.pdf
        """
}