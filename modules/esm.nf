#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
use pretrained protein language model to make input embeddings
*/
process make_input_embeddings {

    label 'lightning'

    publishDir "${out_dir}",
                pattern: "input_embeddings/${esm_model}/*.pt",
                mode: 'copy'

    input:
        path 'input/dataset.tsv'
        val esm_model

    output:
        tuple file("input_embeddings/${esm_model}/x_train.npy"),
              file("input_embeddings/${esm_model}/y_train.npy"),
              file("input_embeddings/${esm_model}/x_valid.npy"),
              file("input_embeddings/${esm_model}/y_valid.npy"),
              val(esm_model)
    script:
    """
    export TORCH_HOME=${torch_cache_dir}

    mkdir -p input_embeddings/${esm_model}

    make_input_embeddings.py \
        input/dataset.tsv \
        ${esm_model} \
        input_embeddings/${esm_model}
    """

}