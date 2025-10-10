#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
Fit best model architecture found in

Timo Reim, Anne Hartebrodt, David B Blumenthal, Judith Bernett, Markus List
Deep learning models for unbiased sequence-based PPI prediction plateau at an accuracy of 0.65
Bioinformatics, July 2025
https://doi.org/10.1093/bioinformatics/btaf192

This model uses ESM per-protein, t33 embeddings (650 million parameters)
*/
process fit_reim_model{

    publishDir "${out_dir}",
                pattern: 'models/reim2025/*',
                mode: 'copy'

    input:
        path "input/X_train.pkl"
        path "input/y_train.pkl"
        path "input/X_val.pkl"
        path "input/y_val.pkl"

    output:
        path "models/reim2025/model.pkl"

    script:
    """
    mkdir -p models/reim2025

    fit_reim_model.py
    """

}