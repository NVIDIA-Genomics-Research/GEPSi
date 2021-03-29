#
# Copyright (c) 2021, NVIDIA CORPORATION.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#
"""Unit tests for Package."""
import os

import pandas as pd
from argparse import Namespace
import numpy as np
import pytest

from scripts.phenotype_simulator import PhenotypeSimulator

def test_phenotype_simulation(CHR=0,data_path="./",data_identifier="test",prefilter="exon",
                                      phenotype_experiement_name="playground_example",interactive_cut=0.2,mask_rate=0.1,
                                      dosage_frac=0.5, max_interaction_coeff=2, causal_snp_mode="gene", noise_scalar=1,
                                      heritability=1, phenotype_threshold=50, n_causal_snps=100, causal_gene_cut=0.005, max_gene_risk=5,
                                      total_snps=1000, total_genes=100, total_people=100):
    """
    Runs end to end test of phenotype simulation. Validates size and typing for all major results.
    """
    
    args = Namespace(chr=CHR, data_path=data_path, data_identifier=data_identifier,
                         prefilter=prefilter, phenotype_experiement_name=phenotype_experiement_name,
                         interactive_cut=interactive_cut, mask_rate=mask_rate, dosage_frac=dosage_frac,
                         max_interaction_coeff=max_interaction_coeff, causal_snp_mode=causal_snp_mode,
                         noise_scalar=noise_scalar, heritability=heritability, phenotype_threshold=phenotype_threshold,
                         n_causal_snps=n_causal_snps, causal_gene_cut=causal_gene_cut, max_gene_risk=max_gene_risk)
    
    pheno_sim = PhenotypeSimulator(args)
    
    genotype_data_test(pheno_sim, total_snps, total_genes, total_people)
    
    phenotype, data, causal_snps_idx, effect_size, interactive_snps = pheno_sim.simulate_phenotype()
    
    phenotype_test(phenotype, total_people)
    
    causality_test(effect_size, causal_snps_idx)
    
    interaction_test(data, interactive_snps, causal_snps_idx)

def interaction_test(data, interactive_snps, causal_snps_idx):
    assert type(interactive_snps) == dict
    for causal_snp in interactive_snps.keys():
        assert causal_snp in causal_snps_idx
        info = interactive_snps[causal_snp]
        assert len(info) == 3
        other_snp, other_risk = info[0], info[2]
        assert other_snp >= 0 and other_snp <= data.shape[0]
        assert data.iloc[other_snp]['Risk Allele'] == other_risk
    
def causality_test(effect_size, causal_snps_idx):
    assert type(effect_size) == dict
    assert type(causal_snps_idx) == list
    for snp in effect_size.keys():
        assert snp in causal_snps_idx
        assert len(effect_size[snp]) == 3
    
    
def phenotype_test(phenotype, total_people):
    assert len(phenotype) == total_people
    assert min(phenotype) == 0
    assert max(phenotype) == 1

    
def genotype_data_test(pheno_sim, snps, genes, people):
    data = pheno_sim.read_genotype_data()
    expected_shape = (snps, people+5)
    assert data.shape == expected_shape
    expected_columns = ['Feature ID', 'Position', 'Allele 1', 'Allele 2', 'Risk Allele']+[str(x) for x in range(people)]
    for dc, ec in zip(data.columns, expected_columns):
        assert dc == ec
        
    def condition(x):
        if type(x)== list:
            if min(x) < 0 or max(x) >= genes:
                return 1
            return 0
        elif x < 0 or x >= genes:
            return 1
        return 0
    
    invalid_gene_mapping = sum(data['Feature ID'].apply(condition))
    assert invalid_gene_mapping == 0
    for person in range(people):
        assert data[str(person)].value_counts().shape[0] == 2
