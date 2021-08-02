#
# Copyright (c) 2021, NVIDIA CORPORATION & AFFILIATES.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#
"""Unit tests for Package."""

from argparse import Namespace
import pandas as pd

from scripts.phenotype_simulator import PhenotypeSimulator


def test_phenotype_simulation(data_path="./sample_data/",data_identifier="chr0_test",
                                      phenotype_experiment_name="playground_example",interactive_cut=0.2,mask_rate=0.1,
                                      dominance_frac=0.1,recessive_frac=0.1, max_interaction_coeff=2, causal_snp_mode="gene",
                                      heritability=1, case_frac=0.5, stratify=False, n_causal_snps=100, causal_gene_cut=0.05, max_gene_risk=5,
                                      total_snps=1000, total_genes=100, total_people=100):
    """
    Runs end to end test of phenotype simulation. Validates size and typing for all major results.
    """
    
    args = Namespace(data_path=data_path, data_identifier=data_identifier,
                         phenotype_experiment_name=phenotype_experiment_name,
                         interactive_cut=interactive_cut, mask_rate=mask_rate,
                         max_interaction_coeff=max_interaction_coeff, causal_snp_mode=causal_snp_mode,recessive_frac=recessive_frac,
                         heritability=heritability, case_frac=case_frac, stratify=stratify, dominance_frac=dominance_frac,
                         n_causal_snps=n_causal_snps, causal_gene_cut=causal_gene_cut, max_gene_risk=max_gene_risk)

    pheno_sim = PhenotypeSimulator(args)
    
    phenotype, causal_snps_idx, effect_size, interactive_snps = pheno_sim.simulate_phenotype()
    
    phenotype_validation(phenotype, total_people)
    
    causality_validation(effect_size, causal_snps_idx)
    
    snplist = pd.read_csv(data_path + "snplist_{}.csv".format(data_identifier), sep=" ")
    
    interaction_validation(snplist, interactive_snps, causal_snps_idx)


def interaction_validation(data, interactive_snps, causal_snps_idx):
    assert type(interactive_snps) == dict
    for causal_snp in interactive_snps.keys():
        assert causal_snp in causal_snps_idx
        info = interactive_snps[causal_snp]
        assert len(info) == 3
        other_snp, other_risk = info[0], info[2]
        assert other_snp >= 0 and other_snp <= data.shape[0]
        assert data.iloc[other_snp]['Risk Allele'] == other_risk


def causality_validation(effect_size, causal_snps_idx):
    assert type(effect_size) == dict
    assert type(causal_snps_idx) == list
    for snp in effect_size.keys():
        assert snp in causal_snps_idx
        assert len(effect_size[snp]) == 3


def phenotype_validation(phenotype, total_people):
    assert len(phenotype) == total_people
    assert min(phenotype) == 0
    assert max(phenotype) == 1
