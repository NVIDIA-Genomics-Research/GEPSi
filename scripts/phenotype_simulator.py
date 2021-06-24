#
# Copyright (c) 2021, NVIDIA CORPORATION & AFFILIATES.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#
import json
import os
import pickle
import subprocess
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyranges as pr
import scipy
import scipy.sparse
import tables
import tqdm

class PhenotypeSimulator():
    def __init__(self, args):
        self.data_path = args.data_path
        self.data_identifier = args.data_identifier
        self.phenotype_experiement_name = args.phenotype_experiement_name
        self.interactive_cut = args.interactive_cut 
        self.mask_rate = args.mask_rate 
        self.dominance_frac = args.dominance_frac  
        self.recessive_frac = args.recessive_frac
        self.max_interaction_coeff = args.max_interaction_coeff 
        self.causal_snp_mode = args.causal_snp_mode
        self.heritability = args.heritability 
        self.phenotype_threshold = args.phenotype_threshold
        if self.causal_snp_mode == "random":
            self.n_causal_snps = args.n_causal_snps 
        if self.causal_snp_mode == "gene":
            self.causal_gene_cut = args.causal_gene_cut 
            self.max_gene_risk = args.max_gene_risk 
        
    def get_number_patients(self):
        with tables.open_file(self.data_path + "genotype_{}.h5".format(self.data_identifier), "r") as f:
            self.num_patients = f.root.data.shape[0]
        return self.num_patients
                         
    def read_genotype_data(self):
        """
        Reads in the annotated genotype csv file.
        """
        df = pd.read_csv(self.data_path + "snplist_{}.csv".format(self.data_identifier), sep=" ", usecols=['Risk Allele','Feature ID'], dtype={'Feature ID': object})
        df['Feature ID'] = df['Feature ID'].apply(lambda x: json.loads(x) if type(x) == str else x)
        self.risk_alleles = df['Risk Allele'].to_list()
        self.feature_id = df['Feature ID'].to_list()
        self.num_snps = len(self.risk_alleles)
        print("SNPs: {} People: {}".format(len(self.risk_alleles), self.get_number_patients()))
        return self.risk_alleles, self.feature_id
    
    def save_file(self, fname, data):
        with open(self.data_path + "{}_{}_{}.pkl".format(fname, self.data_identifier, self.phenotype_experiement_name), 'wb') as f:
            pickle.dump(data, f)
    
    def simulate_phenotype(self):
        """
        Simulates causal SNP sampling and effect size creation if not done already.
        Generates the phenotype scores and returns resulting phenotype for the input patients

        This can be run directly after generate_genotype_file

        Parameters
        ----------
        Define kwargs based on initialized CAUSAL SNP MODE
        
        interactive_cut: Fraction of causal SNPs to have an interacting pair
        mask_rate: Fraction of interactive realtions that are masking
        dosage_frac: Fraction of effect sizes that would be dosage dependant vs absolute
        max_interaction_coeff: Max value for interaction coefficient uniform distribution

        Results
        -------
        Takes resulting genotype data once the genes have been mapped and simulates causal snps, interactive relations and returns
        a list of resulting Phenotypes
        """
        if self.causal_snp_mode == "random":
            effect_size, causal_snps_idx = self.simulate_causal_snps_random()
        elif self.causal_snp_mode == "gene":
            effect_size, causal_snps_idx = self.simulate_causal_snps_gene()
        phenotype_scores, interactive_snps = self.generate_phenotypes_scores(effect_size, causal_snps_idx)
        phenotype_scores = self.heritability_injection(phenotype_scores)
        phenotype_cutoff = np.percentile(phenotype_scores, self.phenotype_threshold)
        phenotype = [1 if x >= phenotype_cutoff else 0 for x in phenotype_scores]
        self.save_file("phenotype", phenotype)
        print("Success!")
        return phenotype, causal_snps_idx, effect_size, interactive_snps
    
    def heritability_injection(self, scores):
        """
        Given phenotype scores normalize and apply heritability function.
        Creates a score distrobution histogram plot as well as return the updated scores.
        """
        #Normalize scores
        scores = (scores - np.mean(scores))/np.std(scores)
        #Apply g' = h*g + sqrt(1-h^2)*N(0,1)
        phenotype_scores = self.heritability * scores + np.random.randn(len(scores)) * np.sqrt(1 - self.heritability * self.heritability)
        self.get_distribution(phenotype_scores, title = "Phenotype Scores with Heredity {} {}".format(self.heritability, self.phenotype_experiement_name), ylabel="Number of People", xlabel="Genetic Risk Score")
        return phenotype_scores
    
    def patient_level_score_injection(self, scores, patients, coefficients):
        #coefficients is a dictionary of labels to coefficients
        patient_level_bias = [coefficients[patient] for patient in patients]
        new_score = score + patient_level_bias
        return new_score
    
    def patient_level_func_injection(self, scores, patients, coefficients):
        #coefficients is a dictionary of labels to functions of genetic risk score
        patient_level_bias = [coefficients[patient](scores[i]) for i,patient in enumerate(patients)]
        new_score = score + patient_level_bias
        return new_score
        
    def gen_snp_interaction(self, causal_snps_idx):
        """
        Generates SNP-SNP interaction Coeffecients for Phentype score calculation.

        Parameters
        ----------
        data: Dataframe that can be supplied directly from generate_genotype_file
        n_causal_snps: Number of desired causal SNPs
        cut: Fraction of causal SNPs to have an interacting pair
        mask_rate: Fraction of interactive realtions that are masking
        max_interaction_coeff: Max value for interaction coefficient uniform distribution

        Results
        -------
        Returns the interactive_snps dictionary for pairing and interaction coeficients
        interactive_snps: Dictionary that maps causal snp indexes to a length 3 list [Interactive SNP Index Pair, Interaction Coefficeint, Partner Risk Allele]
        """
        interactive_snps = {}
        if self.interactive_cut == 0: #No epistasis
            return interactive_snps
        interactive_snps_idx = np.random.choice(causal_snps_idx, size = max(1,int(self.interactive_cut*len(causal_snps_idx))), replace=False)
        coeff = []
        for idx in interactive_snps_idx:
            interaction = []
            possible_snps = [i for i in range(len(self.risk_alleles)) if i != idx]
            partner = np.random.choice(possible_snps, size = 1, replace=False)[0]
            interaction.append(partner)
            if np.random.random() > self.mask_rate:
                interaction_coeff = np.random.uniform(0, self.max_interaction_coeff)
                interaction.append(interaction_coeff)
                coeff.append(interaction_coeff)
            else:
                interaction.append(0)
                coeff.append(0)
            interaction.append(self.risk_alleles[partner])
            interactive_snps[idx] = interaction
        self.save_file("interactive_snps", interactive_snps)
        self.get_distribution(coeff, title = "{} Interactive Coefficients".format(self.phenotype_experiement_name), ylabel="Number of SNPs", xlabel="Interaction Coefficients")
        return interactive_snps

    def get_score(self, person, effect_size, interactive_snps):
        """
        Calculates the Phenotype Score of a input person.

        Parameters
        ----------
        person: Person row of Genotype file
        effect_size: Dictionary that maps causal snp indexes to length 3 list that has the effect size per genotyple value {0, 1, 2} accesed via index
        interactive_snps: Dictionary that maps causal snp indexes to a length 3 list [Interactive SNP Index Pair, Interaction Coefficeint, Partner Risk Allele]

        Results
        -------
        Returns the sum of effective sizes as the phenotype score for the given person.
        """
        score = 0
        for idx in effect_size.keys():
            interaction_coeff = 1
            if idx in interactive_snps.keys():
                other_gen = person[interactive_snps[idx][0]]
                other_risk = interactive_snps[idx][2]
                if (other_risk == 0 and other_gen != 2) or (other_risk == 1 and other_gen != 0):
                    interaction_coeff = interactive_snps[idx][1]
                #  if the partner risk allele is present in 1 or 2 copies activate the effect
            score += effect_size[idx][person[idx]]*interaction_coeff
        return score

    def generate_phenotypes_scores(self, effect_size, causal_snps_idx):
        """
         Parameters
        ----------
        effect_size: Dictionary that maps causal snp indexes to length 3 list that has the effect size per genotyple value {0, 1, 2} accesed via index
        interactive_snps: Dictionary that maps causal snp indexes to a length 2 list [Interactive SNP Index Pair, Interaction Coefficeint]
        causal_snps_idx: Dictionary of SNP indices that are chosen as the causal SNPs to gene risk coefficients if present.

        Results
        -------
        Generates interactive relations if not given already.
        Returns the calculated phenotype scores.
        """
        interactive_snps =  self.gen_snp_interaction(causal_snps_idx)
        with tables.open_file(self.data_path + "genotype_{}.h5".format(self.data_identifier), "r") as f:
            patients = f.root.data
            pheno_scores = [self.get_score(patients[idx,:], effect_size, interactive_snps) for idx in range(patients.shape[0])]                        
        self.get_distribution(pheno_scores, title = "{} Phenotype Scores".format(self.phenotype_experiement_name), ylabel="Number of People", xlabel="Genetic Risk Score")
        print(len(pheno_scores), "Phenotype Scores Closed")
        return pheno_scores, interactive_snps
                                  
    def simulate_causal_snps_random(self):
        """
        Simulates causal SNP selection and generates effect sizes

        Parameters
        ----------
        data_path: Path to File created via generate_genotype_file to be read
        CHR: Chromosome Number
        n_causal_snps: Number of desired causal SNPs
        data: Dataframe that can be supplied directly from generate_genotype_file
        dosage_frac: Fraction of effect sizes that would be dosage dependant vs absolute

        Results
        -------
        Simulates the causal SNPs and samples the effect sizes for each.
        Returns the data, the effect_size dictionary, list of causal snp indices
        """
        self.read_genotype_data()
        possible_snps = list(range(len(self.risk_alleles)))
        causal_snps_id = np.random.choice(possible_snps, size = self.n_causal_snps, replace=False)
        causal_snps_idx = {idx: 1 for idx in causal_snps_id}
        effect_size = self.simulate_effect_sizes(causal_snps_idx)
        causal_genes = self.get_causal_genes(causal_snps_id)
        self.save_file("causal_genes", causal_genes)
        self.save_file("causal_snp_idx", causal_snps_idx)
        return effect_size, causal_snps_id

    def get_causal_genes(self, snp_ids):
        """
        Grabs genes associated for random causal snp selection.
        """
        causal_genes = {}
        for snp in snp_ids:
            genes = self.feature_id[snp]
            if type(genes) == list:
                for g in genes:
                    causal_genes[g] = 1
            else:
                causal_genes[genes] = 1
        return causal_genes
                    
    def get_unique_genes(self):
        """
        Given dataframe return list of unique genes present
        Can use gene to snps map to count all genes such that have non empty lists
        """
        gene_list = self.feature_id
        gene_set = set()
        for g in gene_list:
            if type(g) == int:
                if g < 0:
                    continue
                gene_set.add(g)
            else:
                gene_set.update(g)
        return list(gene_set)

    def get_snps_from_gene(self, gene_id):
        """
        Given dataframe and gene_id returns list of snp indices that map to that gene
        """
        gene_list = self.feature_id
        snp_idx_list = list(range(len(self.risk_alleles)))
        gsnps = []
        for snp, gene in zip(snp_idx_list,gene_list):
            if type(gene) == int:
                if gene_id == gene:
                    gsnps.append(snp)
            else:
                if gene_id in gene:
                    gsnps.append(snp)
        return gsnps
    
    def get_distribution(self, values, title=None, xlabel = None, ylabel = None, binsize = 0.1):
        plt.hist(values, bins=np.arange(np.min(values)-binsize, np.max(values)+binsize, binsize))
        plt.title(title)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.savefig(self.data_path +'{}.png'.format(title))
        plt.close()
        
    def simulate_causal_snps_gene(self):
        """
        Simulate causal genes and snps and return effect_size, list of causal snp indices
        """
        self.read_genotype_data()
        possible_genes = self.get_unique_genes()
        causal_genes_idx = np.random.choice(possible_genes, size = max(1,int(len(possible_genes)*self.causal_gene_cut)), replace=False)
        causal_snps_idx = {}
        causal_genes = {}
        for cgene in causal_genes_idx:
            gene_risk = np.random.uniform(0, self.max_gene_risk)
            causal_genes[cgene] = gene_risk
            causal_fraction = np.random.random()
            possible_snps = self.get_snps_from_gene(cgene)
            causal_snps = np.random.choice(possible_snps, size = max(1,int(len(possible_snps)*causal_fraction)), replace=False)
            for cs in causal_snps:
                causal_snps_idx[cs] = gene_risk    
        effect_size = self.simulate_effect_sizes(causal_snps_idx)
        self.save_file("causal_snp_idx", causal_snps_idx)
        self.save_file("causal_genes", causal_genes)
        return effect_size, list(causal_snps_idx.keys())
    
    def simulate_effect_sizes(self, causal_snps_idx):
        """
        Simulates causal snps effect size where we apply a sampled fraction to the gene risk
        """
        effect_size = {}
        self.dosage_frac = 1 - self.dominance_frac - self.recessive_frac
        dominant = np.random.choice([k for k in causal_snps_idx.keys()], size = int(self.dominance_frac*len(causal_snps_idx.keys())), replace=False)
        leftover = list(set(causal_snps_idx.keys()).difference(set(dominant)))
        recessive = np.random.choice([k for k in leftover], size = int(self.recessive_frac*len(causal_snps_idx.keys())), replace=False)
        dosage = list(set(leftover).difference(set(recessive)))
        
        for idx in dosage:
            a= np.random.random()
            b= 2*a
            if self.risk_alleles[idx] == 1:
                effect_size[idx] = [0, a*causal_snps_idx[idx], b*causal_snps_idx[idx]]
            else:
                effect_size[idx] = [b*causal_snps_idx[idx], a*causal_snps_idx[idx], 0]
        for idx in dominant:
            a= np.random.normal()
            if self.risk_alleles[idx] == 1:
                effect_size[idx] = [0, a*causal_snps_idx[idx], a*causal_snps_idx[idx]]
            else:
                effect_size[idx] = [a*causal_snps_idx[idx], a*causal_snps_idx[idx], 0]
        for idx in recessive:
            a= np.random.normal()
            if self.risk_alleles[idx] == 1:
                effect_size[idx] = [0, 0, a*causal_snps_idx[idx]]
            else:
                effect_size[idx] = [a*causal_snps_idx[idx], 0, 0]
        self.save_file("effect_size", effect_size)
        return effect_size
                                                                                                
