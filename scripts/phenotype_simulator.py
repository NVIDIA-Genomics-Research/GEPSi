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
import tqdm

class PhenotypeSimulator():
    def __init__(self, args):
        self.chr = args.chr
        self.data_path = args.data_path
        self.data_identifier = args.data_identifier
        self.prefilter = args.prefilter
        self.phenotype_experiement_name = args.phenotype_experiement_name
        self.interactive_cut = args.interactive_cut
        self.mask_rate = args.mask_rate
        self.dosage_frac = args.dosage_frac
        self.max_interaction_coeff = args.max_interaction_coeff
        self.causal_snp_mode = args.causal_snp_mode
        self.noise_scalar = args.noise_scalar
        self.heritability  = args.heritability 
        self.phenotype_threshold = args.phenotype_threshold
        if self.causal_snp_mode == "random":
            self.n_causal_snps = args.n_causal_snps
        if self.causal_snp_mode == "gene":
            self.causal_gene_cut = args.causal_gene_cut
            self.max_gene_risk = args.max_gene_risk
        self.patient_chunk = args.patient_chunk
        
    def get_number_patients(self):
        #known to have 5 metadata columns
        df = pd.read_csv(self.data_path + "chr{}/chr{}_genotype_{}_{}.csv".format(self.chr, self.chr, self.data_identifier, self.prefilter), sep=" ", nrows=1)
        self.num_patients = df.shape[1]-5
        return self.num_patients
                         
    def get_person(self, person):
        if type(person) == list:
            patient = pd.read_csv(self.data_path + "chr{}/chr{}_genotype_{}_{}.csv".format(self.chr, self.chr, self.data_identifier, self.prefilter), sep=" ", usecols=[str(pear) for pear in person])
        else:
            patient = pd.read_csv(self.data_path + "chr{}/chr{}_genotype_{}_{}.csv".format(self.chr, self.chr, self.data_identifier, self.prefilter), sep=" ", usecols=[str(person)])
        return patient

    def read_genotype_data(self):
        """
        Reads in the annotated genotype csv file.
        """
        df = pd.read_csv(self.data_path + "chr{}/chr{}_genotype_{}_{}.csv".format(self.chr, self.chr, self.data_identifier, self.prefilter), sep=" ", usecols=['Risk Allele','Feature ID'], dtype={'Feature ID': object})
        df['Feature ID'] = df['Feature ID'].apply(lambda x: json.loads(x) if type(x) == str else x)
        self.risk_alleles = df['Risk Allele'].to_list()
        self.feature_id = df['Feature ID'].to_list()
        print("SNPs: {} People: {}".format(len(self.risk_alleles), self.get_number_patients()))
        return self.risk_alleles, self.feature_id
    
    def save_file(self, fname, data):
        with open(self.data_path + "chr{}/chr{}_{}_{}_{}_{}.pkl".format(self.chr, self.chr, fname, self.data_identifier, self.prefilter, self.phenotype_experiement_name), 'wb') as f:
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
        return phenotype, causal_snps_idx, effect_size, interactive_snps
    
    def heritability_injection(self, scores):
        """
        Given phenotype scores normalize and apply heritability function.
        Creates a score distrobution histogram plot as well as return the updated scores.
        """
        #Normalize scores
        scores = (scores - np.mean(scores))/np.std(scores)
        #Apply g' = h*g + k*sqrt(1-h^2)*N(0,1)
        phenotype_scores = self.heritability * scores + self.noise_scalar*np.random.randn(len(scores)) * np.sqrt(1 - self.heritability * self.heritability)
        self.get_distribution(phenotype_scores, title = "Chr {} Phenotype Scores with Heredity {}".format(self.chr, self.heritability), ylabel="Number of People", xlabel="Genetic Risk Score")
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
        interactive_snps_idx = np.random.choice(causal_snps_idx, size = max(1,int(self.interactive_cut*len(causal_snps_idx))), replace=False)
        interactive_snps = {}
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
        self.get_distribution(coeff, title = "Chr {} Interactive Coefficients".format(self.chr), ylabel="Number of SNPs", xlabel="Interation Coefficients")
        return interactive_snps

    def get_score(self, person, effect_size, interactive_snps):
        """
        Calculates the Phenotype Score of a input person.

        Parameters
        ----------
        person: Person column of Genotype file
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
        data: Dataframe that can be supplied directly from generate_genotype_file
        effect_size: Dictionary that maps causal snp indexes to length 3 list that has the effect size per genotyple value {0, 1, 2} accesed via index
        interactive_snps: Dictionary that maps causal snp indexes to a length 2 list [Interactive SNP Index Pair, Interaction Coefficeint]
        causal_snps_idx: List of SNP indices that are chosen as the causal SNPs
        cut: Fraction of causal SNPs to have an interacting pair
        mask_rate: Fraction of interactive realtions that are masking
        max_interaction_coeff: Upper bound of Uniform Distribution for interaction coefficient Sampling

        Results
        -------
        Generates interactive relations if not given already.
        Returns the calculated phenotype scores.
        """
        interactive_snps =  self.gen_snp_interaction(causal_snps_idx)
        num_patients = self.get_number_patients()
#         pheno_scores = [self.get_score(self.get_person(person), effect_size, interactive_snps) for person in range(num_patients)]
        pheno_scores = []
        for patient in range(0, num_patients,  self.patient_chunk):
            patient_list = list(range(patient, patient + self.patient_chunk))
            patients = self.get_person(patient_list)
            inner_pheno_scores = [self.get_score(patients[str(person)].to_list(), effect_size, interactive_snps) for person in patient_list]
            pheno_scores.extend(inner_pheno_scores)
        if num_patients % self.patient_chunk != 0:
            patient_list = list(range(num_patients-(num_patients % self.patient_chunk), num_patients))
            patients = self.get_person(patient_list)
            inner_pheno_scores = [self.get_score(patients[str(person)].to_list(), effect_size, interactive_snps) for person in patient_list]
            pheno_scores.extend(inner_pheno_scores)

        self.get_distribution(pheno_scores, title = "Chr {} Phenotype Scores".format(self.chr), ylabel="Number of People", xlabel="Genetic Risk Score")
        print(len(pheno_scores), "Phenotype scores")
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
        possible_snps = range(self.risk_alleles.shape[0])
        causal_snps_idx = np.random.choice(possible_snps, size = self.n_causal_snps, replace=False)
        effect_size = self.simulate_effect_sizes_random(causal_snps_idx)
        return effect_size, causal_snps_idx
    
    def simulate_effect_sizes_random(self, causal_snps_idx):
        """
        Simulates the effect size of causal SNPs.
        Half the causal SNPs get different values accoring to their genotype [0, 1, 2].

        Parameters
        ----------
        data: Dataframe that can be supplied directly from generate_genotype_file
        causal_snps_idx: List of SNP indices that are chosen as the causal SNPs
        data_path: Path to File created via generate_genotype_file to be read
        CHR: Chromosome number
        dosage_frac: Fraction of effect sizes that would be dosage dependant vs absolute

        Results
        -------
        Returns effect size dictionary for mapping of causal snps to index valued sampled effect sizes
        """
        effect_size = {}
        dosage = np.random.choice(causal_snps_idx, size = int(self.dosage_frac*len(causal_snps_idx)), replace=False)
        absolute = list(set(causal_snps_idx).difference(set(dosage)))
        for idx in dosage:
            a= np.random.normal()
            b= np.random.normal()
            if self.risk_allele[idx] == 1:
                effect_size[idx] = [0, min(a,b), max(a,b)]
            else:
                effect_size[idx] = [max(a,b), min(a,b), 0]
        for idx in absolute:
            a= np.random.normal()
            if self.risk_allele[idx] == 1:
                effect_size[idx] = [0, a, a]
            else:
                effect_size[idx] = [a, a, 0]

        self.save_file("effect_size", effect_size)
        return effect_size

    def get_unique_genes(self):
        """
        Given dataframe return list of unique genes present
        Can use gene to snps map to count all genes such that have non empty lists
        """
        gene_list = self.feature_id
        gene_set = set()
        for g in gene_list:
            if type(g) == int:
                gene_set.add(g)
            else:
                gene_set.update(g)
        return list(gene_set)

    def get_snps_from_gene(self, gene_id):
        """
        Given dataframe and gene_id returns list of snp indices that map to that gene
        can be replaced by second return of get gene map func
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
        plt.savefig(self.data_path +'chr{}/{}.png'.format(self.chr, title))
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
        effect_size = self.simulate_effect_sizes_gene(causal_snps_idx)
        self.save_file("causal_snp_idx", causal_snps_idx)
        self.save_file("causal_genes", causal_genes)
        return effect_size, list(causal_snps_idx.keys())

    def simulate_effect_sizes_gene(self, causal_snps_idx):
        """
        Simulates causal snps effect size where we apply a sampled fraction to the gene risk
        """
        effect_size = {}
        dosage = np.random.choice([k for k in causal_snps_idx.keys()], size = int(self.dosage_frac*len(causal_snps_idx.keys())), replace=False)
        absolute = list(set(causal_snps_idx.keys()).difference(set(dosage)))
        for idx in dosage:
            a= np.random.random()
            b= np.random.random()
            if self.risk_alleles[idx] == 1:
                effect_size[idx] = [0, min(a,b)*causal_snps_idx[idx], max(a,b)*causal_snps_idx[idx]]
            else:
                effect_size[idx] = [max(a,b)*causal_snps_idx[idx], min(a,b)*causal_snps_idx[idx], 0]
        for idx in absolute:
            a= np.random.normal()
            if self.risk_alleles[idx] == 1:
                effect_size[idx] = [0, a*causal_snps_idx[idx], a*causal_snps_idx[idx]]
            else:
                effect_size[idx] = [a*causal_snps_idx[idx], a*causal_snps_idx[idx], 0]

        self.save_file("effect_size", effect_size)
        return effect_size
                                                                                                