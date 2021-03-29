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

class GenotypeSimulator():
    def __init__(self, args):
        self.chr = args.chr
        self.data_path = args.data_path
        self.data_identifier = args.data_identifier
        self.annotation_name = args.annotation_name
        self.features = args.features
        self.gfilter = "_".join(self.features)
        self.matrix_name = args.matrix_name
        self.snplist_name = args.snplist_name

    def simulate_genotype(self):
        print("Simulating Genotypes")
        return self.generate_genotype_file()

    def clean_annotations(self, features = ["gene", "transcript", "exon"] ):
        """
        Generates Annotation file for specifed chromosome with the desired features
        Parameters. If file exists read it and and return else create
        ----------
        data_path: Path to Annotations file
        annotation_name: Name of Annotation File
        CHR: Chromosome Number
        features: List of desired annotion types. Default results in no feature pruning
        Results
        -------
        Reads in general whole genome annotation file and filters it for the desired chromosome and featuures.
        The result is saved for later Gene Mapping.
        """
        print("Cleaning Annotations")
        annotation_path = self.data_path + self.annotation_name
        chromosome = "chr{}".format(self.chr)
        feature_names = "_".join(features)
        result_name = "{}_annotations_{}.csv".format(self.chr, feature_names)
        if os.path.isfile(self.data_path + result_name):
            return pd.read_csv(self.data_path + result_name, sep =" ")
        annotation = pr.read_gtf(annotation_path, as_df=True)
        ant = annotation.loc[annotation['Chromosome'] == chromosome]
        ant = ant.loc[ant['Feature'].isin(features)]
        ant = ant.reset_index(drop=True)
        ant.drop(ant.columns[range(5,25)], axis = 1, inplace = True)
        ant.insert(5, "Feature ID", [i for i in range(ant.shape[0])])
        ant.drop(ant.columns[[1]], axis = 1, inplace = True)
        ant.to_csv(self.data_path + result_name, sep=" ", index_label=False)
        print("Saved {}".format(self.data_path + result_name))
        return ant
    
    def generate_genotype_file(self, feature_map = "gene"):
        """
        This function takes in the paths of the resulting PLINK bash script products. 
        It saves and returns a gene mapped and filtered dataframe.
        Once the HAPGEN2 command is run. Use this function to receive a usable dataframe for phenptype simulation.
        Parameters
        ----------
        data_path: Path to Resutling Script Data Files
        annotation_name: Name to Annotation File
        CHR: Chromosome Number
        features: List of desired annotion types. Default results in no feature pruning
        save_name: optional parameter to have different versions of data from same simulation
        Results
        -------
        Reads in direct results of PLINK script
        Reads and cleans Gene annotations of desired chromosome and features.
        Generates and applies Gene Mapping.
        Saves the resutling file to be used for Modeling + Analysis
        """
        print("Generating Genotype Files")
        T = time.time()
        data, risk_alleles = self.script_result_clean()
        print("Read Time: {}".format((time.time()-T)/60))
        T = time.time()

        data.insert(3, column = "Risk Allele", value = risk_alleles)
        print("Risk Time: {}".format((time.time()-T)/60))
        data['Position'] = data['Position'].apply(lambda x: json.loads(x) if type(x) == str else x)
        # data type conversion str to int
        T = time.time()

        ant_result_name = "{}_annotations_{}.csv".format(self.chr, feature_map)
        if not os.path.isfile(self.data_path + ant_result_name):
            annotations = self.clean_annotations([feature_map])
        else:
            annotations = pd.read_csv(self.data_path + ant_result_name, sep=" ")
            
        print("Annotation Time: {}".format((time.time()-T)/60))
        T = time.time()
        gene_map = self.get_feature_map(data, annotations, feature = feature_map)
        data.insert(0, "Feature ID", range(data.shape[0]))
        data["Feature ID"] = data["Feature ID"].astype('object') #so we can use lists
        for snp_id in range(data.shape[0]):
            genes = gene_map[snp_id]
            if len(genes) == 1:
                data.at[snp_id, "Feature ID"] = genes[0]
            else:
                 data.at[snp_id, "Feature ID"] = genes
        print("Map Time: {}".format((time.time()-T)/60))
        T = time.time()
        subprocess.check_call("mkdir -p {}chr{}".format(self.data_path, self.chr))
        data.to_csv(self.data_path + "chr{}/chr{}_genotype_{}_{}.csv".format(self.chr, self.chr, self.data_identifier, self.gfilter), sep=" ", index_label=False)                   
        print("Saved Result {}".format(self.data_path + "chr{}/chr{}_genotype_{}_{}.csv".format(self.chr, self.chr, self.data_identifier, self.gfilter)))     
        return data
    
    def script_result_clean(self):
        """
        Results
        -------
        Matrix name is a person X SNP matrix
        """
        df = pd.read_csv(self.data_path+matrix_name, sep=' ')
        print(df.head())
        risk_alleles = [1 if df.iloc[:,idx].value_counts().to_list()[0] > df.shape[0]//2 else 0 for idx in range(0, df.shape[1])]
        df = df.transpose()
        df.reset_index(drop=True, inplace=True)
        sf = pd.read_csv(self.data_path+snplist_name, sep=' ', header = None)
        sf.columns = ['CHR','Position', 'Allele 1', 'Allele 2']
        sf.drop(sf.columns[[0]], axis=1,inplace=True)
        data = pd.concat([sf, df], axis = 1)
        print("Combined Script Results", data.shape)
        print(data.head())
        return data, risk_alleles
    
    def get_feature_map(self, data, annotations, feature = "gene"):
        """
        Creates a Gene map or list of all genes a SNP at a specific position falls into
        Parameters
        ----------
        data: Dataframe that can be supplied directly from generate_genotype_file
        annotations: Annotations dataframe 
        Results
        -------
        Returns gene map that is a list such that it contains lists for the corresponding genes with same indices as data's snps
        """
        annotations = annotations.loc[annotations['Feature'] == feature]
        snps_to_genes = []
        T = time.time()
        for snp_ind in range(data.shape[0]):
            pos=data["Position"][snp_ind]
            entry = []
            for gene_ind in range(annotations.shape[0]):
                gstart = annotations["Start"][gene_ind]
                gend = annotations["End"][gene_ind]
                if pos >= gstart and pos <= gend:
                    entry.append(gene_ind)
                elif gstart > pos:
                    break
            if len(entry) == 0:
                entry.append(-1)
            snps_to_genes.append(entry)
        return snps_to_genes
