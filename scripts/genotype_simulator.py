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
import numpy as np
import os
import pandas as pd
import pyranges as pr
import subprocess
import tables
import tqdm

class GenotypeSimulator():
    def __init__(self, args):
        self.data_path = args.data_path
        self.data_identifier = args.data_identifier
        self.annotation_name = args.annotation_name
        self.features = args.features
        self.gfilter = "_".join(self.features)
        self.matrix_name = args.matrix_name
        self.snplist_name = args.snplist_name
        self.risk_rare = args.risk_rare
        self.ignore_gene_map = args.ignore_gene_map
        self.memory_cautious = args.memory_cautious
        self.matrix_chunk_size = args.matrix_chunk_size
        self.separator = args.separator
        
    def simulate_genotype(self):
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
        feature_names = "_".join(features)
        result_name = "annotations_{}.csv".format(feature_names)
        if os.path.isfile(self.data_path + result_name):
            return pd.read_csv(self.data_path + result_name, sep =" ")
        annotation = pr.read_gtf(annotation_path, as_df=True)
        annotation = annotation.loc[annotation['Feature'].isin(features)]
        annotation = annotation.reset_index(drop=True)
        annotation.drop(annotation.columns[range(5,25)], axis = 1, inplace = True)
        annotation.insert(5, "Feature ID", [i for i in range(annotation.shape[0])])
        annotation.drop(annotation.columns[[1]], axis = 1, inplace = True)
        annotation.to_csv(self.data_path + result_name, sep=" ", index_label=False)
        print("Saved {}".format(self.data_path + result_name))
        return annotation
    
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
        sf, risk_alleles = self.script_result_clean()
        sf.insert(5, column = "Risk Allele", value = risk_alleles)
        sf['Position'] = sf['Position'].apply(lambda x: json.loads(x) if type(x) == str else x) # data type conversion str to int
        sf.insert(6, "Feature ID", [-1]*sf.shape[0])
        if not self.ignore_gene_map:
            sf["Feature ID"] = sf["Feature ID"].astype('object') #so we can use lists
            ant_result_name = "annotations_{}.csv".format(feature_map)
            if not os.path.isfile(self.data_path + ant_result_name):
                annotations = self.clean_annotations([feature_map])
            else:
                annotations = pd.read_csv(self.data_path + ant_result_name, sep=" ")
            gene_map = self.get_feature_map(sf, annotations, feature = feature_map)
            for snp_id in range(sf.shape[0]):
                genes = gene_map[snp_id]
                if len(genes) == 1:
                    sf.at[snp_id, "Feature ID"] = genes[0]
                else:
                     sf.at[snp_id, "Feature ID"] = genes
        sf.to_csv(self.data_path + "snplist_{}_{}.csv".format(self.data_identifier, self.gfilter), sep=" ", index_label=False)
        print("Success!")
        return sf
    
    def load_dtypes(self):
        """
        Parse through header of RAW file to assign dtypes to decrease memory consuption of Pandas calls.
        """
        subprocess.check_call("head -n 1 {} > {}".format(self.data_path + self.matrix_name, self.data_path + "matrix_header.txt"), shell=True)
        header = pd.read_csv(self.data_path +"matrix_header.txt", sep=self.separator, header = None)
        types = [str, str]
        types += [int]*(header.shape[1]-2)
        dtypes = {}
        for name, typ in zip(header.iloc[0], types):
            dtypes["{}".format(name)] = typ
        if self.memory_cautious: # if memory cautious we read in the first column to give us our number of patients
            header = pd.read_csv(self.data_path + self.matrix_name, sep=self.separator, usecols=[0])
            self.num_patients = header.shape[0]
        return dtypes
    
    def script_result_clean(self):
        """
        Matrix name is a person X SNP RAW file. Note input snplist must have same order as raw file.
        Load in the RAW file and convert it to H5 format for Phenotype Simulation.
        Load in snplist and prepare for augmentation.
        """
        dtypes = self.load_dtypes()
        num_feat = len(dtypes) - 6
        f = tables.open_file(self.data_path + "genotype_{}_{}.h5".format(self.data_identifier, self.gfilter), mode='a')
        risk_alleles = []
        if self.memory_cautious:
            array_c = f.create_earray(f.root, 'data', tables.IntCol(), (0,num_feat), expectedrows=self.num_patients,filters=tables.Filters(complib='zlib', complevel=1))
            for df in pd.read_csv(self.data_path + self.matrix_name, iterator=True, sep=self.separator, chunksize=self.matrix_chunk_size, dtype=dtypes):
                df.drop(['IID', 'FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis = 1, inplace = True)
                self.h5_append(df,f)
                if self.risk_rare:
                    if len(risk_alleles) == 0:
                        risk_alleles = [[] for _ in range(df.shape[1])]
                    for idx in range(df.shape[1]):
                        val = 0
                        if df.iloc[:,idx].value_counts().to_list()[0] > df.shape[0]//2:
                            val = 1
                        risk_alleles[idx].append(val)
                else:
                    risk_alleles = [1 for idx in range(df.shape[1])]
            if self.risk_rare:
                for idx in range(len(risk_alleles)):
                    total = sum(risk_alleles[idx])
                    if total >= len(risk_alleles[idx]) // 2:
                        risk_alleles[idx] = 1
                    else:
                        risk_alleles[idx] = 0              
        else:
            df = pd.read_csv(self.data_path + self.matrix_name,  sep=self.separator, dtype=dtypes)
            df.drop(['IID', 'FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis = 1, inplace = True)
            num_pat = df.shape[0]
            num_feat = df.shape[1]
            array_c = f.create_earray(f.root, 'data', tables.IntCol(), (0,num_feat), expectedrows=num_pat,filters=tables.Filters(complib='zlib', complevel=1))
            self.h5_append(df,f)
            if self.risk_rare:
                risk_alleles = [1 if df.iloc[:,idx].value_counts().to_list()[0] > df.shape[0]//2 else 0 for idx in range(df.shape[1])]
            else:
                risk_alleles = [1 for idx in range(df.shape[1])]
            
        f.close()   
        sf = pd.read_csv(self.data_path + self.snplist_name, sep=self.separator, header = None)
        sf.drop(sf.columns[[2]], inplace = True, axis=1) # drop centimorgan value of .bim format
        sf.columns = ['Chromosome','ID','Position', 'Allele 1', 'Allele 2']
        return sf, risk_alleles
    
    def h5_append(self, df, f):
        """
        vectorize and append dataframe 1 patient at a time
        """
        num_pat = df.shape[0]
        for pat in tqdm.tqdm(range(num_pat)):
            a = df.iloc[pat,:].to_numpy()
            a=np.reshape(a, (1,-1))
            f.root.data.append(a)
            
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
        for snp_ind in range(data.shape[0]):
            pos=data["Position"][snp_ind]
            chromosome = data["Chromosome"][snp_ind]
            entry = []
            annotation = annotations[annotations["Chromosome"] == "chr{}".format(chromosome)]
            for gene_ind in annotation['Feature ID']:
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
