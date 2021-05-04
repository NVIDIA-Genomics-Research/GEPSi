import json
import os
import pickle
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyranges as pr
import scipy
import scipy.sparse
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
        print("Generating Genotype Files")
        sf, risk_alleles = self.script_result_clean()
        sf.insert(5, column = "Risk Allele", value = risk_alleles)
        sf['Position'] = sf['Position'].apply(lambda x: json.loads(x) if type(x) == str else x) # data type conversion str to int
        sf.insert(6, "Feature ID", [-1]*data.shape[0])
        if not self.ignore_gene_map:
            sf["Feature ID"] = sf["Feature ID"].astype('object') #so we can use lists
            ant_result_name = "annotations_{}.csv".format(feature_map)
            if not os.path.isfile(self.data_path + ant_result_name):
                annotations = self.clean_annotations([feature_map])
            else:
                annotations = pd.read_csv(self.data_path + ant_result_name, sep=" ")
            gene_map = self.get_feature_map(sf, annotations, feature = feature_map)
            for snp_id in range(data.shape[0]):
                genes = gene_map[snp_id]
                if len(genes) == 1:
                    sf.at[snp_id, "Feature ID"] = genes[0]
                else:
                     sf.at[snp_id, "Feature ID"] = genes
        sf.to_csv(self.data_path + "snplist_{}_{}.csv".format(self.data_identifier, self.gfilter), sep=" ", index_label=False)
        return data
    
    def load_dtypes(self):
        subprocess.check_call("head -n 1 {} > {}".format(self.data_path + self.matrix_name, self.data_path + "matrix_header.txt"), shell=True)
        header = pd.read_csv(self.data_path +"matrix_header.txt", sep=self.separator, header = None)
        types = [object, object]
        types += [int]*(header.shape[1]-2)
        dtypes = {}
        for name, typ in zip(header.iloc[0], types):
            dtypes["{}".format(name)] = typ
        return dtypes
    
    def script_result_clean(self):
        """
        Results
        -------
        Matrix name is a person X SNP matrix. Note input snplist must have same order as raw file.
        """
        dtypes = self.load_dtypes()
        if self.memory_cautious:
            tp = pd.read_csv(self.data_path + self.matrix_name, iterator=True, sep=self.separator, chunksize=self.matrix_chunk_size, dtype=dtypes)
            df = pd.concat(tp, ignore_index=True)
        else:
            df = pd.read_csv(self.data_path + self.matrix_name,  sep=self.separator, dtype=dtypes) #w
        print(df.head())
        df.drop(['IID', 'FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'], axis = 1, inplace = True)
        data = df.to_numpy()
        num_pat = data.shape[0]
        num_feat = data.shape[1]
        print("Genotype snps and patients", num_feat, num_pat)
        f = tables.open_file(self.data_path + "genotype_{}_{}.h5".format(self.data_identifier, self.gfilter), mode='w')
        array_c = f.create_earray(f.root, 'data', tables.IntCol(), (0,num_feat ), expectedrows=num_pat,filters=tables.Filters(complib='zlib', complevel=1))
        f.close()
        f = tables.open_file(self.data_path + "genotype_{}_{}.h5".format(self.data_identifier, self.gfilter), mode='a')
        for pat in tqdm.tqdm(range(num_pat)):
            a = data[pat,:]
            a=np.reshape(a, (1,-1))
            f.root.data.append(a)
        f.close()
        if self.risk_rare:
            risk_alleles = [1 if df.iloc[:,idx].value_counts().to_list()[0] > df.shape[0]//2 else 0 for idx in range(df.shape[1])]
        else:
            risk_alleles = [1 for idx in range(df.shape[1])]
        sf = pd.read_csv(self.data_path + self.snplist_name, sep=self.separator, header = None)
        sf.drop(sf.columns[[2]], inplace = True, axis=1) # drop centimorgan value
        sf.columns = ['CHR','ID','Position', 'Allele 1', 'Allele 2']
        return sf, risk_alleles
    
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
