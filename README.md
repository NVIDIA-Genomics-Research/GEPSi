# GWAS Data Simulation

**TODO: Name the Package**

PACKAGE is a toolkit to simulate phenotypes for GWAS analysis, given input genotype data for a population.

## Installation

### System requirements

* Python 3.6+

## Build from Source

### 1. Clone repository

#### Latest released version
This will clone the repo to the `main` branch, which contains code for latest released version
and hot-fixes.

```
git clone --recursive -b master https://github.com/clara-genomics/gwas-data-simulation-public.git
```

### 2. Install dependencies

Install Package and its associated dependencies from requirements.txt

```
pip install .
```
    
### 3. Tests

Run unit tests to verify that installation was successful

    ```
    python -m pytest tests/
    ```
    
## Workflow

### 1. Formatting genotype data

Genotype data should be supplied in a `.raw` format along with a `.snplist` file. PACKAGE gives us the ability to format the genotype data matrix and associated annotations into an annotated csv file.
    
```
    gwas-sim genotype -chr21 -data_path /DLGWAS/data/ --matrix_name genotype.raw --snplist_name genotype.snplist
```
    
Results in the creation of a `.csv` file containing an annotated SNP X Person matrix with Genotype Values of 0,1,2. This data matrix is needed to run the phenotype simulation. 

![Alt text](./annotated_matrix_result.png?raw=true "Title")

The `.raw` and `.snplist` files can be produced from other formats using [PLINK](https://www.cog-genomics.org/plink/). PLINK can also be used to filter SNPs within selected regions (exons, transcripts, or genes) as well as filter SNPs based on their allele frequencies. 

For example, we used the following PLINK v1.9 command to filter and format genotype data for human chromosome 21:

```
/plink \
  --gen gensim_chr21_100k.controls.gen.gz \
  --sample gensim_chr21_100k.sample \
  --maf 0.01 \
  --extract range <BED file containing exon positions for chr21> \
  --allow-no-sex \
  --snps-only \
  --write-snplist \
  --recode A \
  --oxford-single-chr 21 \
  --out genotype
    
```
Resulting in the creation of

    /DLGWAS/data/genotype.raw: a Person X SNP Genotype Matrix
    /DLGWAS/data/genotype.snplist: Meta data for each SNP


### 2. Generating Phenotypes

Create Phenotypes for generated phenotypes using default values.
```
gwas-sim phenotype -dp /DLGWAS/data/ -chr 21 --data_identifier 100k --prefilter exon --phenotype_experiement_name example_name
```
    
Results in the creation of

    /DLGWAS/data/chr21/chr21_phenotype_100k_exon_example_name.pkl
    /DLGWAS/data/chr21/chr21_effect_size_100k_exon_example_name.pkl
    /DLGWAS/data/chr21/chr21_interactive_snps_100k_exon_example_name.pkl
    /DLGWAS/data/chr21/chr21_causal_snp_idx_100k_exon_example_name.pkl
    /DLGWAS/data/chr21/chr21_causal_genes_100k_exon_example_name.pkl
    
**chr21_phenotype_100k_exon_example_name.pkl**: a list of binary phenotypes for each person defined by the Genotype Matrix <br />
**chr21_effect_size_100k_exon_example_name.pkl**: a dictionary with key SNP index and value a list of the dosage dependent effect sizes  <br />
**chr21_interactive_snps_100k_exon_example_name.pkl**: a dictionary that maps causal snp indices to a list of length 3 [Interactive SNP Index Pair, Interaction Coefficient, Partner Risk Allele] <br />
**chr21_causal_snp_idx_100k_exon_example_name.pkl**: a dictionary mapping SNP ID to its mapped Gene Risk <br />
**chr21_causal_genes_100k_exon_example_name.pkl**: a dictionary mapping Gene Feature ID to Gene Risk Score <br />

Histograms of the sampling distributions are created and saved for every major statistical choice.


## Parameter Documentation

| Genotype Parameters | Default Value | Definition |
| ---  | --- | --- |
| -h --help | None | List all parameters |
| -dp --data_path | /DLGWAS/data/ | path to 1000 GP Data |
| -chr | 21 | Chromosome Number |
| -data --data_identifier | 100k | data size identifier |
| -ant --annotation_name | gencode.v19.annotation.gtf | Name of Annotations file for gene/exon mapping |
| -f --features | ["gene", "transcript", "exon"] | List of features for filtering |
| -mtx --matrix_name | genotype.raw | Genotype Matrix (0,1,2) |
| -snplist --snplist_name | genotype.snplist | SNP meta data |


| Phenotype Parameters | Default Value | Definition |
| ---  | --- | --- |
| -h --help | None | List all parameters |
| -dp --data_path | /DLGWAS/data/ | path to data |
| -chr | 21 | Chromosome Number |
| -hd --heritability | 1 | Heritability of phenotype |
| -data --data_identifier | 100k | data size identifier |
| -pname               <br />--phenotype_experiement_name | "" | Name of phenotype simulation experiment |
| -pf --prefilter | "exon" | feature names of genotype filtering |
| -cut --interactive_cut | 0.2 | Fraction of causal SNPs to experience epistatic effects |
| -mask --mask_rate | 0.1 | Fraction of inter-SNP interactions that are masking |
| -df --dosage_frac | 0.5 | Fraction of causal SNPs whose effects are dosage dependent |
| -mic --max_interaction_coeff | 2 | Upper bound for Interaction Coefficient between two SNPs|
| --causal_snp_mode | "gene" | Method to select causal SNPs {gene, random} |
| -num_snps --n_causal_snps | 100 | Number of Causal SNPs <br /> **required for random mode** |
| -cgc --causal_gene_cut | 0.05 | Fraction of Causal Genes <br /> **required for gene mode** |
| -mgr --max_gene_risk | 5 | Upper bound for Gene Risk Coefficient <br /> **required for gene mode** |


## Results
**TODO**
Overview of paper and LINK


## Simulation Playground

[Exploratory Notebook](Example_Simulation_Playground.ipynb) details the custom genotype data creation process for phenotype simulation.

Utilizing randomly generated SNPs, the notebook walks through how to form custom genotype datasets for phenotype simulation. Generated outputs are stored in the Chromosome 0 directory and are used to test the validity of the package.

The command below can be run inside the PACKAGE directoryto create sample data for testing purposes.
```
gwas-sim-public phenotype -dp ./ -chr 0 --data_identifier test --prefilter exon --phenotype_experiement_name playground_example
```