# GEPSi: GWAS Epistatic Phenotype Simulator

GEPSi is a toolkit to simulate phenotypes for GWAS analysis, given input genotype data for a population.

## Installation

### System requirements

* Python 3.6+

## Build from Source

### 1. Clone repository

#### Latest released version
This will clone the repo to the `main` branch, which contains code for latest released version
and hot-fixes.

```
git clone --recursive -b master https://github.com/clara-genomics/GEPSi.git
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

Genotype data should be supplied in a `.raw` format along with a `.bim` snplist file. GEPS gives us the ability to format the genotype data matrix and associated annotations into an annotated csv file.
    
```
    gepsi genotype -data_path /GWAS/data/chr21/ --matrix_name genotype.raw --snplist_name full_snplist.bim
```
    
Results in the creation of a `.h5` file containing a Person X SNP matrix with Genotype Values of 0,1,2 and and annotated snplist `.csv` that is needed to run the phenotype simulation. The snplist has columns for Chromosome, Feature ID, Position, Allele 1, Allele 2, and Risk Allele. 

The `.raw` and `.bim` files can be produced from other formats using [PLINK](https://www.cog-genomics.org/plink/). PLINK can also be used to filter SNPs within selected regions (exons, transcripts, or genes) as well as filter SNPs based on their allele frequencies. 

For example, we used the following PLINK v1.9 command to filter and format genotype data for human chromosome 21:

```
/plink \
  --gen gensim_chr21_100k.controls.gen.gz \
  --sample gensim_chr21_100k.sample \
  --maf 0.01 \
  --extract range <BED file containing exon positions for chr21> \
  --allow-no-sex \
  --snps-only \
  --recode A \
  --oxford-single-chr 21 \
  --out genotype
  
/plink \
  --gen gensim_chr21_100k.controls.gen.gz \
  --sample gensim_chr21_100k.sample \
  --maf 0.01 \
  --extract range <BED file containing exon positions for chr21> \
  --allow-no-sex \
  --snps-only \
  --oxford-single-chr 21 \
  --make-just-bim \
  --out full_snplist
    
```
Resulting in the creation of

    /GWAS/data/genotype.raw: a Person X SNP Genotype Matrix
    /GWAS/data/full_snplist.bim: Meta data for each SNP


### 2. Generating Phenotypes

Create Phenotypes for generated phenotypes using default values.
```
gepsi phenotype --data_path /GWAS/data/chr21/ --data_identifier chr21_100k --prefilter exon --phenotype_experiement_name example_name
```
    
Results in the creation of

    /DLGWAS/data/chr21/phenotype_chr21_100k_exon_example_name.pkl
    /DLGWAS/data/chr21/effect_size_chr21_100k_exon_example_name.pkl
    /DLGWAS/data/chr21/interactive_snps_chr21_100k_exon_example_name.pkl
    /DLGWAS/data/chr21/causal_snp_idx_chr21_100k_exon_example_name.pkl
    /DLGWAS/data/chr21/causal_genes_chr21_100k_exon_example_name.pkl
    
**phenotype_chr21_100k_exon_example_name.pkl**: a list of binary phenotypes for each person defined by the Genotype Matrix <br />
**effect_size_chr21_100k_exon_example_name.pkl**: a dictionary with key SNP index and value a list of the genotype indexed effect sizes  <br />
**interactive_snps_chr21_100k_exon_example_name.pkl**: a dictionary that maps causal snp indices to a list of length 3 [Interactive SNP Index Pair, Interaction Coefficient, Partner Risk Allele] <br />
**causal_snp_idx_chr21_100k_exon_example_name.pkl**: a dictionary mapping SNP ID to its mapped Gene Risk <br />
**causal_genes_chr21_100k_exon_example_name.pkl**: a dictionary mapping the causal Gene Feature IDs to Gene Risk Scores <br />

Histograms of the sampling distributions are created and saved for every major statistical product.


## Parameter Documentation

| Genotype Parameters | Default Value | Definition |
| ---  | --- | --- |
| -h --help | None | List all parameters |
| -dp --data_path | /GWAS/data/ | path to 1000 GP Data |
| -data --data_identifier | chr1_100k | genotype file name identifier |
| -ant --annotation_name | gencode.v19.annotation.gtf | Name of Annotations file for gene/exon mapping |
| -f --features | ["gene", "transcript", "exon"] | List of features for filtering |
| -rr --risk_rare | False | Use the rare allele as the risk allele |
| -sep --separator | \t | Genetic file separator |
| -ign_map --ignore_gene_map | False | Skip Gene Mapping |
| -low_mem --memory_cautious | False | Use batched reading of Matrix raw file |
| -chunk --matrix_chunk_size | 1000 | Chunk size for low memory matrix read |
| -mtx --matrix_name | genotype.raw | Genotype Matrix (0,1,2) |
| -snplist --snplist_name | genotype.snplist | SNP meta data |


| Phenotype Parameters | Default Value | Definition |
| ---  | --- | --- |
| -h --help | None | List all parameters |
| -dp --data_path | /GWAS/data/ | path to data |
| -hd --heritability | 1 | Heritability of phenotype |
| -data --data_identifier | chr1_100k | genotype file name identifier |
| -pname               <br />--phenotype_experiement_name | "" | Name of phenotype simulation experiment |
| -cut --interactive_cut | 0.2 | Fraction of causal SNPs to experience epistatic effects |
| -mask --mask_rate | 0.1 | Fraction of inter-SNP interactions that are masking |
| -df --dominance_frac | 0.1 | Fraction of causal SNPs whose effects are dominant |
| -rf --recessive_frac | 0.1 | Fraction of causal SNPs whose effects are recessive |
| -mic --max_interaction_coeff | 2 | Upper bound for Interaction Coefficient between two SNPs|
| -pthresh --phenotype_threshold | 50 | Percentile for Phenotype case/control determination |
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

The command below can be run inside the GEPS directoryto create sample data for testing purposes.
```
gepsi phenotype -dp ./sample_data/ --data_identifier chr0_test --phenotype_experiement_name playground_example
```