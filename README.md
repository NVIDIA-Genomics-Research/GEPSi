# GWAS Data Simulation

PACKAGE is a toolkit for the simulation of diverse genotypes and phenotypes for GWAS analysis.

## Installation

### System requirements

* **TODO: Not sure yet**

### PyPI installation
**TODO: To install atacworks in your environment, run the following in your terminal**
```
pip install gwas-sim
```

### Docker Installation
**TODO: Docker based installation**

If you'd like to skip all installation and use a pre-installed docker image instead,
follow the instructions [here](Dockerfile.md), section "Pre-installed AtacWorks".

**We have the docker file that clones the repo**

## Build from Source
Follow the instructions below if you are interested in running gwas-sim outside
of docker.

### 1. Clone repository

#### Latest released version
This will clone the repo to the `main` branch, which contains code for latest released version
and hot-fixes.

```
git clone --recursive -b master https://github.com/clara-genomics/gwas-data-simulation-public.git
```

### 2. Install dependencies

* Install gwas-sim and its associated dependencies from requirements.txt

    ```
    pip install .
    ```
    
### 3. Tests

Run unit tests to verify that installation was successful

    ```
    python -m pytest tests/
    ```
    
## Data
Genomewide Data is found here: [1000 Genome Project](https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html) (Note: in older format. HAPGEN2 does not work directly with .vcf input)

[HAPGEN2](https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html) supplied with 1000 Genome Project data is able to produce variable number of genetic diverse people that will later be used to determine casual SNPs and phenotype. HAPGEN2 was chosen due to its understandable sampling logic as well as focus on retention of ld patterns in input data. It allows us to create populations of any size with SNP allele frequencies found in true data. 

[PLINK](https://www.cog-genomics.org/plink/) is used to filter the genotypes by any combination of features (exon, transcript, or gene) as well as covert to the Minor Allele frequency matrix. This can give as a matrix of people X SNP with the value being 0,1,2 correspodning to the MAF.

The above features are taken from [Gencode Genome Annotations](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz). This gives our mapping of snp postion to features and allows us to annotate our SNP data.


## Workflow
GWAS Data simulation relies on Hapgen2 to simulate a diverse genetic population of the desired size. It then simulates the phenotypes for this population for downstream analysis.

### 1. Generating Genotypes

Example of how we generated 100000 Individuals from the Chromosome 21 using the 1000 Genome Project.
**TODO: Update CLI to handle pre and post protected binaries**

Utilize HAPGEN2 to create diverse genotype Controls to be used for the Phenotype simulation of desired population size. Please see HAPGEN2 Documentation for detailed parameter walkthough.

See the Exploratory Notebook for an example of the data manipulations for custom input genotypes.
```
/hapgen2 \
       -h /DLGWAS/data/1000GP_Phase3/1000GP_Phase3_chr21.hap \
       -l /DLGWAS/data/1000GP_Phase3/1000GP_Phase3_chr21.legend \
       -m /DLGWAS/data/1000GP_Phase3/genetic_map_chr21_combined_b37.txt \
       -n 100000 100 \
       -dl 100011 1 1 1 \
       -o /DLGWAS/data/gensim_chr21_100k.gz \
       -no_haps_output
```
Results in the creation of

    /DLGWAS/data/gensim_chr21_100k.controls.gen.gz
    /DLGWAS/data/gensim_chr21_100k.sample

Create the filtered annotation via gwas-sim genotype --data_path /DLGWAS/data/ --annotation_name annotations -chr 21 - features exon

    /DLGWAS/data/21_annotations_exon.csv

Utlizie PLINK to filter and format the data into tabular form for down stream compuational analysis.
```
cut -f2,4- -d" " ${data_path}21_annotations_${feature_names}.csv | tail -n +2 | sed 's/^chr//' > ${data_path}${feature_names}_ranges_chr21.txt
# PLINK v1.9
/plink \
    --gen ${data_path}gensim_chr21_100k.controls.gen.gz \
    --sample ${data_path}gensim_chr21_100k.sample \
	--maf 0.01 \
	--extract range ${data_path}${feature_names}_ranges_chr21.txt \
	--allow-no-sex \
	--snps-only \
	--write-snplist \
	--recode A \
	--oxford-single-chr 21 \
    --out ${data_path}genotype
    
sed -i 's/:/ /g' ${data_path}genotype.snplist
# Cleanup snplist formatting for downstream usage
```
Results in the creation of

    /DLGWAS/data/genotype.raw a Person X SNP Genotype Matrix
    /DLGWAS/data/genotype.snplist Meta data for each SNP


PACKAGE gives us the ability to format the data matrix and associated annotations into a gene mapped annotated csv file.
    
```
    gwas-sim genotype -chr21 -data_path /DLGWAS/data/ --matrix_name genotype.raw --snplist_name genotype.snplist
```
    
Results in the creation of an annotated matrix. Users are encouraged to use PLINK for their data conversions if not using 1000 Genome Project. PACKAGE provides easy annotation cleaning and genotype data manipulation. Phenotype simulation requires custom pandas dataframe so as long as other formats can be read in via PLINK phenotype simulation is possible.

    /DLGWAS/data/chr21_genotype_100k_exon.csv

chr21_genotype_100k_exon.csv an annotated SNP X Person matrix with Genotype Values of 0,1,2
**This is what is needed to run the phenotype simulation**
![Alt text](./annotated_matrix_result.png?raw=true "Title")

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
    
**chr21_phenotype_100k_exon_example_name.pkl** a list of binary phenotypes for each person defined by the Genotype Matrix <br />
**chr21_effect_size_100k_exon_example_name.pkl** is a dictionary with key SNP index and value a list of the dosage dependent effect sizes  <br />
**chr21_interactive_snps_100k_exon_example_name.pkl** is a dictionary that maps causal snp indexes to a length 3 list [Interactive SNP Index Pair, Interaction Coefficeint, Partner Risk Allele] <br />
**chr21_causal_snp_idx_100k_exon_example_name.pkl** is a dictionary mapping SNP ID to its mapped Gene Risk <br />
**chr21_causal_genes_100k_exon_example_name.pkl** is a dictionary mapping Gene Feature ID to Gene Risk Score <br />

### Results
**TODO**
Here I would like to showcase a figure or 2 from the paper or just link to paper as it is a good read.

### Parameter Documentation

| Annotation Parameters | Default Value | Definition |
| ---  | --- | --- |
| -h --help | None | List all parameters |
| -dp --data_path | /DLGWAS/data/ | path to 1000 GP Data |
| -chr | 21 | Chromosome Number |
| -data --data_identifier | 100k | data size identifier |
| -ant --annotation_name | gencode.v19.annotation.gtf | Name of Annotations file for gene/exon mapping |
| -f --features | ["gene", "transcript", "exon"] | List of features for HAPGEN2 filtering |

| Genotype Parameters | Default Value | Definition |
| ---  | --- | --- |
| -h --help | None | List all parameters |
| -dp --data_path | /DLGWAS/data/ | path to 1000 GP Data |
| -chr | 21 | Chromosome Number |
| -data --data_identifier | 100k | data size identifier |
| -ant --annotation_name | gencode.v19.annotation.gtf | Name of Annotations file for gene/exon mapping |
| -f --features | ["gene", "transcript", "exon"] | List of features for HAPGEN2 filtering |
| -mtx --matrix_name | genotype.raw | PLINK Recode A Genotype Matrix (0,1,2) |
| -snplist --snplist_name | genotype.snplist | PLINK Recode A SNP meta data |


| Phenotype Parameters | Default Value | Definition |
| ---  | --- | --- |
| -h --help | None | List all parameters |
| -dp --data_path | /DLGWAS/data/ | path to 1000 GP Data |
| -chr | 21 | Chromosome Number |
| -hd --heredity | 1 | Heredity coeficient |
| -data --data_identifier | 100k | data size identifier |
| -pname               <br />--phenotype_experiement_name | "" | Name of phenotype simulation experiement |
| -pf --prefilter | "exon" | feature names of genotype filtering |
| -cut --interactive_cut | 0.2 | Fraction of Causal SNPS to experience epistatic effects |
| -mask --mask_rate | 0.1 | Fraction Interaction relations that are masking |
| -df --dosage_frac | 0.5 | Fraction of Causal SNPS to experience dosage dependent effect size |
| -mic --max_interaction_coeff | 2 | Upper bound for Interaction Coefficient Uniform Distribution |
| --causal_snp_mode | "gene" | Phenotype Algorithm Selector {gene, random} |
| -num_snps --n_causal_snps | 100 | Number of Causal SNPs <br /> **required for random mode** |
| -cgc --causal_gene_cut | 0.05 | Fraction of Causal Genes <br /> **required for gene mode** |
| -mgr --max_gene_risk | 5 | Upper bound for Gene Risk Coefficient Uniform Distribution <br /> **required for gene mode** |

