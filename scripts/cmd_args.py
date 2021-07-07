#
# Copyright (c) 2021, NVIDIA CORPORATION & AFFILIATES.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#

"""Module to parse command line arguments of GEPSi"""

import os
import configargparse


def add_common_options(parser):
    """Add common options to the parser.
    Args:
        parser : parser to add the arguments to.
    Return:
        parser : After adding the arguments.
    """
    # Common Arguments
    parser.add_argument('--data_path', '-dp', dest='data_path', required=False, type=str, nargs=None, action = 'store', default="/GWAS/data/")
    parser.add_argument('--data_identifier','-data', dest='data_identifier', required=False, type=str, nargs=None, action = 'store', default= "file_name")
    
def add_annotation_options(parser):
    """Add model options to the parser.
    Args:
        parser : parser to add the arguments to.
    Return:
        parser : After adding the arguments.
    """
    add_common_options(parser)
    parser.add_argument('--annotation_name', '-an', dest='annotation_name', required=False, type=str, nargs=None, action = 'store', default="gencode.v19.annotation.gtf")
    parser.add_argument('--features', '-f', dest='features', required=False, type=str, nargs="+", action = 'store', default=["gene", "transcript", "exon"])
    parser.add('--config', required=False, is_config_file=True, help='config file path')

def add_genotype_options(parser):
    """Add genotype options to the parser.
    Args:
        parser : parser to add the arguments to.
    Return:
        parser : After adding the arguments.
    """
    add_common_options(parser)
    parser.add_argument('--annotation_name', '-ant', dest='annotation_name', required=False, type=str, nargs=None, action = 'store', default="gencode.v19.annotation.gtf")
    parser.add_argument('--features', '-f', dest='features', required=False, type=str, nargs="+", action = 'store', default=["gene", "transcript", "exon"])
    parser.add_argument('--matrix_name', '-mtx', dest='matrix_name', required=False, type=str, nargs=None, action = 'store', default="genotype.raw")
    parser.add_argument("--risk_rare",'-rr', default=False, action="store_true", help="Flag to use rare allele as risk allele")
    parser.add_argument("--ignore_gene_map",'-ign_map', default=False, action="store_true", help="Flag to ignore the creation of a gene map")
    parser.add_argument('--snplist_name', '-snplist', dest='snplist_name', required=False, type=str, nargs=None, action = 'store', default="genotype.snp_list")
    parser.add_argument('--memory_cautious', '-low_mem', dest='memory_cautious', default=False, action="store_true", help="Flag to use a batched readin of genotype matrix")
    parser.add_argument('--separator', '-sep', dest='separator', required=False, type=str, nargs=None, action = 'store', default="\t")
    parser.add_argument('--matrix_chunk_size', '-chunk', dest='matrix_chunk_size', type=int, nargs=None, action = 'store', default=1000)        
    parser.add('--config', required=False, is_config_file=True, help='config file path')

def add_phenotype_options(parser):
    """Add phenotype options to the parser.
    Args:
        parser : parser to add the arguments to.
    Return:
        parser : After adding the arguments.
    """
    add_common_options(parser)
    parser.add_argument('--phenotype_experiement_name', '-pname', dest='phenotype_experiement_name', type=str, nargs=None, action = 'store', default="")
    parser.add_argument('--interactive_cut', '-cut', dest='interactive_cut', required=False, type=float, nargs=None, action = 'store', default=0.2)
    parser.add_argument('--mask_rate', '-mask', dest='mask_rate', required=False, type=float, nargs=None, action = 'store', default=0.1)
    parser.add_argument('--dominance_frac', '-df', dest='dominance_frac', required=False, type=float, nargs=None, action = 'store', default=0.1)
    parser.add_argument('--recessive_frac', '-rf', dest='recessive_frac', required=False, type=float, nargs=None, action = 'store', default=0.1)
    parser.add_argument('--max_interaction_coeff', '-mic', dest='max_interaction_coeff', required=False, type=float, nargs=None, action = 'store', default=2)
    parser.add('--causal_snp_mode', required=False, type=str, choices=['gene', 'random'],help="causal snp generation method")
    parser.add_argument('--n_causal_snps', '-num_snps', dest='n_causal_snps', type=int, nargs=None, action = 'store', default=100)
    parser.add_argument('--causal_gene_cut', '-cgc', dest='causal_gene_cut', type=float, nargs=None, action = 'store', default=0.05)
    parser.add_argument('--heritability', '-hrd', dest='heritability', type=float, nargs=None, action = 'store', default=1)
    parser.add_argument('--stratify', '-st', dest='stratify', default=False, action="store_true", help="Flag to stratify individuals based on given groups. Group and Group coefficient files must be included in --data_path.")
    parser.add_argument('--phenotype_threshold', '-pthresh', dest='phenotype_threshold', type=float, nargs=None, action = 'store', default=50)
    parser.add_argument('--max_gene_risk', '-mgr', dest='max_gene_risk', type=float, nargs=None, action = 'store', default=5)
    parser.add('--config', required=False, is_config_file=True, help='config file path')

def parse_args(root_dir):
    """Parse command line arguments.
    Args:
        root_dir : Path to the root directory,
        where the configs folder can be found.
    Return:
        args : parsed argument object.
    """
    parser = configargparse.ArgParser()
    subparsers = parser.add_subparsers(dest="mode")
    # =========================================================================
    # annotation args
    annotation_config_path = os.path.join(root_dir, 'configs', 'genotype.yaml')
    parser_annotation = subparsers.add_parser(
        'annotation',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        default_config_files=[annotation_config_path])
    add_annotation_options(parser_annotation)
    # =========================================================================
    # genotype args
    geno_config_path = os.path.join(root_dir, 'configs', 'genotype.yaml')
    parser_geno = subparsers.add_parser(
        'genotype',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        default_config_files=[geno_config_path])
    add_genotype_options(parser_geno)
    # =========================================================================
    # phenotype args
    pheno_config_path = os.path.join(root_dir, 'configs', 'phenotype.yaml')
    parser_pheno = subparsers.add_parser(
        'phenotype',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        default_config_files=[pheno_config_path])
    add_phenotype_options(parser_pheno)
    # =========================================================================
    args, extra = parser.parse_known_args()
    if args.mode == "genotype":
        # check genotype params
        pass
    elif args.mode == "phenotype":
        params = {
            "interactive_cut":args.interactive_cut,
            "mask_rate":args.mask_rate,
            "dominance_frac":args.dominance_frac,
            "recessive_frac":args.recessive_frac,
            "dosage_frac": 1 - args.dominance_frac - args.recessive_frac,
            "heritability":args.heritability }
        validate_parameters(params, case = 0)
        if args.causal_snp_mode == "random":
             params = {
                "max_interaction_coeff":args.max_interaction_coeff,
                "phenotype_threshold":args.phenotype_threshold,
                "n_causal_snps":args.n_causal_snps}
        if args.causal_snp_mode == "gene":
            params = {
                "max_interaction_coeff":args.max_interaction_coeff,
                "phenotype_threshold":args.phenotype_threshold,
                "causal_gene_cut":args.causal_gene_cut,
                "max_gene_risk":args.max_gene_risk}
        validate_parameters(params, case = 1)
    elif args.mode == "annotation":
        # check annotation params
        pass
    elif not args.mode:
        parser.print_help()
        parser.error("Please choose simulation mode")
    return args


def validate_parameters(params, case):
    if case == 0:
        for name, param in params.items():
            if param < 0 or param > 1:
                raise ValueError('Invalid value for {}. Acceptible values [0,1]'.format(name))
    elif case == 1:
        for name, param in params.items():
            if param <= 0:
                raise ValueError('Invalid value for {}. Acceptible values > 0'.format(name)) 
