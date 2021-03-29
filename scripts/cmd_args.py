#!/usr/bin/env python

#
# Copyright (c) 2021, NVIDIA CORPORATION.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#

"""Module to parse command line arguments of GWAS Data Simulation"""

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
    parser.add_argument('-chr', dest='chr', required=True, type=int, nargs=None, action = 'store', default=21)
    parser.add_argument('--data_path', '-dp', dest='data_path', required=True, type=str, nargs=None, action = 'store', default="/DLGWAS/data/")
    parser.add_argument('--data_identifier','-data', dest='data_identifier', required=True, type=str, nargs=None, action = 'store', default= "100k")
    
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
    parser.add_argument('--snplist_name', '-snplist', dest='snplist_name', required=False, type=str, nargs=None, action = 'store', default="genotype.snp_list")
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
    parser.add_argument('--prefilter', '-pf', dest='prefilter', type=str, nargs=None, action = 'store', default="exon")
    parser.add_argument('--interactive_cut', '-cut', dest='interactive_cut', required=True, type=float, nargs=None, action = 'store', default=0.2)
    parser.add_argument('--mask_rate', '-mask', dest='mask_rate', required=True, type=float, nargs=None, action = 'store', default=0.1)
    parser.add_argument('--dosage_frac', '-df', dest='dosage_frac', required=True, type=float, nargs=None, action = 'store', default=0.5)
    parser.add_argument('--max_interaction_coeff', '-mic', dest='max_interaction_coeff', required=True, type=float, nargs=None, action = 'store', default=2)
    parser.add('--causal_snp_mode', required=True, type=str, choices=['gene', 'random'],help="causal snp generation method")
    parser.add_argument('--n_causal_snps', '-num_snps', dest='n_causal_snps', type=int, nargs=None, action = 'store', default=100)
    parser.add_argument('--causal_gene_cut', '-cgc', dest='causal_gene_cut', type=float, nargs=None, action = 'store', default=0.05)
    
    parser.add_argument('--noise_scalar', '-noise', dest='noise_scalar', type=float, nargs=None, action = 'store', default=1)
    parser.add_argument('--heritability', '-hrd', dest='heritability', type=float, nargs=None, action = 'store', default=1)
    parser.add_argument('--phenotype_threshold', '-pthresh', dest='phenotype_threshold', type=float, nargs=None, action = 'store', default=50)
    parser.add_argument('--max_gene_risk', '-mgr', dest='max_gene_risk', type=float, nargs=None, action = 'store', default=5)
    parser.add('--config', required=False, is_config_file=True, help='config file path')
#     if arguments are required will it error if they are specified in the config self.phenotype_threshold

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
    # model args
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
#     should we have nested subparsers?
    parser_pheno = subparsers.add_parser(
        'phenotype',
        config_file_parser_class=configargparse.YAMLConfigFileParser,
        default_config_files=[pheno_config_path])
    add_phenotype_options(parser_pheno)
    # =========================================================================
    args, extra = parser.parse_known_args()

    return args
