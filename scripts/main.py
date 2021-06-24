#
# Copyright (c) 2021, NVIDIA CORPORATION & AFFILIATES.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#
import os, sys
from scripts.phenotype_simulator import PhenotypeSimulator
from scripts.genotype_simulator import GenotypeSimulator
from scripts.cmd_args import parse_args

def main():
    """Main Method. Dispatch args to desired Simulation Step."""
    root_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), ".."))
    args = parse_args(root_dir)
    print("Main: ", args)
    if args.mode == "genotype":
        genotype_sim = GenotypeSimulator(args)
        data = genotype_sim.simulate_genotype()
    elif args.mode == "phenotype":
        phenotype_sim = PhenotypeSimulator(args)
        phenotype, causal_snps_idx, effect_size, interactive_snps = phenotype_sim.simulate_phenotype()
    elif args.mode == "annotation":
        genotype_sim = GenotypeSimulator(args)
        genotype_sim.clean_annotations(args.features)
        
if __name__ == '__main__':
    main()
