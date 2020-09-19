#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:29:05 2019

@author: wenrchen
"""

import merge_psi_coordinate
import generate_all_se_db
import generate_all_mutation_db
import generate_all_mutation_combination_db
import generate_all_se_combination_db
import argparse


def main():
    parser = argparse.ArgumentParser(description='Generate the protein database with variations. ',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-s','--splicing', help='If the parameter is added, the splicing variations will be added to the customized protein database', action='store_true')
    parser.add_argument('-f','--file', help='The output file from Annovar,for example,ex1.exonic_variant_function', required=False)
    parser.add_argument('-f1','--file1', help='The output file from rMATS,for example,SE.MATS.JCEC.txt', required=False)
    parser.add_argument('-f2','--file2', help='The output file from rMATS, for example, fromGTF.SE.txt', required=False)
    parser.add_argument('-g','--gff',help='The annotation file in gff3 format',required=True)
    parser.add_argument('-r','--rna', help='The fasta file of transcript sequences', required=True)
    parser.add_argument('-o','--output',help='The name of output file', default="customized_db",required=False)
    parser.add_argument('-t','--het', help='The number of heterozygous genetic variants', default="0",required=True)
    parser.add_argument('-e','--exclude', help='Exclude the original protein sequence', action='store_true')
    
    args = parser.parse_args()
    het=int(args.het)
    if(args.splicing==True):
        if(args.file1=="" or args.file2=="" or args.file=="" or args.gff=="" or args.rna==""):
            print("Some parameters required for this pattern are missing!")
        else:
            merge_result=merge_psi_coordinate.merge_file(args.file1,args.file2)
            if(het==0):
                generate_all_se_db.generate_se_sequence(merge_result,args.file,args.gff,args.rna,args.exclude,args.output)
            elif(het==1 or het==2):
                generate_all_se_combination_db.generate_se_sequence_combination(merge_result,args.file,args.gff,args.rna,het,args.exclude,args.output)
                

        
    else:   
        
        if((args.file=="") or (args.gff=="") or (args.rna=="")):
            print("Some parameters required for this pattern are missing!")
        else:
            if((het==1) or (het==0)):
                generate_all_mutation_db.get_new_sequence(args.file,args.gff,args.rna,het,args.exclude,args.output)
            elif(het==2):
                generate_all_mutation_combination_db.get_new_sequence_combination(args.file,args.gff,args.rna,args.exclude,args.output)
                
    
        

if(__name__ == "__main__"):
    main()
    
