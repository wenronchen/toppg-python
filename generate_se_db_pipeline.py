#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:29:05 2019

@author: wenrchen
"""

import sys
import merge_psi_coordinate
import match_trans_exon
import generate_se_db
import argparse

#df_SE_name=sys.argv[1]##"../data/SW480/SW480_as_bam/"+se_type+".MATS.JCEC.txt"
#gtf_SE_name=sys.argv[2]##"../data/SW480/SW480_as_bam/fromGTF."+se_type+".txt"
#dbname=sys.argv[3]##"../data/gencode.v28.basic.annotation.gff3"
#se_type=sys.argv[4]##"SE","RI","A3SS","A5SS","MXE"

def main():
    parser = argparse.ArgumentParser(description='Generate the protein database with splicing events. ',formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-f1','--file1', help='The output file from rMATS,for example,SE.MATS.JCEC.txt', required=True)
    parser.add_argument('-f2','--file2', help='The output file from rMATS, for example, fromGTF.SE.txt', required=True)
    parser.add_argument('-g','--gff',help='The annotation file in gff3 format',required=True)
    parser.add_argument('-t','--type',help='The type of splicing events',required=True)
    parser.add_argument('-d','--data',help='The name of data set',required=True)
    
    args = parser.parse_args()
    
    if len(sys.argv) <6 :
        parser.print_help()
    else:
        merge_result=merge_psi_coordinate.merge_file(args.type,args.file1,args.file2)
        match_result=match_trans_exon.match_exon(merge_result,args.gff,args.type)
        generate_se_db.generate_db(args.type,args.gff,match_result,args.data)
                


if(__name__ == "__main__"):
    main()
    