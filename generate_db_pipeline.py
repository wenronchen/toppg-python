#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 10:26:56 2019

@author: wenrchen
"""

import pandas as pd
import get_protein_id 
import export_to_xml
import generate_db_with_xml
import get_new_fsSequence
from Bio import SeqIO
import sys
import argparse

#dfname="../data/ex1/ex1.exonic_variant_function"
#dbname_annotation1="../data/GCF_000001405.38_GRCh38.p12_rna.gbff"
#dbname_annotation2="../data/gencode.v28.com.annotation.gff3"
#dbname_sequence="../data/gencode.v28.swissprot.fasta"

#dfname=sys.argv[1]
#dbname_annotation=sys.argv[2]
#dbname_sequence=sys.argv[3]
#rna_db="../data/gencode.v28.pc_transcripts.fa"

def get_database_type(dbname_sequence):
    version=dbname_sequence
    ver=version.find("v")
    version=version[ver:]
    ver=version.find(".")
    version=version[ver+1:]
    ver=version.find(".")
    version=version[:ver]
    
    return version



def pipeline_point_mutation(dfname,dbname_annotation,dbname_sequence,dataset_name):
    database_type=get_database_type(dbname_sequence)

    protein_map_result=get_protein_id.get_protein_id_nonsynonymous(dfname,dbname_annotation,dataset_name)
    protein_map_result.to_csv('../data/'+dataset_name+'/'+database_type+'/'+dataset_name+'_'+database_type+'_map_result.tsv',sep='\t',index=False)
    
    
    export_to_xml.generate_xml(protein_map_result,dataset_name,database_type)
    
    xml_name='../data/'+dataset_name+'/'+database_type+'/'+dataset_name+'_'+database_type+'_protein_mutation.xml'
    
    generate_db_with_xml.customized_db_generate(xml_name,dbname_sequence,dataset_name)


def pipeline_nonfs_fs(dfname,dbname_annotation,dbname_sequence,rna_db,dataset_name):
    database_type=get_database_type(dbname_sequence)
    
    protein_map_result=get_protein_id.get_protein_id_nonfs(dfname,dbname_annotation,dataset_name)
    protein_map_result.to_csv('../data/'+dataset_name+'/'+database_type+'/'+dataset_name+'_'+database_type+'_nonfs_map_result.tsv',sep='\t',index=False)
    
    export_to_xml.generate_xml_nonfs(protein_map_result,dataset_name,database_type)
    
    xml_name='../data/'+dataset_name+'/'+database_type+'/'+dataset_name+'_'+database_type+'_nonfs_protein_mutation.xml'
    rna_map_protein=get_protein_id.prepare_gencode_database(dbname_annotation)
    
    #rna_map_protein=pd.read_csv('../data/protein_gencode_db.tsv',sep='\t')
    nonfs_seqs=generate_db_with_xml.nonfs_db_generate(xml_name,dbname_sequence,dataset_name)
    fs_seqs=get_new_fsSequence.get_protein_sequence_fs(rna_db,rna_map_protein,dfname,dbname_sequence, dataset_name)
    
    seqs=nonfs_seqs+fs_seqs
    print("The number of sequences in total is "+ str(len(seqs)))
    handle=open("../data/"+dataset_name+"/"+database_type+'/'+dataset_name+"_nonfs_fs_"+database_type+".fasta","w")
    for sequence in seqs:
        SeqIO.write(sequence,handle,"fasta") 

def pipeline_with_combination(dfname,dbname_annotation,dbname_sequence,rna_db,dataset_name):
    database_type=get_database_type(dbname_sequence)
    
    protein_map_result=get_protein_id.get_protein_id_nonsynonymous(dfname,dbname_annotation,dataset_name)
    protein_nonfs_map_result=get_protein_id.get_protein_id_nonfs(dfname,dbname_annotation,dataset_name)
    
    export_to_xml.generate_xml(protein_map_result,dataset_name,database_type)
    export_to_xml.generate_xml_nonfs(protein_nonfs_map_result,dataset_name,database_type)
    
    xml_name='../data/'+dataset_name+'/'+database_type+'/'+dataset_name+'_'+database_type+'_protein_mutation.xml'
    nonfs_xml_name='../data/'+dataset_name+'/'+database_type+'/'+dataset_name+'_'+database_type+'_nonfs_protein_mutation.xml'
    
    no_fs_seqs=generate_db_with_xml.customized_db_generate_with_combination(xml_name,nonfs_xml_name,dbname_sequence,dataset_name)
    
    rna_map_protein=get_protein_id.prepare_gencode_database(dbname_annotation)
    
    fs_seqs=get_new_fsSequence.get_protein_sequence_fs_with_combination(rna_db,rna_map_protein,dfname,dbname_sequence,dataset_name,xml_name,nonfs_xml_name)
    seqs=no_fs_seqs+fs_seqs
    
    print("The number of sequences in total is "+ str(len(seqs)))
    handle=open("../data/"+dataset_name+"/"+database_type+"/"+dataset_name+"_combination_"+database_type+".fasta","w")
    for sequence in seqs:
        SeqIO.write(sequence,handle,"fasta") 

def pipeline_with_combination_only(dfname,dbname_annotation,dbname_sequence,rna_db,dataset_name):
    database_type=get_database_type(dbname_sequence)
    
    protein_map_result=get_protein_id.get_protein_id_nonsynonymous(dfname,dbname_annotation,dataset_name)
    protein_nonfs_map_result=get_protein_id.get_protein_id_nonfs(dfname,dbname_annotation,dataset_name)
    
    export_to_xml.generate_xml(protein_map_result,dataset_name,database_type)
    export_to_xml.generate_xml_nonfs(protein_nonfs_map_result,dataset_name,database_type)
    
    xml_name='../data/'+dataset_name+'/'+database_type+'/'+dataset_name+'_'+database_type+'_protein_mutation.xml'
    nonfs_xml_name='../data/'+dataset_name+'/'+database_type+'/'+dataset_name+'_'+database_type+'_nonfs_protein_mutation.xml'
    
    no_fs_seqs=generate_db_with_xml.customized_db_generate_with_combination_only(xml_name,nonfs_xml_name,dbname_sequence,dataset_name)
    
    rna_map_protein=get_protein_id.prepare_gencode_database(dbname_annotation)
    
    fs_seqs=get_new_fsSequence.get_protein_sequence_fs_with_combination_only(rna_db,rna_map_protein,dfname,dbname_sequence,dataset_name,xml_name,nonfs_xml_name)
    seqs=no_fs_seqs+fs_seqs
    
    print("The number of sequences in total is "+ str(len(seqs)))
    handle=open("../data/"+dataset_name+"/"+database_type+"/"+dataset_name+"_combination_only_"+database_type+".fasta","w")
    for sequence in seqs:
        SeqIO.write(sequence,handle,"fasta") 

def main():
    parser = argparse.ArgumentParser(description='Generate the protein database with mutations. ',formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-f','--file', help='The output file from Annovar,for example,ex1.exonic_variant_function', required=True)
    parser.add_argument('-g','--gff',help='The annotation file in gff3 format',required=True)
    parser.add_argument('-s','--sequence', help='The fasta file of protein sequences', required=True)
    parser.add_argument('-t','--type',help='The type of mutations when generating database, \
                        either f(frameshift and non-frameshift) or p(point mutations) or c(with combination)',required=True)##f and p and c and o(combination only))
    parser.add_argument('-r','--rna', help='The fasta file of transcript sequences', required=False)
    args = parser.parse_args()
    
    if len(sys.argv) <5 :
        parser.print_help()
    else:
        dataset_name=''
        df_name=args.file
        while(df_name.find('/')!=-1):
            label=df_name.find('/')
            df_name=df_name[label+1:]
        label1=df_name.find('.')
        label2=df_name.find('_')
        if(label1<label2):
            dataset_name=df_name[:label1]
        else:
            dataset_name=df_name[:label2]
        if(args.type=='p'):
            pipeline_point_mutation(args.file,args.gff,args.sequence,dataset_name)
        if(args.type=='f'):
            pipeline_nonfs_fs(args.file,args.gff,args.sequence,args.rna,dataset_name)
        if(args.type=='c'):
            pipeline_with_combination(args.file,args.gff,args.sequence,args.rna,dataset_name)
        if(args.type=='o'):
            pipeline_with_combination_only(args.file,args.gff,args.sequence,args.rna,dataset_name)

if(__name__ == "__main__"):
    main()