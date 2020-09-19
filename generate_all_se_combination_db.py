#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: wenrchen
"""

import sys
import pandas as pd
import time
#from gffutils.iterators import DataIterator
from Bio import SeqIO
#from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import pybedtools
import generate_all_se_db
import generate_all_mutation_db

import warnings
warnings.filterwarnings("ignore")

argv=sys.argv
begin=time.time()

fasta = pybedtools.example_filename('GRCh38.primary_assembly.genome.fa')

def get_version(rna_db):
    version=rna_db
    ver=version.find("trans")
    version=version[ver:]
    ver=version.find(".")
    version=version[ver+1:]
    ver=version.find(".")
    version=version[:ver]
    
    return version


def get_variant_pos_from_des(des):
    start_flag=des.find(':')
    end_flag=des.find('-')
    return int(des[start_flag+1:end_flag])

def generate_se_sequence_combination(df_SE,dfname,dbname,rna_db,het,exclude,output_name):
    
    version=get_version(rna_db)
    if(version=='swissprot'):
        db='sp'
    else:
        db=version
     
    trans_seq_dict,trans_coding_dict=generate_all_se_db.get_trans_seq_dict(rna_db)
    print("transcript sequence dict ready")
    print(len(list(trans_seq_dict.items())))

    
    df=pd.read_csv(dfname,sep='\t',header=None)
    change_df=generate_all_mutation_db.extract_transcript_change(df) 
    trans_index_dict={}
    for i in range(0,change_df.shape[0]):
        if(change_df.iloc[i]['mrna'] not in trans_index_dict.keys()):
            trans_index_dict[change_df.iloc[i]['mrna']]=[i]
        else:
            trans_index_dict[change_df.iloc[i]['mrna']].append(i)
   
    trans_records,strand_dict,protein_id_dict=generate_all_se_db.get_trans_records(dbname,trans_seq_dict)
    
    splicing_event_dict=generate_all_se_db.get_splicing_event_dict(df_SE,trans_records)
    print("splicing events dict ready")
        
    trans_exon_dict=generate_all_se_db.get_transcript_exon_dict(trans_records)
    print("transcript exon dict ready")
    
      
    
    
    middle=time.time()
    print("---- %s minutes ----" % ((middle-begin)/60))
    
    trans_db_annotation_dict={}##output
    cnt=0
    trans_records=trans_records.sort_values(by=['exon_id'],axis=0)

    
    for k in splicing_event_dict.keys():
        cnt+=1
        if(cnt%1000==0):
            print(cnt)
        down_id=splicing_event_dict[k][0]
        skip_id=splicing_event_dict[k][1]
        up_id=splicing_event_dict[k][2]
        for key in trans_exon_dict.keys():
            T=trans_exon_dict[key]
            if((down_id in T) and (up_id in T)):
                #print(key)
                down_index=T.index(down_id)
                #print(down_index)
                up_index=T.index(up_id)
                low_index=min(down_index,up_index)
                high_index=max(down_index,up_index)
                
                tmp_exon_coordinate_dict={}
                tmp_exon_position_dict={}
                tmp_pos=0
                for item in range(1,len(T)):
                    ei=generate_all_se_db.find_key_index('exon_id',T[item],trans_records)
                    if(ei==-1):
                        print("the exon id cannot be found")
                    else:
                        tmp_start=trans_records.iloc[ei]['start']
                        tmp_end=trans_records.iloc[ei]['end']
                        tmp_exon_coordinate_dict[T[item]]=(trans_records.iloc[ei]['chr'],tmp_start,tmp_end)
                        
                        tmp_exon_position_dict[item]=(tmp_pos,tmp_pos+tmp_end-tmp_start+1)
                        tmp_pos+=tmp_end-tmp_start+1
                            
                tmp_seq=trans_seq_dict[key]
                coding_start=trans_coding_dict[key][0]
                coding_end=trans_coding_dict[key][1]
                
                het_list=[]
                hom_position_list=[]
                if((key in trans_index_dict.keys()) and (key in strand_dict.keys())):
                    shift=0
                    transcript=tmp_seq
                    hom_des=[]
                    
                    for i in trans_index_dict[key]:
                        if(change_df.iloc[i]['mutation_type'].find('snv')!=-1):
                            if(change_df.iloc[i]['snp_type']=='hom'):
                                hom_position_list.append((int(change_df.iloc[i]['c_start']),change_df.iloc[i]['mutation_type']))
        
                                if(strand_dict[key]=='-'):
                                    tmp=Seq(str(tmp_seq),IUPAC.ambiguous_dna).complement()
                                    tmp=str(tmp)
                                    transcript=generate_all_mutation_db.change_seq(tmp,int(change_df.iloc[i]['c_start'])+shift+coding_start-1,\
                                                          int(change_df.iloc[i]['c_end'])+shift+coding_start-1,str(Seq(str(change_df.iloc[i]['c_content'])).complement()),\
                                                          change_df.iloc[i]['mutation_type'])
                                else:
                                    transcript=generate_all_mutation_db.change_seq(tmp_seq,int(change_df.iloc[i]['c_start'])+shift+coding_start-1,\
                                                          int(change_df.iloc[i]['c_end'])+shift+coding_start-1,change_df.iloc[i]['c_content'],\
                                                          change_df.iloc[i]['mutation_type'])
                                hom_des.append(change_df.iloc[i]['mutation_type']+":"+str(change_df.iloc[i]['c_start'])+'-'+\
                                str(change_df.iloc[i]['c_end'])+str(change_df.iloc[i]['c_content'])+'_')
                                
#                                if(change_df.iloc[i]['mutation_type'].find('del')!=-1):
#                                    shift-=(int(change_df.iloc[i]['c_end'])-int(change_df.iloc[i]['c_start'])+1)
#                                elif(change_df.iloc[i]['mutation_type'].find('ins')!=-1):
#                                    shift+=(int(change_df.iloc[i]['c_end'])-int(change_df.iloc[i]['c_start'])+1)
                            else:
                                het_list.append(i)
                    
                    
                    
                    het_number=len(het_list)
                    het_des=[]
                    het_seqs=[]
                    if(len(hom_position_list)!=0):
                        if(strand_dict[key]=='-'):
                            tmp=Seq(str(transcript),IUPAC.ambiguous_dna).complement()
                            tmp=str(tmp)
                        else:
                            tmp=transcript
                        het_seqs.append(tmp)
                    
                    else:
                        if(strand_dict[key]=='-'):
                            transcript=Seq(str(transcript),IUPAC.ambiguous_dna).complement()
                            transcript=str(transcript)
                        het_seqs.append(transcript)
                    het_des.append("")
                    
                    if(het==1):
                        for n in range(0,het_number):
                            tmp_het_des=[]
                            new_sequence=""
                            splicing_position=0
                            if((high_index-low_index)==1):
                                splicing_position=tmp_exon_position_dict[low_index][1]
                            elif(T[low_index+1]==skip_id):
                                splicing_position=tmp_exon_position_dict[low_index+1][0]
                            if(abs(int(change_df.iloc[het_list[n]]['c_start'])-splicing_position)<=900):
                                if((int(change_df.iloc[het_list[n]]['c_start'])>=coding_start-1) & (int(change_df.iloc[het_list[n]]['c_start'])<=coding_end)):
                                    if(strand_dict[key]=='-'):
                                        new_sequence=generate_all_mutation_db.change_seq(transcript,int(change_df.iloc[het_list[n]]['c_start'])+shift+coding_start-1,\
                                                              int(change_df.iloc[het_list[n]]['c_end'])+shift+coding_start-1,\
                                                              str(Seq(str(change_df.iloc[het_list[n]]['c_content'])).complement()),\
                                                              change_df.iloc[het_list[n]]['mutation_type'])
                                    else:
                                        new_sequence=generate_all_mutation_db.change_seq(transcript,int(change_df.iloc[het_list[n]]['c_start'])+shift+coding_start-1,\
                                                              int(change_df.iloc[het_list[n]]['c_end'])+shift+coding_start-1,change_df.iloc[het_list[n]]['c_content'],\
                                                              change_df.iloc[het_list[n]]['mutation_type'])
                                    tmp_het_des.append(change_df.iloc[het_list[n]]['mutation_type']+":"+\
                                    str(change_df.iloc[het_list[n]]['c_start'])+'-'+str(change_df.iloc[het_list[n]]['c_end'])+\
                                    str(change_df.iloc[het_list[n]]['c_content']))
                                    
    #                                if(change_df.iloc[i]['mutation_type'].find('del')!=-1):
    #                                    shift-=(int(change_df.iloc[i]['c_end'])-int(change_df.iloc[i]['c_start'])+1)
    #                                elif(change_df.iloc[i]['mutation_type'].find('ins')!=-1):
    #                                    shift+=(int(change_df.iloc[i]['c_end'])-int(change_df.iloc[i]['c_start'])+1)
                                het_des.append(tmp_het_des)
                                if(strand_dict[key]=='-'):
                                    new_seq=Seq(str(new_sequence),IUPAC.ambiguous_dna).complement()
                                    new_seq=str(new_seq)
                                else:
                                    new_seq=new_sequence
                                het_seqs.append(new_seq)
                    elif(het==2):
                        for n in range(0,het_number):
                            tmp_het_des=[]
                            new_sequence=""
                            splicing_position=0
                            if((high_index-low_index)==1):
                                splicing_position=tmp_exon_position_dict[low_index][1]
                            elif(T[low_index+1]==skip_id):
                                splicing_position=tmp_exon_position_dict[low_index+1][0]
                            if(abs(int(change_df.iloc[het_list[n]]['c_start'])-splicing_position)<=900):
                                if((int(change_df.iloc[het_list[n]]['c_start'])>=coding_start-1) & (int(change_df.iloc[het_list[n]]['c_start'])<=coding_end)):
                                    if(strand_dict[key]=='-'):
                                        new_sequence=generate_all_mutation_db.change_seq(transcript,int(change_df.iloc[het_list[n]]['c_start'])+shift+coding_start-1,\
                                                              int(change_df.iloc[het_list[n]]['c_end'])+shift+coding_start-1,\
                                                              str(Seq(str(change_df.iloc[het_list[n]]['c_content'])).complement()),\
                                                              change_df.iloc[het_list[n]]['mutation_type'])
                                    else:
                                        new_sequence=generate_all_mutation_db.change_seq(transcript,int(change_df.iloc[het_list[n]]['c_start'])+shift+coding_start-1,\
                                                              int(change_df.iloc[het_list[n]]['c_end'])+shift+coding_start-1,change_df.iloc[het_list[n]]['c_content'],\
                                                              change_df.iloc[het_list[n]]['mutation_type'])
                                    tmp_het_des.append(change_df.iloc[het_list[n]]['mutation_type']+":"+\
                                    str(change_df.iloc[het_list[n]]['c_start'])+'-'+str(change_df.iloc[het_list[n]]['c_end'])+\
                                    str(change_df.iloc[het_list[n]]['c_content']))
                                    
    #                                if(change_df.iloc[i]['mutation_type'].find('del')!=-1):
    #                                    shift-=(int(change_df.iloc[i]['c_end'])-int(change_df.iloc[i]['c_start'])+1)
    #                                elif(change_df.iloc[i]['mutation_type'].find('ins')!=-1):
    #                                    shift+=(int(change_df.iloc[i]['c_end'])-int(change_df.iloc[i]['c_start'])+1)
                                het_des.append(tmp_het_des)
                        
                            for j in range(n,len(het_list)):
                                tmp_new_sequence=new_sequence
                                tmp_new_des=tmp_het_des

                                if(j!=n):
                                    if(abs(int(change_df.iloc[het_list[j]]['c_start'])-splicing_position)<=900):
                                        if((int(change_df.iloc[het_list[j]]['c_start'])>=coding_start-1) & (int(change_df.iloc[het_list[j]]['c_start'])<=coding_end)):
                                            if(strand_dict[key]=='-'):
                                                tmp_new_sequence=generate_all_mutation_db.change_seq(tmp_new_sequence,int(change_df.iloc[het_list[j]]['c_start'])+shift,\
                                                                      int(change_df.iloc[het_list[j]]['c_end'])+shift,\
                                                                      str(Seq(str(change_df.iloc[het_list[j]]['c_content'])).complement()),\
                                                                      change_df.iloc[het_list[j]]['mutation_type'])
                                            else:
                                                tmp_new_sequence=generate_all_mutation_db.change_seq(tmp_new_sequence,int(change_df.iloc[het_list[j]]['c_start'])+shift,\
                                                                      int(change_df.iloc[het_list[j]]['c_end'])+shift,change_df.iloc[het_list[j]]['c_content'],\
                                                                      change_df.iloc[het_list[j]]['mutation_type'])
                                            tmp_new_des=tmp_new_des.append(change_df.iloc[het_list[j]]['mutation_type']+":"+\
                                                str(change_df.iloc[het_list[j]]['c_start'])+'-'+str(change_df.iloc[het_list[j]]['c_end'])+\
                                                str(change_df.iloc[het_list[j]]['c_content']))
                                        het_des.append(tmp_new_des)
                                
                                if(len(tmp_new_sequence)!=0):
                                    if(strand_dict[key]=='+'):
                                        new_seq=tmp_new_sequence
                                    else:
                                        new_seq=str(Seq(str(tmp_new_sequence),IUPAC.ambiguous_dna).complement())
                                    het_seqs.append(new_seq)
                    
                    
                    
                            
                    for c in range(0, len(het_seqs)):##tmp_cnt is the number of modified sequences 
                        tmp_seq=het_seqs[c]
                    
                        if((high_index-low_index)==1):
                            
                            position=tmp_exon_position_dict[low_index][1]
                            if(position>=coding_start):
                                db_seq=tmp_seq[coding_start-1:]
                                if(db_seq not in trans_db_annotation_dict.keys()):
                                    trans_db_annotation_dict[db_seq]=[key,str(k)+"_exclusive_from_"+str(position-coding_start+1)]
                                else:
                                    des_tmp=str(k)+"_exclusive_from_"+str(position-coding_start+1)
                                    if(des_tmp not in trans_db_annotation_dict[db_seq]):
                                        trans_db_annotation_dict[db_seq].append(des_tmp)
                                
                                tmp_anno=""
                                for d in hom_des:
                                    tmp_anno+=d
                                for d in het_des[c]:
                                    tmp_anno+=d
                                trans_db_annotation_dict[db_seq].append(tmp_anno) 
                                #tmp_info=exon_coordinate_dict[skip_id]
                                tmp_skip_id_index=generate_all_se_db.find_key_index('exon_id',skip_id,trans_records)
                                tmp_info=(trans_records.iloc[tmp_skip_id_index]['chr'],trans_records.iloc[tmp_skip_id_index]['start'],\
                                          trans_records.iloc[tmp_skip_id_index]['end'])
        
                                skip_seq=generate_all_se_db.fetch_exon_seq(tmp_info[0],tmp_info[1],tmp_info[2],T[0])
                                trans_prime_seq=generate_all_se_db.add_exon(tmp_seq,skip_seq,position)
                                db_prime_seq=trans_prime_seq[coding_start-1:]
                                if(db_prime_seq not in trans_db_annotation_dict.keys()):
                                    trans_db_annotation_dict[db_prime_seq]=[key,"modified_"+str(k)+"_inclusive_from_"\
                                                            +str(position-coding_start+1)+"_to_"+\
                                                            str(position-coding_start+1+len(skip_seq))]
                                else:
                                    des_tmp="modified_"+str(k)+"_inclusive_from_"+\
                                                            str(position-coding_start+1)+"_to_"+\
                                                            str(position-coding_start+1+len(skip_seq))
                                    if(des_tmp not in trans_db_annotation_dict[db_prime_seq]):
                                        trans_db_annotation_dict[db_prime_seq].append(des_tmp)
                                
                                tmp_anno=""
                                for d in hom_des:
                                    tmp_anno+=d
                                for d in het_des[c]:
                                    tmp_anno+=d
                                trans_db_annotation_dict[db_prime_seq].append(tmp_anno)
                                
                        elif(T[low_index+1]==skip_id):
    
                            skip_exon_start=tmp_exon_position_dict[low_index+1][0]
                            skip_exon_end=tmp_exon_position_dict[low_index+1][1]
                            
                            if(skip_exon_start>=coding_start):
                                db_seq=tmp_seq[coding_start-1:]
                                if(db_seq not in trans_db_annotation_dict.keys()):
                                    trans_db_annotation_dict[db_seq]=[key,str(k)+"_inclusive_from_"+str(skip_exon_start-coding_start+1)\
                                                            +"_to_"+str(skip_exon_end-coding_start+1)]
                                else:
                                    des_tmp=str(k)+"_inclusive_from_"+str(skip_exon_start-coding_start+1)\
                                                            +"_to_"+str(skip_exon_end-coding_start+1)
                                    if(des_tmp not in trans_db_annotation_dict[db_seq]):
                                        trans_db_annotation_dict[db_seq].append(des_tmp)
                                
                                tmp_anno=""
                                for d in hom_des:
                                    tmp_anno+=d
                                for d in het_des[c]:
                                    tmp_anno+=d
                                trans_db_annotation_dict[db_seq].append(tmp_anno)
                                
                                trans_prime_seq=generate_all_se_db.remove_exon(tmp_seq,skip_exon_start,skip_exon_end)
                                db_prime_seq=trans_prime_seq[coding_start-1:]
                                if(db_prime_seq not in trans_db_annotation_dict.keys()):
                                    trans_db_annotation_dict[db_prime_seq]=[key,"modified_"+str(k)+"_exclusive_from_"+str(skip_exon_start-coding_start+1)]
                                else:
                                    des_tmp="modified_"+str(k)+"_exclusive_from_"+\
                                                            str(skip_exon_start-coding_start+1)
                                    if(des_tmp not in trans_db_annotation_dict[db_prime_seq]):
                                        trans_db_annotation_dict[db_prime_seq].append(des_tmp)
                                
                                tmp_anno=""
                                for d in hom_des:
                                    p=get_variant_pos_from_des(d)
                                    if p not in range(skip_exon_start,skip_exon_end):
                                        tmp_anno+=d
                                for d in het_des[c]:
                                    p=get_variant_pos_from_des(d)
                                    if p not in range(skip_exon_start,skip_exon_end):
                                        tmp_anno+=d
                                trans_db_annotation_dict[db_prime_seq].append(tmp_anno)
                    
                if((len(hom_position_list)==0) and (len(het_list)==0)):
                    coding_start=trans_coding_dict[key][0]
                
                    tmp_exon_coordinate_dict={}
                    tmp_exon_position_dict={}
                    tmp_pos=0
                    for item in range(1,len(T)):
                        ei=generate_all_se_db.find_key_index('exon_id',T[item],trans_records)
                        if(ei==-1):
                            print("the exon id cannot be found")
                        else:
                            tmp_start=trans_records.iloc[ei]['start']
                            tmp_end=trans_records.iloc[ei]['end']
                            tmp_exon_coordinate_dict[T[item]]=(trans_records.iloc[ei]['chr'],tmp_start,tmp_end)
                            
                            tmp_exon_position_dict[item]=(tmp_pos,tmp_pos+tmp_end-tmp_start+1)
                            tmp_pos+=tmp_end-tmp_start+1
                    
                    if((high_index-low_index)==1):
                        position=tmp_exon_position_dict[low_index][1]
                        if(position>=coding_start):
                            db_seq=tmp_seq[coding_start-1:]
                            if(db_seq not in trans_db_annotation_dict.keys()):
                                trans_db_annotation_dict[db_seq]=[key,str(k)+"_exclusive_from_"+str(position-coding_start+1)]
                            else:
                                des_tmp=str(k)+"_exclusive_from_"+str(position-coding_start+1)
                                if(des_tmp not in trans_db_annotation_dict[db_seq]):
                                    trans_db_annotation_dict[db_seq].append(des_tmp)
    
                            tmp_skip_id_index=generate_all_se_db.find_key_index('exon_id',skip_id,trans_records)
                            tmp_info=(trans_records.iloc[tmp_skip_id_index]['chr'],trans_records.iloc[tmp_skip_id_index]['start'],\
                                      trans_records.iloc[tmp_skip_id_index]['end'])
    
                            skip_seq=generate_all_se_db.fetch_exon_seq(tmp_info[0],tmp_info[1],tmp_info[2],T[0])
                            trans_prime_seq=generate_all_se_db.add_exon(tmp_seq,skip_seq,position)
                            db_prime_seq=trans_prime_seq[coding_start-1:]
                            if(db_prime_seq not in trans_db_annotation_dict.keys()):
                                trans_db_annotation_dict[db_prime_seq]=[key,"modified_"+str(k)+"_inclusive_from_"\
                                                        +str(position-coding_start+1)+"_to_"+\
                                                        str(position-coding_start+1+len(skip_seq))]
                            else:
                                des_tmp="modified_"+str(k)+"_inclusive_from_"+\
                                                        str(position-coding_start+1)+"_to_"+\
                                                        str(position-coding_start+1+len(skip_seq))
                                if(des_tmp not in trans_db_annotation_dict[db_prime_seq]):
                                    trans_db_annotation_dict[db_prime_seq].append(des_tmp)
                    elif(T[low_index+1]==skip_id):
                        skip_exon_start=tmp_exon_position_dict[low_index+1][0]
                        skip_exon_end=tmp_exon_position_dict[low_index+1][1]
                        
                        if(skip_exon_start>=coding_start):
                            db_seq=tmp_seq[coding_start-1:]
                            if(db_seq not in trans_db_annotation_dict.keys()):
                                trans_db_annotation_dict[db_seq]=[key,str(k)+"_inclusive_from_"+str(skip_exon_start-coding_start+1)\
                                                        +"_to_"+str(skip_exon_end-coding_start+1)]
                            else:
                                des_tmp=str(k)+"_inclusive_from_"+str(skip_exon_start-coding_start+1)\
                                                        +"_to_"+str(skip_exon_end-coding_start+1)
                                if(des_tmp not in trans_db_annotation_dict[db_seq]):
                                    trans_db_annotation_dict[db_seq].append(des_tmp)
                            
                            trans_prime_seq=generate_all_se_db.remove_exon(tmp_seq,skip_exon_start,skip_exon_end)
                            db_prime_seq=trans_prime_seq[coding_start-1:]
                            if(db_prime_seq not in trans_db_annotation_dict.keys()):
                                trans_db_annotation_dict[db_prime_seq]=[key,"modified_"+str(k)+"_exclusive_from_"+str(skip_exon_start-coding_start+1)]
                            else:
                                des_tmp="modified_"+str(k)+"_exclusive_from_"+\
                                                        str(skip_exon_start-coding_start+1)
                                if(des_tmp not in trans_db_annotation_dict[db_prime_seq]):
                                    trans_db_annotation_dict[db_prime_seq].append(des_tmp)
                    
                    
                    
                    
                    
    trans_records=trans_records.sort_values(by=['trans_id'])
    my_seqs=[] 
   # my_transcripts=[]  
    modified_trans=[]
    mutation_seq_cnt=0                   
    for k in trans_db_annotation_dict.keys():
        annotation_list=trans_db_annotation_dict[k]
        trans_id=annotation_list[0]
        modified_trans.append(trans_id)
        #coding_start=trans_coding_dict[trans_id][0]
        #print(trans_id)
        pid=protein_id_dict[trans_id][0]
        gid=protein_id_dict[trans_id][1]
        
        tmp_id=db+"|"+pid+'|'+trans_id+'|'+gid
        
        des=''
        flag=0
        for an in range(1,len(annotation_list)):
            if(annotation_list[an].find('snv')!=-1):
                flag=1
            des+=annotation_list[an]+'|'
#            my_transcripts.append(SeqRecord(Seq(str(k),IUPAC.ambiguous_dna),id=tmp_id,description=des))
        if(flag==1):
            mutation_seq_cnt+=1
        new_seq=str(Seq(str(k),IUPAC.ambiguous_dna).transcribe().translate(to_stop=True))
            
        while(new_seq.find('None')!=-1):
            new_seq=new_seq.replace('None','')
        my_seqs.append(SeqRecord(Seq(str(new_seq),IUPAC.protein),id=tmp_id,description=des))
            
    print("The number of sequences with mutation generated is "+str(mutation_seq_cnt))         

    
    if(exclude!=True):
        original_cnt=0
        for t in trans_seq_dict.keys():
            if t not in modified_trans:
                
                coding_start=int(trans_coding_dict[t][0])
                coding_end=int(trans_coding_dict[t][1])
                trans_seq=str(trans_seq_dict[t])
                trans_seq=trans_seq[coding_start-1:]
                if(trans_seq not in trans_db_annotation_dict.keys()):
                
                    pid=protein_id_dict[t][0]
                    gid=protein_id_dict[t][1]
                    new_sequence=trans_seq
                    
                    new_seq=str(Seq(str(new_sequence),IUPAC.ambiguous_dna).transcribe().translate(to_stop=True))
        #            else:
        #                new_seq=str(Seq(str(new_sequence),IUPAC.ambiguous_dna).complement().transcribe().translate(to_stop=True))
                    
                    my_seqs.append(SeqRecord(Seq(str(new_seq),IUPAC.protein),\
                                                 id=db+'|'+pid+'|'+t+'|'+gid+'_0:'+str(coding_start)+'-'+str(coding_end),\
                                             description="no splicing"))
                    original_cnt+=1
        
        print("The number of original sequences is "+str(original_cnt))
    
    print("The number of sequences generated is "+str(len(my_seqs)))
    
    handle=open(output_name+".fasta","w")
    #handle=open(dataset_name+"_all_se_db"+".fasta","w")
    for sequence in my_seqs:
        SeqIO.write(sequence,handle,"fasta")
    handle.close()
#    handle_trans=open("../data/"+dataset_name+"/"+dataset_name+"_all_se_transcripts"+".fasta","w")
#    for sequence in my_transcripts:
#        SeqIO.write(sequence,handle_trans,"fasta")
#    handle_trans.close()
    #return my_seqs

#dbname="../data/gencode.v28.basic.annotation.gff3"   
#SE_name=argv[1] 
#dfname=argv[2]  
#rna_db=argv[3] 
#output_name=argv[4]
#
#df_SE=pd.read_csv(SE_name,sep='\t') 
#generate_se_sequence_combination(df_SE,dfname,dbname,rna_db,output_name)      
##
#finish=time.time()
#print("---- %s minutes ----" % ((finish-begin)/60))               