#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: wenrchen
"""
import sys
import pandas as pd
import time
from gffutils.iterators import DataIterator
from Bio import SeqIO
#from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import pybedtools
import random
import generate_all_mutation_db

import warnings
warnings.filterwarnings("ignore")

argv=sys.argv
begin=time.time()

fasta = pybedtools.example_filename('GRCh38.primary_assembly.genome.fa')
def fetch_exon_seq(chr_,start,end,strand):
    as_str = ' '.join([chr_, str(start-1), str(end)])
    a = pybedtools.BedTool(as_str, from_string=True)
    a = a.sequence(fi=fasta)
    out=''
    handle=open(a.seqfn)
    if(strand==1):
        out=handle.read().split('\n')[1]
    else:
        s=handle.read().split('\n')[1]
        reverse_s=s[len(s)::-1]
        out=str(Seq(reverse_s,IUPAC.ambiguous_dna).complement())
    handle.close()
    return out
        

def get_trans_records(dbname,trans_seq_dict):
    records=DataIterator(dbname)
    
    e=0
    chromosome=[]
    strand=[]
    start=[]
    end=[]
    trans_id=[]
    exon_number=[]
    exon_id=[]
    gene_id=[]
    protein_id=[]
    strand_dict={}
    protein_id_dict={}
    
    for record in records:
        if('transcript_type' in record.attributes):
            if(record.attributes['transcript_id'][0] in trans_seq_dict.keys()):
                strand_dict[record.attributes['transcript_id'][0]]=record[6]
                protein_id_dict[record.attributes['transcript_id'][0]]=(record.attributes['protein_id'][0],\
                                   record.attributes['gene_id'][0])
                if(record[2]=='exon'):
                    chromosome.append(record[0])
                    strand.append(record[6])
                    start.append(int(record[3]))
                    end.append(int(record[4]))
                    tmp=record.attributes['ID'][0]
                    flag=tmp.find(':')
                    tmp=tmp[flag+1:]
                    s=tmp.split(':',1)
                    if(s[0].find('_')!=-1):
                        flag=s[0].find('_')
                        s[0]=s[0][:flag]
                    trans_id.append(s[0])
                    exon_number.append(int(s[1]))
                    gene_id.append(record.attributes['gene_id'][0])
                    exon_id.append(record.attributes['exon_id'][0])
                    protein_id.append(record.attributes['protein_id'][0])
                    
                    e+=1
    
    #print(e)
    exon_records={'chr':chromosome,'strand':strand,'start':start,'end':end,'trans_id':trans_id,\
                      'exon_number':exon_number,'exon_id':exon_id,'protein_id':protein_id,'gene_id':gene_id}
    trans_records=pd.DataFrame(exon_records,columns=['chr','strand','start','end','trans_id','exon_number',\
                                                     'exon_id','protein_id','gene_id'])
    trans_records=trans_records.sort_values(by=['trans_id','exon_number'],axis=0)
    
    return trans_records,strand_dict,protein_id_dict

def find_key_index(key,key_id,df):
    #df=df.sort_values(by=[key],axis=0)
    low=0
    high=df.shape[0]-1
    while(low<=high):
        mid=(low+high)//2
        if(df.iloc[mid][key]<key_id):
            low=mid+1
        elif(df.iloc[mid][key]>key_id):
            high=mid-1
        else:
            return mid
    return -1

def get_splicing_event_dict(df_SE,records): #df of SE merge result
    df_SE=df_SE.sort_values(by=['chr','upstreamES'],axis=0)
    
    records=records.sort_values(by=['chr','start'],axis=0)
    
    chromosome=[]
    strand=[]
    start=[]
    end=[]
    #ID=[]
    exon_id=[]
    gene_id=[]
    protein_id=[]
    
    
    for i in range(0,records.shape[0]):
        chromosome.append(records.iloc[i]['chr'])   
        strand.append(records.iloc[i]['strand'])
        start.append(records.iloc[i]['start'])
        end.append(records.iloc[i]['end'])
        #ID.append(records.iloc[i]['ID'])
        exon_id.append(records.iloc[i]['exon_id'])
        gene_id.append(records.iloc[i]['gene_id'])
        protein_id.append(records.iloc[i]['protein_id'])
    exon_sorted_records={'chr':chromosome,'strand':strand,'start':start,'end':end,'exon_id':exon_id,\
                      'gene_id':gene_id,'protein_id':protein_id}
    print('Annotation records are ready.')
    chr_name=set(chromosome)
    chr_count=dict((c,chromosome.count(c)) for c in chr_name)
    #print(chr_count)
    chr_count_name=[n for n in sorted(chr_count.keys())]
    chr_count_list=[chr_count[l] for l in sorted(chr_count.keys())]
    chr_index=dict((c,(0,0)) for c in chr_name)
    for ch in range(0,len(chr_count_name)):
        if(ch==0):
            chr_index[chr_count_name[ch]]=(0,chr_count_list[ch]-1)
        else:
            chr_index[chr_count_name[ch]]=(sum(chr_count_list[h] for h in range(0,ch)),\
                                     sum(chr_count_list[h] for h in range(0,ch+1))-1)
    del chr_index['chrM']
    #print(chr_index)
    
    def chr_binary_search(chro,strand,start,end):
        low=0
        high=0
        mid=0
        if chro in chr_index:
            low,high=chr_index[chro]
        else:
            print("The chromosome name is incorrect.")
            return -1
        while(low <= high):
               # Find the midpoint of the sequence.
            mid = (high + low) // 2
            if((start == exon_sorted_records['start'][mid]) and (end==exon_sorted_records['end'][mid])):
                if(exon_sorted_records['strand'][mid]==strand):
                    return mid
                else:
                    return -1
            elif(start<exon_sorted_records['start'][mid]):
                high = mid-1
            else:
                low = mid+1
        return -1
    
    result_dict={}
    for i in range(0,df_SE.shape[0]):
        tmp_id=df_SE.iloc[i]['ID']
        skip_target=chr_binary_search(df_SE.iloc[i]['chr'],df_SE.iloc[i]['strand'],\
                                          df_SE.iloc[i]['exonStart'],df_SE.iloc[i]['exonEnd'])
        up_target=chr_binary_search(df_SE.iloc[i]['chr'],df_SE.iloc[i]['strand'],\
                                      df_SE.iloc[i]['upstreamES'],df_SE.iloc[i]['upstreamEE'])
        down_target=chr_binary_search(df_SE.iloc[i]['chr'],df_SE.iloc[i]['strand'],\
                                      df_SE.iloc[i]['downstreamES'],df_SE.iloc[i]['downstreamEE'])
        if(skip_target!=-1 and up_target!=-1 and down_target!=-1):
            result_dict[tmp_id]=(exon_sorted_records['exon_id'][down_target],\
                       exon_sorted_records['exon_id'][skip_target],exon_sorted_records['exon_id'][up_target])
        
    return result_dict

def get_transcript_exon_dict(records):
    trans_exon_dict={}
    #trans_strand_dict={}
    #exon_coordinate_dict={}
    #trans_exon_position_dict={}
    
    i=0
    while(i<records.shape[0]):
#        if(i%10000==0):
#            print(i)
        k=records.iloc[i]['trans_id']
        
        if(records.iloc[i]['strand']=='+'):
            trans_exon_dict[k]=[1]
        else:
            trans_exon_dict[k]=[2]
        while((i<records.shape[0]) and (records.iloc[i]['trans_id']==k)):
            trans_exon_dict[k].append(records.iloc[i]['exon_id'])
            i+=1
        
#    for i in range(0,records.shape[0]):
#        k=records.iloc[i]['trans_id']
#        
#        if(records.iloc[i]['strand']=='+'):
#            #trans_strand_dict[k]=1
#            trans_exon_dict[k]=[1]
#        else:
#            #trans_strand_dict[k]=2
#            trans_exon_dict[k]=[2]
#        while(i<records.shape[0] and records.iloc[i]['trans_id']==k):
##            if(k not in trans_exon_dict.keys()):
##                trans_exon_dict[k]=[records.iloc[i]['exon_id']]
##                    
##            else:
#            trans_exon_dict[k].append(records.iloc[i]['exon_id'])
#            #exon_coordinate_dict[records.iloc[i]['exon_id']]=(records.iloc[i]['chr'],records.iloc[i]['start'],\
#                                    #records.iloc[i]['end'])
#            i+=1
            
#    for k in trans_exon_dict.keys():
#        tmp_pos=0
#        tmp_list=trans_exon_dict[k]
#        for v in range(1,len(tmp_list)):
#            tmp_exon=exon_coordinate_dict[tmp_list[v]]
#            if(tmp_pos==0):
#                
#                trans_exon_position_dict[k]=[(0,(tmp_exon[2]-tmp_exon[1]+1))]
#                #tmp_pos+=tmp_exon[2]-tmp_exon[1]+1
#            else:
#                trans_exon_position_dict[k].append((tmp_pos,(tmp_pos+tmp_exon[2]-tmp_exon[1]+1)))
#            tmp_pos+=tmp_exon[2]-tmp_exon[1]+1
            
    #return trans_exon_dict,exon_coordinate_dict,trans_exon_position_dict
    return trans_exon_dict


def get_trans_seq_dict(rna_db):
#    rna_seq=SeqIO.parse(rna_db,'fasta')
# 
#    for correct in rna_seq:
#        tmp=correct.id
#        flag=tmp.find("|")
#        tmp=tmp[flag+1:]
#        cds=tmp.find('CDS:')
#        if(cds==-1):
#            print("The format of file of parameter -r(--rna) is incorrect!")
#            return 0
#        tmp=tmp[cds:]
#        cds_end=tmp.find('|')
#        tmp=tmp[tmp.find(':')+1:cds_end]
#        
#        split_flag=tmp.find('-')
#        if(split_flag==-1):
#            print("The format of file of parameter -r(--rna) is incorrect!")
#            return 0

    rna_seqs=SeqIO.parse(rna_db,'fasta') 
       
    trans_seq_dict={}
    trans_coding_dict={}#while translating, start from coding_start-1
    
    for seq in rna_seqs:
        tmp=seq.id
        flag=tmp.find("|")
        mrna_id=tmp[:flag]

        tmp=tmp[flag+1:]        
        cds=tmp.find('CDS:')
        tmp=tmp[cds:]
        cds_end=tmp.find('|')
        tmp=tmp[tmp.find(':')+1:cds_end]
        
        split_flag=tmp.find('-')
        coding_start=int(tmp[:split_flag])
        coding_end=int(tmp[split_flag+1:])
        
        trans_seq_dict[mrna_id]=seq.seq
        trans_coding_dict[mrna_id]=(coding_start,coding_end)
        
    
    return trans_seq_dict,trans_coding_dict

def add_exon(transcript,skip_seq,position): #position is the end of the low index exon
    return transcript[:position]+skip_seq+transcript[position:]

def remove_exon(transcript,start,end):
    return transcript[:start]+transcript[end:]

def get_variant_pos_from_des(des):
    start_flag=des.find(':')
    end_flag=des.find('-')
    return int(des[start_flag+1:end_flag])


def generate_random_splicing(trans_seq,splicing_type,length,position):
    alphabet='AGCT'
    #alphabet="ARNDCEQGHILKMFPSTWYV"
    
    if(len(trans_seq)!=0):
        random_pos=position
        while(random_pos==position):
            random_pos=random.randint(0,len(trans_seq)-1)
        
        if(splicing_type=="remove"):
            new_seq=trans_seq[:random_pos]+trans_seq[random_pos+length:]
            random_des="random_exclusive_from_"+str(random_pos)
        elif(splicing_type=='add'):
            random_seq=""
            for i in range(0,int(length/3)):
                random_codon=""
                while(len(random_codon)!=3):
                    random_codon+=alphabet[random.randint(0,3)]
                    if((random_codon=="TAA") or (random_codon=="TAG") or (random_codon=="TGA")):
                        random_codon=""
                random_seq=random_seq+random_codon
            
            
            residue=length%3
            for r in range(0,residue):
                random_seq=random_seq+alphabet[random.randint(0,3)]
            new_seq=trans_seq[:random_pos]+random_seq+trans_seq[random_pos:]
            random_des="random_inclusive_from_"+str(random_pos)+"_to_"+str(random_pos+length)
                
        else:
            print("The splicing type for random is incorrect.")
                    
    else:
        new_seq=""
        random_des=""
    return new_seq,random_des

def get_version(rna_db):
    version=rna_db
    ver=version.find("trans")
    version=version[ver:]
    ver=version.find(".")
    version=version[ver+1:]
    ver=version.find(".")
    version=version[:ver]
    
    return version


def generate_se_sequence(df_SE, dfname, dbname,rna_db,exclude,output_name):
    version=get_version(rna_db)
    if(version=='swissprot'):
        db='sp'
    else:
        db=version
    
    
    trans_seq_dict,trans_coding_dict=get_trans_seq_dict(rna_db)
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
    
    trans_records,strand_dict,protein_id_dict=get_trans_records(dbname,trans_seq_dict)
    
    splicing_event_dict=get_splicing_event_dict(df_SE,trans_records)
    print("splicing events dict ready")

        
    trans_exon_dict=get_transcript_exon_dict(trans_records)
    print("transcript exon dict ready")
    trans_records=trans_records.sort_values(by=['exon_id'],axis=0)
    
    middle=time.time()
    print("---- %s minutes ----" % ((middle-begin)/60))
    
    
    trans_db_annotation_dict={}##output
    cnt=0
    random_cnt=0
    
    
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
                tmp_seq=trans_seq_dict[key]
                coding_start=trans_coding_dict[key][0]
                
                
                hom_des=[]
                if((key in trans_index_dict.keys()) and (key in strand_dict.keys())):
                    shift=0
                    transcript=""
                    hom_position_list=[]
                    
                    for i in trans_index_dict[key]:
                   
                        if((change_df.iloc[i]['snp_type']=='hom') and (change_df.iloc[i]['mutation_type'].find('snv')!=-1)):
                            hom_position_list.append((int(change_df.iloc[i]['c_start']),change_df.iloc[i]['mutation_type']))
    
                            if(strand_dict[key]=='-'):
                                tmp=Seq(str(tmp_seq),IUPAC.ambiguous_dna).complement()
                                tmp_seq=str(tmp)
                                transcript=generate_all_mutation_db.change_seq(tmp_seq,int(change_df.iloc[i]['c_start'])+shift+coding_start-1,\
                                                      int(change_df.iloc[i]['c_end'])+shift+coding_start-1,str(Seq(str(change_df.iloc[i]['c_content'])).complement()),\
                                                      change_df.iloc[i]['mutation_type'])
                            else:
                                transcript=generate_all_mutation_db.change_seq(tmp_seq,int(change_df.iloc[i]['c_start'])+shift+coding_start-1,\
                                                      int(change_df.iloc[i]['c_end'])+shift+coding_start-1,change_df.iloc[i]['c_content'],\
                                                      change_df.iloc[i]['mutation_type'])
                            hom_des.append(change_df.iloc[i]['mutation_type']+":"+str(change_df.iloc[i]['c_start'])+'-'+\
                            str(change_df.iloc[i]['c_end'])+str(change_df.iloc[i]['c_content'])+'_')

                        
                        if(len(hom_position_list)!=0):
                            if(strand_dict[key]=='-'):
                                tmp=Seq(str(transcript),IUPAC.ambiguous_dna).complement()
                                tmp_seq=str(tmp)
                            else:
                                tmp_seq=transcript
                        else:
                            tmp_seq=trans_seq_dict[key]
                
                tmp_exon_coordinate_dict={}
                tmp_exon_position_dict={}
                tmp_pos=0
                for item in range(1,len(T)):
                    ei=find_key_index('exon_id',T[item],trans_records)
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
                    #position=trans_exon_position_dict[key][low_index-1][1]
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
                        trans_db_annotation_dict[db_seq].append(tmp_anno)
                        #tmp_info=exon_coordinate_dict[skip_id]
                        tmp_skip_id_index=find_key_index('exon_id',skip_id,trans_records)
                        tmp_info=(trans_records.iloc[tmp_skip_id_index]['chr'],trans_records.iloc[tmp_skip_id_index]['start'],\
                                  trans_records.iloc[tmp_skip_id_index]['end'])

                        skip_seq=fetch_exon_seq(tmp_info[0],tmp_info[1],tmp_info[2],T[0])
                        trans_prime_seq=add_exon(tmp_seq,skip_seq,position)
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
                        trans_db_annotation_dict[db_prime_seq].append(tmp_anno)
                        
                        
                        db_random_seq,random_des=generate_random_splicing(db_seq,"add",len(skip_seq),position)
                        trans_db_annotation_dict[db_random_seq]=[key,random_des]
                        tmp_anno=""
                        for d in hom_des:
                            tmp_anno+=d
                        trans_db_annotation_dict[db_random_seq].append(tmp_anno)
                        random_cnt+=1
                        
                        db_random_prime_seq,random_prime_des=generate_random_splicing(db_prime_seq,"remove",len(skip_seq),position)
                        trans_db_annotation_dict[db_random_prime_seq]=[key,random_prime_des]
                        tmp_anno=""
                        for d in hom_des:
                            tmp_anno+=d
                        trans_db_annotation_dict[db_random_prime_seq].append(tmp_anno)
                        random_cnt+=1
                        
                elif(T[low_index+1]==skip_id):

                    #skip_exon_start=trans_exon_position_dict[key][low_index+1][0]
                    #skip_exon_end=trans_exon_position_dict[key][low_index+1][1]
                    skip_exon_start=tmp_exon_position_dict[low_index+1][0]
                    skip_exon_end=tmp_exon_position_dict[low_index+1][1]
                    
                    #coding_start=trans_coding_dict[key][0]
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
                        trans_db_annotation_dict[db_seq].append(tmp_anno)
                        
                        trans_prime_seq=remove_exon(tmp_seq,skip_exon_start,skip_exon_end)
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
                        trans_db_annotation_dict[db_prime_seq].append(tmp_anno)
                        
                        db_random_seq,random_des=generate_random_splicing(db_seq,"remove",skip_exon_end-skip_exon_start+1,skip_exon_start-coding_start+1)
                        trans_db_annotation_dict[db_random_seq]=[key,random_des]
                        tmp_anno=""
                        for d in hom_des:
                            p=get_variant_pos_from_des(d)
                            if p not in range(skip_exon_start,skip_exon_end):
                                tmp_anno+=d
                        trans_db_annotation_dict[db_random_seq].append(tmp_anno)
                        random_cnt+=1
                        
                        db_random_prime_seq,random_prime_des=generate_random_splicing(db_prime_seq,"add",skip_exon_end-skip_exon_start+1,skip_exon_start-coding_start+1)
                        trans_db_annotation_dict[db_random_prime_seq]=[key,random_prime_des]
                        tmp_anno=""
                        for d in hom_des:
                            p=get_variant_pos_from_des(d)
                            if p not in range(skip_exon_start,skip_exon_end):
                                tmp_anno+=d
                        trans_db_annotation_dict[db_random_prime_seq].append(tmp_anno)
                        random_cnt+=1
    
    trans_records=trans_records.sort_values(by=['trans_id'])
    my_seqs=[] 
    #my_transcripts=[] 
    modified_trans=[]                    
    for k in trans_db_annotation_dict.keys():
        annotation_list=trans_db_annotation_dict[k]
        trans_id=annotation_list[0]
        annotation_list.remove(trans_id)
        annotation_list=list(set(annotation_list))
        annotation_list.sort()
        modified_trans.append(trans_id)
        #coding_start=trans_coding_dict[trans_id][0]
        #print(trans_id)
        #trans_index=find_key_index('trans_id',trans_id,trans_records)
        pid=protein_id_dict[trans_id][0]
        gid=protein_id_dict[trans_id][1]
        
        tmp_id=db+"|"+pid+'|'+trans_id+'|'+gid
        #tmp_strand=trans_exon_dict[trans_id][0]
        des=''
        for an in range(1,len(annotation_list)):
            des+=annotation_list[an]+'|'
        #my_transcripts.append(SeqRecord(Seq(str(k),IUPAC.ambiguous_dna),id=tmp_id,description=des))
        
        new_seq=str(Seq(str(k),IUPAC.ambiguous_dna).transcribe().translate(to_stop=True))
        
        while(new_seq.find('None')!=-1):
            new_seq=new_seq.replace('None','')
        my_seqs.append(SeqRecord(Seq(str(new_seq),IUPAC.protein),id=tmp_id+'_'+des,description=des))
    
               
    print("The number of sequences generated is "+str(len(my_seqs)))
    print("The number of random sequences generated is "+str(random_cnt))
    
    if(exclude!=True):
        original_cnt=0
        for i in trans_seq_dict.keys():
            if i not in modified_trans:
                coding_start=int(trans_coding_dict[i][0])
                coding_end=int(trans_coding_dict[i][1])
                trans_seq=str(trans_seq_dict[i])
                trans_seq=trans_seq[coding_start-1:]
                if(trans_seq not in trans_db_annotation_dict.keys()):
                
                    pid=protein_id_dict[i][0]
                    gid=protein_id_dict[i][1]
                    new_sequence=trans_seq
                    
                    new_seq=str(Seq(str(new_sequence),IUPAC.ambiguous_dna).transcribe().translate(to_stop=True))
                    
                    my_seqs.append(SeqRecord(Seq(str(new_seq),IUPAC.protein),\
                                                 id=db+'|'+pid+'|'+i+'|'+gid+'_0:'+str(coding_start)+'-'+str(coding_end)+"_no_splicing",\
                                             description="no splicing"))
                    original_cnt+=1
        
        print("The number of original sequences is "+str(original_cnt))
    
    handle=open(output_name+".fasta","w")
    for sequence in my_seqs:
        SeqIO.write(sequence,handle,"fasta")
    handle.close()

