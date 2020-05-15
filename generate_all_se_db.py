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
import warnings
warnings.filterwarnings("ignore")

argv=sys.argv
begin=time.time()

#chromosome_dict={'chr1':'NC000001','chr2':'NC000002','chr3':'NC000003','chr4':'NC000004','chr5':'NC000005',\
#                  'chr6':'NC000006','chr7':'NC000007','chr8':'NC000008','chr9':'NC000009','chr10':'NC000010',\
#                  'chr11':'NC000011','chr12':'NC000012','chr13':'NC000013','chr14':'NC000014','chr15':'NC000015',\
#                  'chr16':'NC000016','chr17':'NC000017','chr18':'NC000018','chr19':'NC000019','chr20':'NC000020',\
#                  'chr21':'NC000021','chr22':'NC000021','chrX':'NC000023','chrY':'NC000024'}

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
        

def get_trans_records(dbname):
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
    
    for record in records:
        if('transcript_type' in record.attributes):
            if(record.attributes['transcript_type'][0]=='protein_coding'):
                if(record[2]=='exon'):
                    chromosome.append(record[0])
                    strand.append(record[6])
                    start.append(int(record[3]))
                    end.append(int(record[4]))
                    tmp=record.attributes['ID'][0]
                    flag=tmp.find(':')
                    tmp=tmp[flag+1:]
                    s=tmp.split(':',1)
                    trans_id.append(s[0])
                    exon_number.append(int(s[1]))
                    gene_id.append(record.attributes['gene_id'][0])
                    exon_id.append(record.attributes['exon_id'][0])
                    protein_id.append(record.attributes['protein_id'][0])
                    e+=1
    
    print(e)
    exon_records={'chr':chromosome,'strand':strand,'start':start,'end':end,'trans_id':trans_id,\
                      'exon_number':exon_number,'exon_id':exon_id,'protein_id':protein_id,'gene_id':gene_id}
    trans_records=pd.DataFrame(exon_records,columns=['chr','strand','start','end','trans_id','exon_number',\
                                                     'exon_id','protein_id','gene_id'])
    trans_records=trans_records.sort_values(by=['trans_id','exon_number'],axis=0)
    
    return trans_records

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
    print(chr_count)
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
    print(chr_index)
    
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
        if(i%10000==0):
            print(i)
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
    rna_seqs=SeqIO.parse(rna_db,'fasta')
    
    trans_seq_dict={}
    trans_coding_dict={}#while translating, start from coding_start-1

    for seq in rna_seqs:
        tmp=seq.id
        flag=tmp.find("|")
        mrna_id=tmp[:flag]
        tmp=tmp[flag+1:]
        flag=tmp.find("|")
        
        cds=tmp.find('CDS:')
        tmp=tmp[cds:]
        cds_end=tmp.find('|')
        tmp=tmp[tmp.find(':')+1:cds_end]
        
        split_flag=tmp.find('-')
        coding_start=int(tmp[:split_flag])
        tmp=tmp[split_flag+1:]
        coding_end=int(tmp)
        
        trans_seq_dict[mrna_id]=seq.seq
        trans_coding_dict[mrna_id]=(coding_start,coding_end)
        
    
    return trans_seq_dict,trans_coding_dict

def add_exon(transcript,skip_seq,position): #position is the end of the low index exon
    return transcript[:position]+skip_seq+transcript[position:]

def remove_exon(transcript,start,end):
    return transcript[:start]+transcript[end:]

def generate_se_sequence(SE_name,dbname,rna_db,dataset_name):
    df_SE=pd.read_csv(SE_name,sep='\t')
    
    records=get_trans_records(dbname)
    
    splicing_event_dict=get_splicing_event_dict(df_SE,records)
    print("splicing events dict ready")
    print(len(list(splicing_event_dict.items())))
    for x in list(splicing_event_dict.items())[0:3]:
        print(x)
        
    trans_exon_dict=get_transcript_exon_dict(records)
    print("transcript exon dict ready")
    print(len(list(trans_exon_dict.items())))
    for x in list(trans_exon_dict.items())[0:3]:
        print(x)
    #print(len(list(exon_coordinate_dict.items())))
    #print(len(list(trans_exon_position_dict.items())))
    
    trans_seq_dict,trans_coding_dict=get_trans_seq_dict(rna_db)
    print("transcript sequence dict ready")
    print(len(list(trans_seq_dict.items())))
    
    middle=time.time()
    print("---- %s minutes ----" % ((middle-begin)/60))
    
    trans_db_annotation_dict={}##output
    cnt=0
    records=records.sort_values(by=['exon_id'],axis=0)
    
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
                
                tmp_exon_coordinate_dict={}
                tmp_exon_position_dict={}
                tmp_pos=0
                for item in range(1,len(T)):
                    ei=find_key_index('exon_id',T[item],records)
                    if(ei==-1):
                        print("the exon id cannot be found")
                    else:
                        tmp_start=records.iloc[ei]['start']
                        tmp_end=records.iloc[ei]['end']
                        tmp_exon_coordinate_dict[T[item]]=(records.iloc[ei]['chr'],tmp_start,tmp_end)
                        
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
                    
                        #tmp_info=exon_coordinate_dict[skip_id]
                        tmp_skip_id_index=find_key_index('exon_id',skip_id,records)
                        tmp_info=(records.iloc[tmp_skip_id_index]['chr'],records.iloc[tmp_skip_id_index]['start'],\
                                  records.iloc[tmp_skip_id_index]['end'])
                        #print(tmp_info)
#                        print(T[0],chromosome_dict[tmp_info[0]])
#                        handle=Entrez.efetch(db="nucleotide", id=str(chromosome_dict[tmp_info[0]]), rettype="fasta", \
#                                             strand=int(T[0]), seq_start=tmp_info[1], seq_stop=tmp_info[2])
#                        skip_seq = str(SeqIO.read(handle, "fasta").seq)
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
                        
                        trans_prime_seq=remove_exon(tmp_seq,skip_exon_start,skip_exon_end)
                        db_prime_seq=trans_prime_seq[coding_start-1:]
                        if(db_prime_seq not in trans_db_annotation_dict.keys()):
                            trans_db_annotation_dict[db_prime_seq]=[key,"modified_"+str(k)+"_exclusive_from_"+str(skip_exon_start-coding_start+1)]
                        else:
                            des_tmp="modified_"+str(k)+"_exclusive_from_"+\
                                                    str(skip_exon_start-coding_start+1)
                            if(des_tmp not in trans_db_annotation_dict[db_prime_seq]):
                                trans_db_annotation_dict[db_prime_seq].append(des_tmp)
    
    records=records.sort_values(by=['trans_id'])
    my_seqs=[] 
    my_transcripts=[]                     
    for k in trans_db_annotation_dict.keys():
        annotation_list=trans_db_annotation_dict[k]
        trans_id=annotation_list[0]
        #coding_start=trans_coding_dict[trans_id][0]
        #print(trans_id)
        trans_index=find_key_index('trans_id',trans_id,records)
        if(trans_index!=-1):
            tmp_id="basic|"+records.iloc[trans_index]['protein_id']+'|'+trans_id+'|'+records.iloc[trans_index]['gene_id']
            #tmp_strand=trans_exon_dict[trans_id][0]
            des=''
            for an in range(1,len(annotation_list)):
                des+=annotation_list[an]+'|'
            my_transcripts.append(SeqRecord(Seq(str(k),IUPAC.ambiguous_dna),id=tmp_id,description=des))
            
            new_seq=str(Seq(str(k),IUPAC.ambiguous_dna).transcribe().translate(to_stop=True))
            
            while(new_seq.find('None')!=-1):
                new_seq=new_seq.replace('None','')
            my_seqs.append(SeqRecord(Seq(str(new_seq),IUPAC.protein),id=tmp_id,description=des))
            
               
    print("The number of sequences generated is "+str(len(my_seqs)))
    
    handle=open("../data/"+dataset_name+"/"+dataset_name+"_all_se_db"+".fasta","w")
    for sequence in my_seqs:
        SeqIO.write(sequence,handle,"fasta")
    handle.close()
    handle_trans=open("../data/"+dataset_name+"/"+dataset_name+"_all_se_transcripts"+".fasta","w")
    for sequence in my_transcripts:
        SeqIO.write(sequence,handle_trans,"fasta")
    handle_trans.close()
    #return my_seqs

#dbname="../data/gencode.v28.basic.annotation.gff3"
#rna_db="../data/gencode.v28.pc_transcripts.fa"     
SE_name=argv[1] 
dbname=argv[2]
rna_db=argv[3]
dataset_name=argv[4]        
generate_se_sequence(SE_name,dbname,rna_db,dataset_name)      

finish=time.time()
print("---- %s minutes ----" % ((finish-begin)/60))                    