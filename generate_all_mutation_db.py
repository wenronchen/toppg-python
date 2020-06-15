#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: wenrchen
"""

import pandas as pd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from gffutils.iterators import DataIterator

argv=sys.argv

def change_seq(seq,start,end,content,t):
    if(t=='snv'):
        return seq[:start-1]+content+seq[end:]
    elif(t.find("del")!=-1):
        return seq[:start-1]+seq[end:]
    elif(t.find('ins')!=-1):
        return seq[:start-1]+content+seq[end:]

def get_protein_coding_list_from_db(protein_db):
    protein_coding_list=[]
    seqs = SeqIO.parse(protein_db,"fasta")
    for seq in seqs:
        tmp=seq.id
        flag=tmp.find("ENST")
        tmp=tmp[flag:]
        flag=tmp.find('|')
        protein_coding_list.append(tmp[:flag])
    return protein_coding_list

def extract_transcript_change(df):
    mrna=[]
    c_start=[]
    c_end=[]
    c_content=[]
    mutation_type=[]
    snp_type=[]
    
    for i in range(0,df.shape[0]):
        if(df.iloc[i][1]=='nonframeshift insertion'):
            tmp=df.iloc[i][2]
            tmp=tmp.split(',')
            for t in tmp:
                if(t.find(":")!=-1):
                    mrna_flag=t.find(":")
                    t=t[mrna_flag+1:]
                    exon_flag=t.find(":")
                    mrna.append(t[:exon_flag])
                    t=t[exon_flag+1:]
                    
                    c_flag=t.find(":")
                    t=t[c_flag+1:]
                    
                    ins_flag=t.find('ins')
                    mutation_type.append('nonfs_ins')
                    split_flag=t.find('_')
                    c_start.append(int(t[2:split_flag]))
                    c_end.append(int(t[split_flag+1:ins_flag]))
                    t=t[ins_flag+3:]
                    
                    p_flag=t.find(':')
                    c_content.append(t[:p_flag])
                    
                    snp_type.append(df.iloc[i][8])
                    
                    
        if(df.iloc[i][1]=='nonframeshift deletion'):
            tmp=df.iloc[i][2]
            tmp=tmp.split(',')
            for t in tmp:
                if(t.find(":")!=-1):
                    mrna_flag=t.find(":")
                    t=t[mrna_flag+1:]
                    exon_flag=t.find(":")
                    mrna.append(t[:exon_flag])
                    t=t[exon_flag+1:]
                    
                    c_flag=t.find(":")
                    t=t[c_flag+1:]
                    
                    type_flag=t.find("del")
                    mutation_type.append("nonfs_del")
                    split_flag=t.find('_')
                    c_start.append(t[2:split_flag])
                    c_end.append(t[split_flag+1:type_flag])
                    c_content.append("del")
                   
                    snp_type.append(df.iloc[i][8])
                    
        if(df.iloc[i][1]=='frameshift deletion'):
            tmp=df.iloc[i][2]
            tmp=tmp.split(',')
            for t in tmp:
                if(t.find(":")!=-1):
                    if(t.find("wholegene")==-1):
                        mrna_flag=t.find(":")
                        t=t[mrna_flag+1:]
                        exon_flag=t.find(":")
                        mrna.append(t[:exon_flag])
                        t=t[exon_flag+1:]
                    
                    
                        c_flag=t.find(":")
                        t=t[c_flag+1:]
                        type_flag=t.find("del")
                        mutation_type.append("fs_del")
                        
                        c_pos=t[2:type_flag]
                        if(c_pos.find('_')==-1):
                            c_start.append(int(c_pos))
                            c_end.append(int(c_pos))
                        else:
                            split_flag=c_pos.find('_')
                            c_start.append(int(c_pos[:split_flag]))
                            c_end.append(int(c_pos[split_flag+1:]))
                        c_content.append("del")
                        snp_type.append(df.iloc[i][8])
                        
                    
        if(df.iloc[i][1]=='frameshift insertion'):
            tmp_list=df.iloc[i][2].split(',')
            for tmp in tmp_list:
                
                if(tmp.find(":")!=-1):
                    mrna_flag=tmp.find(":")
                    tmp=tmp[mrna_flag+1:]
                    exon_flag=tmp.find(":")
                    mrna.append(tmp[:exon_flag])
                    tmp=tmp[exon_flag+1:]
                    
                    c_flag=tmp.find(":")
                    tmp=tmp[c_flag+1:]
                    
                    if(tmp.find('dup')!=-1):
                        dup_flag=tmp.find('dup')
                        mutation_type.append("fs_ins")
                        c_start.append(int(tmp[2:dup_flag]))
                        c_end.append(int(tmp[2:dup_flag])+1)
                        tmp=tmp[dup_flag+3:]
                
                    else:
                        ins_flag=tmp.find('ins')
                        mutation_type.append('fs_ins')
                        split_flag=tmp.find('_')
                        c_start.append(int(tmp[2:split_flag]))
                        c_end.append(int(tmp[split_flag+1:ins_flag]))
                        tmp=tmp[ins_flag+3:]
                     
                    flag=tmp.find(":")
                    c_content.append(tmp[:flag])
                    snp_type.append(df.iloc[i][8])
        
        if(df.iloc[i][1]=='nonsynonymous SNV'):
            tmp_list=df.iloc[i][2].split(',')
            for tmp in tmp_list:
                
                if(tmp.find(":")!=-1):
                    mrna_flag=tmp.find(":")
                    tmp=tmp[mrna_flag+1:]
                    exon_flag=tmp.find(":")
                    mrna.append(tmp[:exon_flag])
                    tmp=tmp[exon_flag+1:]
                    c_flag=tmp.find(":")
                    tmp=tmp[c_flag+1:]
                    p_flag=tmp.find(":")
                
                    c_flag_ref=tmp.find(".")
                    #print(tmp[c_flag_ref+2:p_flag-1])
                    c_start.append(int(tmp[c_flag_ref+2:p_flag-1]))
                    c_end.append(int(tmp[c_flag_ref+2:p_flag-1]))
                    c_content.append(tmp[p_flag-1:p_flag])
                    mutation_type.append('snv')
                    snp_type.append(df.iloc[i][8])
                
                    
    merge={'mrna':mrna,'mutation_type':mutation_type,'c_start':c_start,\
           'c_end':c_end,'c_content':c_content,'snp_type':snp_type}
    extract_result=pd.DataFrame(merge,columns=['mrna','mutation_type','c_start','c_end',\
                                               'c_content','snp_type'])
    extract_sorted=extract_result.sort_values(by=['mrna','c_start'])
    #extract_sorted.to_csv('../data/DLD/DLD_aa_change.tsv',sep='\t',index=False) 
    return extract_sorted

def get_version(protein_db):
    version=protein_db
    ver=version.find("v")
    version=version[ver:]
    ver=version.find(".")
    version=version[ver+1:]
    ver=version.find(".")
    version=version[:ver]
    
    return version
 

              
def get_new_sequence(dfname,dbname,rna_db,protein_db,dataset_name):
    version=get_version(protein_db)
    if(version=='swissprot'):
        db='sp'
    else:
        db=version
        
    protein_coding_list=get_protein_coding_list_from_db(protein_db)   
    records=DataIterator(dbname)
    strand_dict={}
    protein_id_dict={}
    for record in records:
        if(record[2]=='transcript'):
            if('transcript_type' in record.attributes):
                if(record.attributes['transcript_type'][0]=='protein_coding'):
                    if(record.attributes['transcript_id'][0] in protein_coding_list):
                        strand_dict[record.attributes['transcript_id'][0]]=record[6]
                        protein_id_dict[record.attributes['transcript_id'][0]]=(record.attributes['protein_id'][0],\
                                   record.attributes['gene_id'][0])
                    
    
    print("protein_id_dict ready")
              
    sequence_dict={}
    rna_seqs=SeqIO.parse(rna_db,'fasta')
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
        sequence_dict[mrna_id]=(coding_start,seq.seq,coding_end)
    
    print("sequence_dict ready")
    
    df=pd.read_csv(dfname,sep='\t',header=None)
    change_df=extract_transcript_change(df) 
    trans_index_dict={}
    for i in range(0,change_df.shape[0]):
        if(change_df.iloc[i]['mrna'] not in trans_index_dict.keys()):
            trans_index_dict[change_df.iloc[i]['mrna']]=[i]
        else:
            trans_index_dict[change_df.iloc[i]['mrna']].append(i)
    
    my_seqs=[]
    k_cnt=0
    hom_cnt=0
    het_cnt=0

    for k in trans_index_dict.keys():
        
        if(k in protein_id_dict.keys()):
            k_cnt+=1
            if(k_cnt%1000==0):
                print(k_cnt)
            pid=protein_id_dict[k][0]
            gid=protein_id_dict[k][1]
            
            if(strand_dict[k]=='+'):
            
                transcript=str(sequence_dict[k][1])
            else:
                tmp=Seq(str(sequence_dict[k][1]),IUPAC.ambiguous_dna).complement()
                transcript=str(tmp)
            coding_start=int(sequence_dict[k][0])-1 
            coding_end=int(sequence_dict[k][2])-1 
            
            transcript=transcript[coding_start:coding_end+1]
                
            shift=0
            des=""
            het_list=[]
            hom_position_list=[]
            for i in trans_index_dict[k]:
               
                if(change_df.iloc[i]['snp_type']=='hom'):
                    hom_cnt+=1
                    hom_position_list.append((int(change_df.iloc[i]['c_start']),change_df.iloc[i]['mutation_type']))
#                    if(change_df.iloc[i]['mutation_type']=='snv'):
#                        transcript=change_seq(transcript,int(change_df.iloc[i]['c_start'])+shift,\
#                                                  int(change_df.iloc[i]['c_end'])+shift,change_df.iloc[i]['c_content'],'snv')
#                        des+="snv:"+str(change_df.iloc[i]['c_start'])+change_df.iloc[i]['c_content']+'_'
#                    else:
                    if(strand_dict[k]=='-'):
                        transcript=change_seq(transcript,int(change_df.iloc[i]['c_start'])+shift,\
                                              int(change_df.iloc[i]['c_end'])+shift,str(Seq(str(change_df.iloc[i]['c_content'])).complement()),\
                                              change_df.iloc[i]['mutation_type'])
                    else:
                        transcript=change_seq(transcript,int(change_df.iloc[i]['c_start'])+shift,\
                                              int(change_df.iloc[i]['c_end'])+shift,change_df.iloc[i]['c_content'],\
                                              change_df.iloc[i]['mutation_type'])
                    des+=change_df.iloc[i]['mutation_type']+":"+str(change_df.iloc[i]['c_start'])+'-'+\
                    str(change_df.iloc[i]['c_end'])+str(change_df.iloc[i]['c_content'])+'_'
                    if(change_df.iloc[i]['mutation_type'].find('del')!=-1):
                        shift-=(int(change_df.iloc[i]['c_end'])-int(change_df.iloc[i]['c_start'])+1)
                    elif(change_df.iloc[i]['mutation_type'].find('ins')!=-1):
                        shift+=(int(change_df.iloc[i]['c_end'])-int(change_df.iloc[i]['c_start'])+1)
                else:
                    het_list.append(i)
                    het_cnt+=1
                    
            coding_start=int(sequence_dict[k][0])+shift-1 
            coding_end=int(sequence_dict[k][2])+shift-1
            count=int(len(transcript)/900)
            if(len(transcript)<=900):
                count=1
            cnt=0
            
            for l in range(0,count):
                l=l*900
                start=l
                if(start+1799<len(transcript)):
                    stop=start+1799
                else:
                    stop=len(transcript)-1
#                while(coding_start>stop):
#                    l=l+1
#                    start=l*1800
#                    if(l+1799<len(transcript)):
#                        stop=l+1799
#                    else:
#                        stop=len(transcript)-1
#                #if(start>coding_end+shift)  
                
                
                
                
                if(len(hom_position_list)!=0):
                    new_sequence=""
                    new_des=""
                    for p in hom_position_list:
                        if((p[0]-1>=start)&(p[0]-1<=stop)):
                            flag=des.find(str(p[0]))
                            tmp=des[flag:]
                            flag=tmp.find('_')
                            new_des+=str(p[1])+":"+tmp[:flag+1]
                            new_sequence=transcript
                            
                    
                    if(len(new_sequence)!=0):
                        
                        
                        new_sequence=new_sequence[start:stop+1]

                        if(strand_dict[k]=='+'):
                            new_seq=str(Seq(str(new_sequence),IUPAC.ambiguous_dna).transcribe().translate(to_stop=True))
                        else:
                            new_seq=str(Seq(str(new_sequence),IUPAC.ambiguous_dna).complement().transcribe().translate(to_stop=True))
                        cnt+=1
                        while(new_seq.find('None')!=-1):
                            new_seq=new_seq.replace('None','')
                        my_seqs.append(SeqRecord(Seq(str(new_seq),IUPAC.protein),\
                                                 id=db+'|'+pid+'|'+k+'|'+gid+'_'+str(cnt)+':'+str(start+1)+'-'+str(stop+1),\
                                             description=new_des))
                het_number=len(het_list)   
                for n in range(0,het_number):
                    new_sequence=""
                    new_des=""
                    
                    if((int(change_df.iloc[het_list[n]]['c_start'])>=start) & (int(change_df.iloc[het_list[n]]['c_start'])<=stop)):
                        if(strand_dict[k]=='-'):
                            new_sequence=change_seq(transcript,int(change_df.iloc[het_list[n]]['c_start'])+shift,\
                                                  int(change_df.iloc[het_list[n]]['c_end'])+shift,\
                                                  str(Seq(str(change_df.iloc[het_list[n]]['c_content'])).complement()),\
                                                  change_df.iloc[het_list[n]]['mutation_type'])
                        else:
                            new_sequence=change_seq(transcript,int(change_df.iloc[het_list[n]]['c_start'])+shift,\
                                                  int(change_df.iloc[het_list[n]]['c_end'])+shift,change_df.iloc[het_list[n]]['c_content'],\
                                                  change_df.iloc[het_list[n]]['mutation_type'])
                        new_des=des+change_df.iloc[het_list[n]]['mutation_type']+":"+\
                        str(change_df.iloc[het_list[n]]['c_start'])+'-'+str(change_df.iloc[het_list[n]]['c_end'])+\
                        str(change_df.iloc[het_list[n]]['c_content'])
                    if(len(new_sequence)!=0):
                        new_sequence=new_sequence[start:stop+1]
                       
                        if(strand_dict[k]=='+'):
                            new_seq=str(Seq(str(new_sequence),IUPAC.ambiguous_dna).transcribe().translate(to_stop=True))
                        else:
                            new_seq=str(Seq(str(new_sequence),IUPAC.ambiguous_dna).complement().transcribe().translate(to_stop=True))
                        cnt+=1
                        while(new_seq.find('None')!=-1):
                            new_seq=new_seq.replace('None','')
                        my_seqs.append(SeqRecord(Seq(str(new_seq),IUPAC.protein),\
                                                 id=db+'|'+pid+'|'+k+'|'+gid+'_'+str(cnt)+':'+str(start+1)+'-'+str(stop+1),\
                                             description=new_des))
                
                    
    print("The number of proteins related is "+str(k_cnt))            
    print("The number of sequences generated is "+str(len(my_seqs)))
    print("The number of homozygous is "+str(hom_cnt)) 
    print("The number of heterozygous is "+str(het_cnt))
    #return my_seqs
    handle=open("../data/"+dataset_name+"/"+version+"/"+dataset_name+"_all_mutation_"+version+".fasta","w")
    hom_seq=0
    het_seq=0
    hom_het_seq=0
    for sequence in my_seqs:
        if(str(sequence.description)[-1]=='_'):
            hom_seq+=1
        else:
            if(str(sequence.description).find('_')==-1):
                het_seq+=1
            else:
                hom_het_seq+=1
        
        SeqIO.write(sequence,handle,"fasta")
    print("The number of homozygous sequences is "+str(hom_seq)) 
    print("The number of heterozygous sequences is "+str(het_seq)) 
    print("The number of mixed sequences is "+str(hom_het_seq))       
            
#df=pd.read_csv(argv[1],sep='\t',header=None)  
#dataset_name=argv[2]
#db_annotation_name=argv[3]
#rna_db=argv[4]
#protein_db=argv[5]
#get_new_sequence(df,db_annotation_name,rna_db,protein_db,dataset_name)
     