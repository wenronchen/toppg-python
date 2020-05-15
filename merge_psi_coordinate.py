#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 17:13:39 2019

@author: wenrchen
"""

import pandas as pd

def merge_file(se_type,df_SE_name,gtf_SE_name):
    df_SE=pd.read_csv(df_SE_name,sep='\t') 
    gtf_SE=pd.read_csv(gtf_SE_name,sep='\t')
    
    ID=[]
    IJC=[]
    SJC=[]
    IncLen=[]
    SkipLen=[]
    IncLevel=[]
    gene_id=[]
    chromosome=[]
    strand=[]
    upES=[]#1_base
    upEE=[]
    downES=[]#1_base
    downEE=[]
    if(se_type=='MXE'):
        exon1_start=[]#to match the coordinate of gff file using 1_base
        exon1_end=[]
        exon2_start=[]
        exon2_end=[]
        
        j=0
        for i in range(0,df_SE.shape[0]):
            while(df_SE.iloc[i]['ID']>gtf_SE.iloc[j]['ID']):
                j+=1
            if(df_SE.iloc[i]['ID']==gtf_SE.iloc[j]['ID']):
            #if((df_SE.iloc[i]['ID']==gtf_SE.iloc[j]['ID']) & (df_SE.iloc[i]['IJC_SAMPLE_1']>10) & \
               #(df_SE.iloc[i]['SJC_SAMPLE_1']>10) & (df_SE.iloc[i]['IncLevel1']<0.9) & (df_SE.iloc[i]['IncLevel1']>0.1)):
                ID.append(df_SE.iloc[i]['ID'])
                IJC.append(df_SE.iloc[i]['IJC_SAMPLE_1'])
                SJC.append(df_SE.iloc[i]['SJC_SAMPLE_1'])
                IncLen.append(df_SE.iloc[i]['IncFormLen'])
                SkipLen.append(df_SE.iloc[i]['SkipFormLen'])
                IncLevel.append(df_SE.iloc[i]['IncLevel1'])
                tmp=gtf_SE.iloc[j]['GeneID']
                gene_id.append(tmp[0:-1])
                chromosome.append(gtf_SE.iloc[j]['chr'])
                strand.append(gtf_SE.iloc[j]['strand'])
                exon1_start.append(int(gtf_SE.iloc[j][5])+1)
                exon1_end.append(int(gtf_SE.iloc[j][6]))
                exon2_start.append(int(gtf_SE.iloc[j][7])+1)
                exon2_end.append(int(gtf_SE.iloc[j][8]))
                upES.append(int(gtf_SE.iloc[j][9])+1)
                upEE.append(int(gtf_SE.iloc[j][10]))
                downES.append(int(gtf_SE.iloc[j][11])+1)
                downEE.append(int(gtf_SE.iloc[j][12]))
            while(df_SE.iloc[i]['ID']<gtf_SE.iloc[j]['ID']):
                i-=1
        
        merge={'ID':ID,'IJC':IJC,'SJC':SJC,'IncLen':IncLen,'SkipLen':SkipLen,'IncLevel':IncLevel,\
               'GeneID':gene_id,'chr':chromosome,'strand':strand,'exon1Start':exon1_start,'exon1End':exon1_end,'exon2Start':exon2_start,'exon2End':exon2_end,\
               'upstreamES':upES,'upstreamEE':upEE,'downstreamES':downES,'downstreamEE':downEE}
        
        output=pd.DataFrame(merge,columns=['ID','IJC','SJC','IncLen','SkipLen','IncLevel',\
               'GeneID','chr','strand','exon1Start','exon1End','exon2Start','exon2End',\
               'upstreamES','upstreamEE','downstreamES','downstreamEE'])
        output.to_csv('../data/DLD/SE_search/'+se_type+'_merge_no_filter.tsv',sep='\t',index=None)
        return output
    elif(se_type=='SE' or 'RI' or 'A3SS' or 'A5SS'):

        exon_start=[]#to match the coordinate of gff file using 1_base
        exon_end=[]
        
        
        j=0
        for i in range(0,df_SE.shape[0]):
            while(df_SE.iloc[i]['ID']>gtf_SE.iloc[j]['ID']):
                j+=1
            if(df_SE.iloc[i]['ID']==gtf_SE.iloc[j]['ID']):
            #if((df_SE.iloc[i]['ID']==gtf_SE.iloc[j]['ID']) & (df_SE.iloc[i]['IJC_SAMPLE_1']>10) & \
               #(df_SE.iloc[i]['SJC_SAMPLE_1']>10) & (df_SE.iloc[i]['IncLevel1']<0.9) & (df_SE.iloc[i]['IncLevel1']>0.1)):
                ID.append(df_SE.iloc[i]['ID'])
                IJC.append(df_SE.iloc[i]['IJC_SAMPLE_1'])
                SJC.append(df_SE.iloc[i]['SJC_SAMPLE_1'])
                IncLen.append(df_SE.iloc[i]['IncFormLen'])
                SkipLen.append(df_SE.iloc[i]['SkipFormLen'])
                IncLevel.append(df_SE.iloc[i]['IncLevel1'])
                tmp=gtf_SE.iloc[j]['GeneID']
                gene_id.append(tmp[0:-1])
                chromosome.append(gtf_SE.iloc[j]['chr'])
                strand.append(gtf_SE.iloc[j]['strand'])
                exon_start.append(int(gtf_SE.iloc[j][5])+1)
                exon_end.append(int(gtf_SE.iloc[j][6]))
                upES.append(int(gtf_SE.iloc[j][7])+1)
                upEE.append(int(gtf_SE.iloc[j][8]))
                downES.append(int(gtf_SE.iloc[j][9])+1)
                downEE.append(int(gtf_SE.iloc[j][10]))
            while(df_SE.iloc[i]['ID']<gtf_SE.iloc[j]['ID']):
                i-=1
        
        merge={'ID':ID,'IJC':IJC,'SJC':SJC,'IncLen':IncLen,'SkipLen':SkipLen,'IncLevel':IncLevel,\
               'GeneID':gene_id,'chr':chromosome,'strand':strand,'exonStart':exon_start,'exonEnd':exon_end,\
               'upstreamES':upES,'upstreamEE':upEE,'downstreamES':downES,'downstreamEE':downEE}
        
        output=pd.DataFrame(merge,columns=['ID','IJC','SJC','IncLen','SkipLen','IncLevel',\
               'GeneID','chr','strand','exonStart','exonEnd',\
               'upstreamES','upstreamEE','downstreamES','downstreamEE'])
        print(output.shape[0])
        output.to_csv('../data/SW480/SW480_as_bam/'+se_type+'_merge_no_filter.tsv',sep='\t',index=None)
        return output
    
        
        
        
    else:
        print("Wrong Type name! The type name should be one of the SE,RI,A3SS,A5SS,MXE.")
        
    
                
            