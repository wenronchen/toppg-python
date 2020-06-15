# toppg-python

## Generate the customized protein database for MS identification

#### 0. Prerequisite
Python 3.X and corresponding version of Pandas, Biopython, Gffutils and  Bedtools package.

#### 1. Protein database with mutations

##### 1.1 Required input parameters
  - -f [ --file ], the output file from ANNOVAR,for example,ex1.exonic_variant_function
  - -g [ --gff ], the annotation file in gff3 format.
  - -r [ --rna ],  the fasta file of reference transcript sequences.
  - -s [ --sequence ], the fasta file of protein sequences.
  - -n [ --name ], the name of dataset.

##### 1.2 Output

A fasta file contains all generated customized protein sequences with mutations.

##### 1.3 Example

    python generate_db_pipeline.py -f DLD.exonic_variant_function -g gencode.v28.basic.annotation.gff3 -s gencode.v28.basic.fasta -r gencode.v28.pc_transcripts.fa -n DLD


  
  
 #### 2. Protein database with splicing events
 
  You should have the output of rMATS and download the relative gff file.
  #####  2.1 Required parameters
  - -f1 [ --file1 ], the output file from rMATS,for example,SE.MATS.JCEC.txt.
  - -f2 [ --file2 ], the output file from rMATS, for example, fromGTF.SE.txt.
  - -g [ --gff ], the annotation file in gff3 format.
  - -r [ --rna ],  the fasta file of reference transcript sequences.
  - -n [ --name ], the name of dataset.
  - -t [ --type ], The type of splicing events.
   
##### 2.2 Output

A fasta file contains all generated customized protein sequences with splicing events.

##### 2.3 Example
   
      python python generate_se_db_pipeline.py -f1 X.MATS.JCEC.txt -f2 fromGTF.X.txt -g gencode.v28.basic.annotation.gff3 -r gencode.v28.pc_transcripts.fa -n DLD -t X
   
   

