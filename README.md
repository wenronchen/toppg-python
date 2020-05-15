# toppg-python

## Generate the customized protein database for MS identification

#### 0. Prerequisite
Python 3.X and corresponding version of Pandas, Biopython package.

#### 1. Protein database with mutations

##### Required input parameters
  - -f [ --file ], The output file from ANNOVAR,for example,ex1.exonic_variant_function
  - -g [ --gff ], The annotation file in gff3 format.
  - -s [ --sequence ], The fasta file of protein sequences.
  - -t [ --type ], The type of mutations when generating database, either f(frameshift and non-frameshift) or p(point mutations or c(the protein database with combinations).

##### 1.1 Generate protein database with point mutations(subtitutions)

You need to choose p as the parameter of type to generate protein sequences only including point mutations.
  
For example,

    python generate_db_pipeline.py -f sample.exonic_variant_function -g ../data/gencode.v28.basic.annotation.gff3 -s gencode.v28.basic.fasta -t p


##### 1.2 Generate protein database with frameshift and non-frameshift mutations

You need to choose f as the parameter of type to generate protein sequences only including frameshift or non-frameshift mutations.
Besides, when choosing f as the parameter, you should add the following parameter:
  - -r [ --rna ], The fasta file of transcript sequences

For example,

    python generate_db_pipeline.py -f sample.exonic_variant_function -g gencode.v28.basic.annotation.gff3 -s gencode.v28.basic.fasta -t f -r gencode.v28.pc_transcripts.fa


##### 1.3 Generate protein database including sequences with combinations

The protein sequences with combinations have 6 types:
  - Point mutation and point mutation
  - Point mutation and non-frameshift mutation
  - Point mutation and frameshift mutation
  - Non-frameshift mutation and non-frameshift mutation
  - Non-frameshift mutation and point mutation
  - Non-frameshift mutation and frameshift mutation
  
 You need to use c as the parameter of type to generate protein sequences including above 6 types of combinations.
 Besides, you need to add the r parameter
 
 For example,
 
    python generate_db_pipeline.py -f sample.exonic_variant_function -g gencode.v28.basic.annotation.gff3 -s gencode.v28.basic.fasta -t c -r gencode.v28.pc_transcripts.fa
 
 
  
  
 #### 2. Protein database with splicing events
 
  You should have the output of rMATS and download the relative gff file.
  #####  Required parameters
  - -f1 [ --file1 ], The output file from rMATS,for example,SE.MATS.JCEC.txt.
  - -f2 [ --file2 ], The output file from rMATS, for example, fromGTF.SE.txt.
  - -g [ --gff ], The annotation file in gff3 format.
  - -t [ --type ], The type of splicing events, like SE, RI.
   
   For example, if you want to generate a database with retained intron(RI) splicing events,
   
      python generate_se_db_pipeline.py -f1 RI.MATS.JCEC.txt -f2 fromGTF.RI.txt -g gencode.v28.basic.annotation.gff3 -t RI
   
   

