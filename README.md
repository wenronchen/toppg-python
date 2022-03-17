# toppg-python

## Generate the customized protein database for MS identification

#### 0. Prerequisite
Python 3.X and corresponding version of Pandas, Biopython, Gffutils and  Bedtools package.

#### 1. Parameters

  - -f [ --file ], require the file name of ANNOVAR variants annotation with file extension "exonic_variant_function".
  - -f1 [ --file1 ], require the name of one output file of rMATS in "AS_Event.MATS.JCEC.txt" format.
  - -f2 [ --file2 ], require the name of one output file of rMATS in "fromGTF.AS_Event.txt" format.
  - -g [ --gff ], require the annotation file name in gff3 format.
  - -r [ --rna ],  require the fasta file of reference transcript sequences.
  - -o [ --output ],  the output name for database, default="customized_db".
  - -t [ --het ] <0|1|2>, the number of heterozygous genetic variants, default=1.
  - -s [ --splicing ], add the splicing variations to sequences.
  - -e [ --exclude ], exclude the sequences without any variations


##### 2. Examples

Generate the customized database with at most one heterozyous variant per sequence.

    python generate_db_pipeline.py -f DLD.exonic_variant_function -g gencode.v28.basic.annotation.gff3 -r gencode.v28.transcripts.basic.fasta -h 1 

  
Generate the customized database with splicing variations(Exon-skipping event) and at most one heterozyous variant per sequence.
   
      python generate_db_pipeline.py -f1 SE.MATS.JCEC.txt -f2 fromGTF.SE.txt -f DLD.exonic_variant_function -g gencode.v28.basic.annotation.gff3 -r gencode.v28.transcripts.basic.fasta -h 1 -s
   
   

