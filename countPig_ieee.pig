/*loading a BZip'ed SAM archive*/
almnts= LOAD '/user/hdfs/workspace/raalesir/sam/schneb_IV_half.sam_noconcat1.bz2' as (  name_almnt:chararray,  flag:int, chr_almnt:chararray, start_almnt:int, f5:chararray, f6:chararray, f7:chararray, end_almnt:int, f9:chararray, f10:chararray,f11:chararray,f12:chararray,f13:chararray, f14:chararray);


/*filtering out chromosomes of interest*/
almnts = FILTER almnts BY (chr_almnt matches 'Chr.*') PARALLEL 54;
/*filtering out pair-ended reads properly mapped to the '+' and '-' strands*/
almnts_forward  = FILTER almnts BY (flag == 99)  PARALLEL 54;
almnts_reverse =  FILTER almnts BY (flag == 163)  PARALLEL 54;
   
/*adding '+' and '-' columns in the read tables*/
almnts_forward = FOREACH almnts_forward GENERATE chr_almnt, (chararray)'+' as strand_almnt, start_almnt, end_almnt+100 as end_almnt, name_almnt PARALLEL 54;
almnts_reverse = FOREACH almnts_reverse GENERATE chr_almnt, (chararray)'-' as strand_almnt, start_almnt, end_almnt+100 as end_almnt, name_almnt PARALLEL 54;

/*union the tables*/
almnts = UNION almnts_forward, almnts_reverse PARALLEL 54;


/*load the gene annotation file in GFF format*/
genes = LOAD '/user/hdfs/workspace/raalesir/TAIR10_GFF3_genes.gff ' USING PigStorage('\t') AS (chr_f:chararray, tmp:chararray, feature:chararray,
start_f:int, end_f:int, tmp1:chararray, strand_f:chararray, tmp2:chararray, name_f:chararray);  

/*filter the features of interest and extract the feature name*/
genes = FILTER genes BY (feature == 'gene') PARALLEL 54;
genes = FOREACH genes GENERATE chr_f, start_f, end_f, strand_f, REGEX_EXTRACT(name_f, '(AT[0-9][A-Z][0-9]*)',1) as name_f;


/* registering a small Python function*/
register 'myudf.py' using jython as myudf ;


/* partition genes into groups, each of which oppucies a stretch of a 10K*/
genesWithBox = foreach genes generate flatten(myudf.foo(chr_f, strand_f, start_f, end_f, name_f, 10000));     

/* joining the tables for genes and mapped reads by the chromosome and left
 * border of the partition interval*/
cogrp = COGROUP genesWithBox  BY (chr_f, start_f/10000+box) , almnts BY (chr_almnt, start_almnt/10000) PARALLEL 54;

/*for each group create cross product of gene-reads couples, so each read in
 * the group will be coupled with a each gene form the same group */
test = FOREACH cogrp GENERATE FLATTEN(genesWithBox.(name_f, strand_f, start_f, end_f)), FLATTEN(almnts.(name_almnt, strand_almnt, start_almnt, end_almnt)) PARALLEL 54;

/* sort the resulting after the cross product table by the read name */
grpAlmnt = GROUP test BY name_almnt PARALLEL 54;

/*for each read in each group filter those genes to which the read is mapped to, and count the genes */
countRead = FOREACH grpAlmnt {   fltt = FILTER test BY (start_almnt >= start_f ) AND (end_almnt < end_f) AND (strand_f == strand_almnt); 
GENERATE group, COUNT (fltt) as count, flatten(fltt.name_f);} ;  

/*filter those reads which are covered just by one gene*/
countRead = FILTER countRead BY $1<2;                  

/* group the data by the gene name*/                                                                                        
grpGene = GROUP countRead BY $2 PARALLEL 54;

/* count the reads mapped to each gene */ 
countGene = foreach grpGene generate group, COUNT(countRead); 

STORE countGene  INTO '/user/hdfs/workspace/raalesir/countGeneIV_half/3';  
 
