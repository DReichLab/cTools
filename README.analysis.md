# How to run cteam-lite analysis tools

This document describes how to run cascertain, cpulldown and cpoly. All three tools are designed to run on compressed .hetfa and .mask files directly.


##.-1  OVERVIEW

Download instructions for the main dataset is in Section ## 0. Then, the following tools may be used:

cascertain: pulls down the SNPs that match the ascertain rule(s) (Section: ## 1)
cpulldown out puts the genotypes of the given individuls at the given SNP loci from a set of hetfa files (see Section: ## 2)
cpoly pulls down all the SNP sites that are polymorphic in sample list from one or multiple hetfas. In default mode only bases with no  missing (see Section: ## 3).




##.0  GENERAL 
Please do the following:

1) Download/extract the SGDP dataset. This contains the sgdp hetfa and mask files. Download the tar file, and extract using tar xvf <tar>. Put this in a convenient place in your file hierarchy.
2) Modify files: .../filter/hetfa.dblist and .../filter/mask.dblist to point to your copy of the SGDP  
3) For all 3 programs the parameter file should have a fragmemt like: 

DIR:           /home/np29/biology/sgdp/info 
dbhetfa:       DIR/hetfa.dblist
dbmask:        DIR/mask.dblist

Format of the .dblist files are 3 columns.  The first column (sample name) must match 
the sample ID in the .ind input file.  Last columns is full path name.
For the mask.dblist NULL is valid => no mask.  


## 1. cascertain
cascertain pulls down the SNPs that match the ascertain rule(s).

```
Usage:   cascertain -p <parameter file> [options]

Options:
  -V	Print more information while the program is running.
  -v	Show version information.
  -?    show the instruction.
*** pointers to hetfa.dblist and mask.dblist as above  

Input: hetfa files, parameter file
Outputs: .snp file (Reich lab format)
```

Sample parameter_file:

```
snpname: ssout.snp
ascertain: S_Yoruba-1::1:2,Altai::1:1;S_Yoruba-1::1:2,Denisova::1:1
noascertain: Altai::1:2,Denisova::1:2
```
In the above sample parameter_file,

snpname: the output .snp file.

ascertain: There is a mini-language for ascertainment.  Allele 1 = derived, 0 =
ancestral (where Chimp allele is ancestral) In the ascertainment string we have
substrings separated by ';'  each is an ascertainment rule.  Then each rule has
substrings separated by ',' ; Each substring is of form `Sample_ID :: a : b`
This means we require sample to have a derived allele(s) out of b. a is the
number of derived alleles; b is the total number of alleles in the sample.

noascertain: has the same syntax.   If this "hits" then the SNP is rejected.
Thus the example above means ascertain (out put the SNP) if `S_Yoruba_1` is a het and either
Altai or Denisova has a derived allele (chosen at random) EXCEPT don't
ascertain if both Altai and Denisova are hets.

Assume the Chimp allele is the ancestral allele, here are some examples of the acertainment marks:

||		Chimp|	Human| note |
|----------|---------|--------|-----|
|0:1	|	A	|	A | The sample allele is randomly picked, and this allele is ancestral. |
|1:1	|	A	|	T | The sample allele is randomly picked, and this allele is derived. |
|1:2	|	A	|	A/T | Alleles on both chromosomes are picked, and one of them is derived. |
|2:2 | A | T/T | Alleles on both chromosomes are picked, and both of them are derived. |

Note: In cascertain we ignore bases which are not biallelic in the samples.

A full parameter list of the parameter_file:

| parameter     | description |
|---------|--------------------------------------------------|
| chrom       | chromosome name; The value can be one of [1 - 22, X]. If chrom is set, please do not set minchrom or maxchrom. [default: all the chromosomes] |
| lopos         | the beginning coordinate of the region [default: the beginning of the chromosome] |
| hipos         | the ending coordinate of the region [default: the end of the chromosome] |
| snpname       | the output snp file name (.snp) |
| minfilterval  | the base quality threshold for taking the genotype information; The quality values are in the mask file. C team has base quality in range (0-9) or no value (N/?) => don't use. Select bases with minfilterval: 3 (say). This selects bases with base quality >=3. [1 is default and recommended for most applications]. Note that the extended C-team files, such as Altai have manifesto filters (made in Leipzig) that are just 0, 1. If you are using such files, do not set minfilterval greater than 1, as all data will be masked out. |
| ascertain     | If the SNP matches the rule(s) here, this SNP will be out put.  |
| noascertain   | If the SNP matches the rule(s) here, it will NOT be out put.  |
| dbhetfa       | .dblist file that specify the hetfa file, refrence.fa and chimp.fa location. This parameter is mandatory for users not in Reich Lab.  |
| dbmask        | .dblist file that specify the mask file location. Tthis parameter is mandatory for users not in Reich Lab.|
| transitions   | work on transitions. [default: Yes]  |
| transversions | work on transversions. [default: Yes] |
|  minchrom:	| the beginning chromosome; If minchrom and maxchrom are set, please do not set chrom. |
|  maxchrom:	| the ending chromosome |
|pagesize|cascertain "pages" through the genome in chunks of size pagesize bases. The default is 20M bases, but this can be overwritten.  Larger pages will run faster (and use more memory). Please use a number only to set this parameter. |
|seed|For hets and some ascertainments a random allele must be chosen.  This is picked by a pseudo-random generator.  Set seed:  SEED where SEED is a generator if you wish the runs to be reproducible.|

Notes: You can write comments in the parameter file like this:  # this is a comment  

## 2. cpulldown
cpulldown out puts the genotypes of the given individuls at the given SNP loci from a set of hetfa files.

```
Usage:   cpulldown -p <parameter file> [options]

Options:
  -V  Print more information while the program is running.
  -c  checkmode
  -v  Show version information.
  -?  Show the instruction.


Inputs: hetfa files, parameter_file, .snp file, .ind file (Reich lab format)
Outputs: .snp file, .ind file, .geno file

```

Sample parameter_file1:

```
indivname:    mix1.ind
snpname:      input.snp
indivoutname:     output.ind
snpoutname:       output.snp
genotypeoutname:  output.geno
outputformat:     eigenstrat
maxchrom:         22
```

Sample parameter_file2:

```
D1:          /home/np29/biology/cteam/mixdir
D2:          D1
S2:           sc
indivname:    D1/mix1.ind
snpname:      D1/s01:01.snp
indivoutname:     D2/S2.ind
snpoutname:       D2/S2.snp
genotypeoutname:  D2/S2.geno
outputformat:     eigenstrat
maxchrom:         22
```
indivname should be a .ind file (Reich lab format). Any .snp file
is OK here, it needn't be output from cascertain. In "sample parameter_file2", D1, D2 and S2 are symbols to replace long strings. For example, here D1 = D2 = /home/np29/biology/cteam/mixdir; S2 = sc.


A full parameter list of the parameter_file:

| parameter       | description        |
|-----------------|--------------------|
| snpname         | input .snp file    |
| indivname       | input .ind file    |
| indivoutname    | output .ind file   |
| snpoutname      | output .snp file   |
| genooutname | output .geno file  |
| outputformat    | [eigenstrat]      |
| minfilterval    | similar to the minfilterval of  cascertain, but [default is 0] |
|  minchrom:	| the beginning chromosome |
|  maxchrom:	| the ending chromosome |
| dbhetfa       | .dblist file that specify the hetfa file, refrence.fa and chimp.fa location. This parameter is mandatory for the users not in Reich Lab.  |
| dbmask        | .dblist file that specify the mask file location. This parameter is mandatory for the users not in Reich Lab.|

For example, to use cpulldown to pull out data from Papuan and Dai, you need to:

1. Put all the hetfa and mask files of Papuan and Dai samples into one folder together with Chimp.fa, Chimp.fa.fai, hs37d5.fa and hs37d5.fa.fai.

2. Using the information from all.ind to create your input.ind as following:

```
B_Papuan-15 M Papuan
A_Papuan-16 M Papuan
S_Papuan-6 M Papuan
S_Papuan-11 M Papuan
S_Papuan-14 F Papuan
S_Papuan-4 M Papuan
S_Papuan-3 M Papuan
S_Papuan-8 M Papuan
S_Papuan-7 M Papuan
S_Papuan-9 M Papuan
S_Papuan-2 M Papuan
S_Papuan-12 M Papuan
S_Papuan-10 M Papuan
S_Papuan-5 M Papuan
S_Papuan-13 F Papuan
S_Papuan-1 F Papuan
B_Dai-4 M Dai
A_Dai-5 M Dai
S_Dai-1 F Dai
S_Dai-3 F Dai
S_Dai-2 M Dai

```

3. Create your input.snp file. The format of this .snp file is described in the "Related file formats" section of this document. If you don't know the "genetic location" of the SNP, please put value "0" in that column. In this condition, the "genetic location" column of the output.snp file will NOT be meaningful either.  
4. Create a parameter file example.par as following:

```
indivname:    input.ind
snpname:      input.snp
indivoutname:     output.ind
snpoutname:       output.snp
genotypeoutname:  output.geno
outputformat:     eigenstrat
```

4. `cpulldown -p example.par -d <dir of your hetfa & mask files>`

Notes: Running time: Linear in number of samples in indivname.  10 samples pull down
in about 1 hour on orchestra.

## 3. cpoly
cpoly pulls down all the SNP sites that are polymorphic in sample list from one or multiple hetfas. In default mode only bases with no  missing
data are considered, so you need to be careful if you use many samples.

```
Usage:   cpoly -p <parameter file> [options]

Options:
  -V	Print more information while the program is running.
  -v	Show version information.
  -? 	Show the instruction.

Inputs: hetfa files, parameter_file, .ind file
Outputs:  .ind file, .snp file, .geno file
```

Sample parameter_file:

```
D1:          /home/mz128/cteam/mixdir
D2:          .
S2:           qpoly2
indivname:       D1/mix1.ind
indivoutname:     D2/S2.ind
snpoutname:       D2/S2.snp
genooutname:     D2/S2.geno
##polarize:         Href
## If this is used Href must be present in indivname file
##chrom:            24
##lopos:            46300000
##hipos:            46610000
minfilterval:       1
allowmissing:          NO
## default YES
allowhets:          NO
## transversions:   YES
## transitions:     YES
## dbhetfa:
## dbmask:
```

A full parameter list of the parameter_file:

| parameter       | description        |
|-----------------|--------------------|
| indivname       | input .ind file    |
| indivoutname    | output .ind file   |
| snpoutname      | output .snp file   |
| genooutname | output .geno file  |
| minfilterval    | similar to the minfilterval of  cascertain, but [default is 0] |
| chrom       | chromosome name; The value can be one of [1 - 22, X, Y, MT]. If chrom is set, please do not set minchrom or maxchrom. [default: all the chromosomes] |
|  minchrom:	| the beginning chromosome; If minchrom and maxchrom are set, please do not set chrom. |
|  maxchrom:	| the ending chromosome|
| lopos         | the beginning coordinate of the region [default: the beginning of the chromosome] |
| hipos         | the ending coordinate of the region [default: the end of the chromosome] |
| dbhetfa       | .dblist file that specify the hetfa file, refrence.fa and chimp.fa location. This parameter is mandatory for the users not in Reich Lab.  |
| dbmask        | .dblist file that specify the mask file location. This parameter is mandatory for the users not in Reich Lab.|
|pagesize|cascertain "pages" through the genome in chunks of size pagesize bases. The default is 20M bases, but this can be overwritten.  Larger pages will run faster (and use more memory). Please use one number only for this parameter setting.|
| transitions   | work on transitions; [default: NO]  |
| transversions | work on transversions; [default: NO] |
| allowmissing| The output may contain sites that for some samples the data are missing.|
|polarize| The polarize parameter is optional.  If present the parameter should be a sample name present in the indivname file. Then only homozygotes of this sample are considered, and the first allele of any snp is the base for the poliarize sample. As an example: "polarize: Href", the first allele of every snp in snpoutname will be the Href allele. Usually the polarize sample will be a pseudo-diploid such as Href or Chimp. |

<!--
Written by Nick on 6/15/14
Revised by Mengyao Zhao and Nick
Last revision: 11/14/16  

-->
