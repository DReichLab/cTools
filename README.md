# cTools

This document is the instruction of using the cTools to access the light
version (compressed) Simons genomic diversity project (SGDP lite) files. The
cTools include:

1. uncompress.pl
2. cascertain
3. cpulldown
4. cpoly

##Download the programs
`git clone https://github.com/mengyao/cTools.git`

or click the "Download ZIP" button of this web page:
https://github.com/mengyao/cTools

##Compile cTools
In the src folder:

1. `make clean`
2. `make`
3. `make install`

##Run cTools

To access the sample files using cascertain, cpulldown and cpoly, these
downloaded sample files (`*.ccomp.fa.gz` and `*.ccompmask.fa.gz`) need to be
uncompressed using uncompress.pl first. Moreover, Samtools is needed (accessable
by command "Samtools") to run uncompress.pl properly.

### <a name="uncompress.pl"></a>1. uncompress.pl

uncompress.pl uncompress the `*.ccomp.fa.gz` and `*.ccompmask.fa.gz` files.
```
Usage: uncompress.pl <reference.fa> <sample.ccomp.fa.gz> (compressed hetfa file)

Inputs: .ccomp.fa.gz, .ccompmask.fa.gz
Outputs: .fa (uncompressed hetfa file), .fa.fai, .filter.fa (uncompressed mask file), .filter.fa.fai
```

***Notes***: The  Samtools is
required by uncompress.pl.

### 2. cascertain
cascertain pulls down the SNPs that match the ascertain rule(s).

```
Usage:   cascertain -p <parameter file> [options]

Options:
  -d	directory of the data files (Please set this parameter, if you do not set .dblist files. If this parameter is used to give the data file location, .dblist files will not be used.)
  -V	Print more information while the program is running.
  -v	Show version information.
  -?    show the instruction.

Input: parameter_file
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

A full parameter list of the parfile:

| parameter     | description |
|---------|--------------------------------------------------|
| chrom       | chromosome name; The value can be one of [1 - 22, X]. If chrom is set, please do not set minchrom or maxchrom. [default: all the chromosomes] |
| lopos         | the beginning coordinate of the region [default: the beginning of the chromosome] |
| hipos         | the ending coordinate of the region [default: the end of the chromosome] |
| snpname       | the output snp file name (.snp) |
| minfilterval  | the base quality threshold for taking the genotype information; The,quality values are in the mask file. C team has base quality in range,(0-9) or no value (N/?) => don't use. Select bases with minfilterval: 3 (say). This selects bases with base quality >=3. [1 is default and,recommended for most applications]. Note that the extended C-team files, such as Altai have manifesto filters (made in Leipzig) that are just 0, 1. If you are using extended C-team do not set minfilterval. |
| ascertain     | If the SNP matches the rule(s) here, this SNP will be out put.  |
| noascertain   | If the SNP matches the rule(s) here, it will NOT be out put.  |
| dbhetfa       | .dblist file that specify the hetfa file, refrence.fa and chimp.fa location. If the -d option is not used, this parameter is mandatory.  |
| dbmask        | .dblist file that specify the mask file location. If the -d option is not used, this parameter is mandatory.|
| transitions   | work on transitions; [default: Yes]  |
| transversions | work on transversions; [default: Yes] |
|  minchrom:	| the beginning chromosome; If minchrom and maxchrom are set, please do not set chrom. |
|  maxchrom:	| the ending chromosome |
|pagesize|cascertain "pages" through the genome in chunks of size pagesize bases. The default is 20M bases, but this can be overwritten.  Larger pages will run faster (and use more memory).|
|seed|For hets and some ascertainments a random allele must be chosen.  This is picked be a pseudo-random generator.  Set seed:  SEED where SEED is a generator if you wish the runs to be reproducible.|

Notes: You can write comments in the parameter file by `#comments`.

### 3. cpulldown
cpulldown out puts the genotypes of the given individuls at the given SNP loci from a set of bams.

```
Usage:   cpulldown -p <parameter file> [options]

Options:
  -d  directory of the data files (Please set this parameter, if you do not set .dblist files. If this parameter is used to give the data file location, .dblist files will not be used.)
  -V  Print more information while the program is running.
  -c  checkmode
  -v  Show version information.
  -?  Show the instruction.


Inputs: parameter_file, .snp file, .ind file (Reich lab format)
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
|  minchrom:	| the beginning chromosome
|  maxchrom:	| the ending chromosome


For example, to use cpulldown to pull out data from Papuan and Dai, you need to:

1. Uncompress all the hetfa and mask files of Papuan and Dai samples, and put them into one folder together with Chimp.fa, Chimp.fa.fai, hs37d5.fa and hs37d5.fa.fai.

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

### 4. cpoly
cpoly pulls down all the SNP sites that are polymorphic in sample list from one or multiple bams. In default mode only bases with no  missing
data are considered, so you need to be careful if you use many samples.

```
Usage:   cpoly -p <parameter file> [options]

Options:
  -d	directory of the data files (Please set this parameter, if you do not set .dblist files. If this parameter is used to give the data file location, .dblist files will not be used.)
  -V	Print more information while the program is running.
  -v	Show version information.
  -? 	Show the instruction.

Inputs: parameter_file, .ind file
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
| dbhetfa       | .dblist file that specify the hetfa file, refrence.fa and chimp.fa location. If the -d option is not used, this parameter is mandatory.  |
| dbmask        | .dblist file that specify the mask file location. If the -d option is not used, this parameter is mandatory.|
|pagesize|cascertain "pages" through the genome in chunks of size pagesize bases. The default is 20M bases, but this can be overwritten.  Larger pages will run faster (and use more memory).|
| transitions   | work on transitions; [default: Yes]  |
| transversions | work on transversions; [default: Yes] |
| allowmissing| The output may contain sites that for some samples the data are missing.|
|polarize| The polarize parameter is optional.  If present the parameter should be a sample name present in the indivname file. Then only homozygotes of this sample are considered, and the first allele of any snp is the base for the poliarize sample. As an example: "polarize: Href", the first allele of every snp in snpoutname will be the Href allele. Usually the polarize sample will be a pseudo-diploid such as Href or Chimp. |

##Related file formats
###1. .ind file
example.ind:
```
#FullyPublicGeneral
Altai   F   Altai
Denisova    F   Denisovan
Loschbour   M   WHG
Stuttgart   F   LBK
```
- Every thing after '#' is comment.
- 1st column: sample name
- 2nd column: gender
- 3rd column: population

###2. .snp file
example.snp:
```
X:21_15838960    21        0.158390        15838960 C T
X:21_16081858    21        0.160819        16081858 T C
X:21_16495307    21        0.164953        16495307 G A
...
```
- 1st column: name of the SNP site (if the SNP site is in DBSNP: DBSNP name; else: arbitrary name)
- 2nd column: chromosome name
- 3rd column: genetic location
- 4th column: coordinate
- 5th column: base allele
- 6th column: derived allele

###3. .geno file
example.geno
```
00002
00002
22220
...
```
- .geno file is used together with .ind file and .snp file.
- The count of columns in .geno file is the count of individuals in .ind file. From left to right, each column of .geno corresponds to one individual (from top to bottom) in .ind file.
- The count of lines in .geno file is the count of SNP sites in .snp file. Each line of .geno corresponds to one line in .snp file.
- Each number in .geno denotes the count of base alleles (the allele of the 5th column of .snp file) of the corresponding individual at the corresponding SNP site.

###4. .dblist

hetfa.example.dblist:
```
Href    -   /home/mz128/cteam/usr/data/Href.fa
Chimp   -   /home/mz128/cteam/usr/data/Chimp.fa
S_Albanian-1    -   /home/mz128/cteam/usr/data/S_Albanian-1.fa
S_Dai-1 -   /home/mz128/cteam/usr/data/S_Dai-1.fa
S_Yoruba-1  -   /home/mz128/cteam/usr/data/S_Yoruba-1.fa
```
- 1st column: sample name
- 2nd column: -
- 3rd column: location of the hetfa files (/absolute_path/sample.fa)

mask.example.dblist:
```
Href    -   NULL
Chimp   -   NULL
S_Albanian-1    -   /home/mz128/cteam/usr/data/S_Albanian-1.filter.fa
S_Dai-1 -   /home/mz128/cteam/usr/data/S_Dai-1.filter.fa
S_Yoruba-1  -   /home/mz128/cteam/usr/data/S_Yoruba-1.filter.fa
```
- The 1st and 2nd columns are the same as the .dblist file for hetfa files.
- For chimp and human references - 3rd column: NULL
- For samples - 3rd column: location of the mask files (/absolute_path/sample.filter.fa)

<!--
Written by Nick on 6/15/14
Revised by Mengyao Zhao
Last revision: 12/08/14
-->
