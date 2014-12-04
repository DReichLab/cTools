# cTools

This document is the instruction of using the cTools to access the light
version (compressed) Simon's genomic diversity project (SGDP lite) files. The
cTools include:

1. [cuncompress](#cuncompress)
2. cascertain
3. cpulldown
4. cpoly

SGDP lite includes the full C team data and the following high quality ancient
genomes (extended C-team):

               Altai  F            Altai
            Denisova  F        Denisovan
           Loschbour  M              WHG
           Stuttgart  F              LBK
           Ust_Ishim  M        Ust_Ishim

Currently, there are 340 samples in SGDP lite:

- FullyPublicGeneral: 163 samples (including the above ancient genomes)
- FullyPublicHGDP: 120 samples
- SignedLetterGeneral: 9 samples
- SignedLetterDiRienzo: 4 samples
- SignedLetterTishkoff: 44 samples

Size of the data:

- The size of the downloaded data (compressed):
  - One sample: ~379M
  - All the public data (183 samples in FullyPublicGeneral and FullyPublicHGDP): ~108G
  - All the 340 samples: ~130G


- The size of the data after uncompress:
  - One sample: ~6G
  - All the public data (183 samples in FullyPublicGeneral and FullyPublicHGDP): ~1,700G
  - All the 340 samples: ~2,040G

##How to download the SGDP lite sample files?
##Compile the source code
##Configuration (Only needed when you customized the architecture of your downloaded file folders.)
##How to access the data using cTools?

To access the sample files using cascertain, cpulldown and cpoly, these
downloaded sample files (`*.ccomp.fa.gz` and `*.ccompmask.fa.gz`) need to be
uncompressed using cuncompress first. Moreover, Samtools is needed (accessable
by command "Samtools") to run cuncompress properly.

#### <a name="cuncompress"></a>1. cuncompress -i sgdpsampname [options]

cuncompress uncompress the `*.ccomp.fa.gz` and `*.ccompmask.fa.gz` files.
```
options:
-w working dir / output dir [default: ./]
-? Show the instruction.

Inputs: .ccomp.fa.gz, .ccompmask.fa.gz
Outputs: .hetfa (fasta format), hetfa.fai, .mask.fa, .mask.fa.fai
```

***Notes***: The compressed input files are required in the workdir to allow
cuncompress work properly. The output files are in the workdir. Samtools is
required by cuncompress.

#### 2. cascertain
```
Usage: cascertain -p parfile [options]

Options: -V    verbose
         -v    show version information.
         -?    show the instruction.

Inputs: criteria to ascertain SNP.
Outputs: .snp file (Reich lab format)
```

Sample parfile:

```
snpname: ssout.snp
ascertain: S_Yoruba-1::1:2,Altai::1:1;S_Yoruba-1::1:2,Denisova::1:1
noascertain: Altai::1:2,Denisova::1:2
```
In the above sample parfile,
snpname: the output .snp file.

ascertain: There is a mini-language for ascertainment.  Allele 1 = derived, 0 =
ancestral (where Chimp allele is ancestral) In the ascertainment string we have
substrings separated by ';'  each is an ascertainment rule.  Then each rule has
substrings separated by ',' ; Each substring is of form `Sample_ID :: a : b`
This means we require sample to have a derived allele(s) out of b. a is the
number of derived alleles; b is the total number of alleles in the sample.

noascertain: has the same syntax.   If this "hits" then the SNP is rejected.
Thus the example above means ascertain if `S_Yoruba_1` is a het and either
Altai or Denisova has a derived allele (chosen at random) EXCEPT don't
ascertain if both Altai and Denisova are hets.

A full parameter list of the parfile:

| parameter     | description |
|---------|--------------------------------------------------|
| regname       | chromosome name [default: all the chromosomes] |
| lopos         | the beginning coordinate of the region [default: the beginning of the chromosome] |
| hipos         | the ending coordinate of the region [default: the end of the chromosome] |
| snpname       | the output snp file name (.snp) |
| minfilterval  | the base quality threshold for taking the genotype information; The,quality values are in the mask file. C team has base quality in range,(0-9) or no value (N/?) => don't use. Select bases with minfilterval: 3 (say). This selects bases with base quality >=3. [1 is default and,recommended for most applications]. Note that the extended C-team files,such as Altai have manifesto filters (made in Leipzig) that are just 0, 1. If you are using extended C-team do not set minfilterval. |
| ascertain     | ascertain criteria  |
| noascertain   | rejection criteria  |
| dbhetfa       | table of hetfa files (.dblist); [default: `/home/mz128/cteam/dblist/,hetfa_postmigration.dblist`]  |
| dbmask        | table of mask files (.dblist); [default: `/home/mz128/cteam/dblist/,mask_postmigration.dblis`]  |
| transitions   | work on transitions; [default: Yes]  |
| transversions | work on transversions; [default: Yes] |
|  minchrom:	| the beginning chromosome |
|  maxchrom:	| the ending chromosome |

Notes: You can write comments in the parameter file by `#comments`.

#### 3. cpulldown  -p parfile [options]

options:
 -V	verbose
 -c	checkmode
 -v	Show version information.
 -? 	Show the instruction.


Inputs: .snp file, .ind file (Reich lab format)

Sample parfile

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
This is very similar to the parameter file for convertf which most of you will
have used.  indivname should be a .ind file (Reich lab format). Any .snp file
is OK here, it needn't be output from cascertain.


A full parameter list of the parfile:

| parameter       | description        |
|-----------------|--------------------|
| snpname         | input .snp file    |
| indivname       | input .ind file    |
| indivoutname    | output .ind file   |
| snpoutname      | output .snp file   |
| genotypeoutname | output .geno file  |
| outputformat    | output format      |
| minfilterval    | similar to the minfilterval of  cascertain, but [default is 0] |
|  minchrom:	| the beginning chromosome
|  maxchrom:	| the ending chromosome



Notes: Running time: Linear in number of samples in indivname.  10 samples pull down
in about 1 hour on orchestra.

#### 4. cpoly
DOCUMENTATION NEEDED

##Related file formats

<!--
Written by Nick on 6/15/14
Revised by Mengyao Zhao
Last revision: 11/26/14
-->
