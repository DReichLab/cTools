# cTools

This document is the instruction of using the cTools to access the light
version (compressed) Simons genomic diversity project (SGDP lite) files. The
cTools include:

Analysis tools:

1. cascertain
2. cpulldown
3. cpoly

Tools for handling samples:

<ol start="4">
  <li>uncompress.pl</li>
  <li>ccompress</li>
  <li>cmakefilter</li>
</ol>

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
Samtools may be needed (accessable by command "Samtools") to run cTools properly.

Please see how to run cascertain, cpulldown and cpoly in README.analysis.md; see how to run uncompress.pl, ccompress and cmakefilter in README.prepare.md.

## Related file formats
### 1. .ind file
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

### 2. .snp file
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
- 5th column: base allele. In the output of cascertain (with polarize to Href set), this is the Href allele.
- 6th column: alternative allele. In the output of cascertain, this is the forward strand allele.
- For cpulldown and cpoly, the column 5 and 6 of the output .snp file are the same as that of the input. They can be any user chosen alleles.

### 3. .geno file
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
- The number is the count of the allele in the 5th column of the .snp file. For example, if the .snp file is example.snp, for the 1st snp locus, '0' means homozygous T, '2' means homozygous C.

### 4. .dblist

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


----------------------------------------
*For questions please contact Mengyao Zhao (mengyao_zhao@hms.harvard.edu) or Nick Patterson (nickp@broadinstitute.org).*

<!--
Written by Nick on 6/15/14
Revised by Mengyao Zhao
Last revision: 10/22/15
-->
