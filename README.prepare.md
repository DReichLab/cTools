# The tools for preparing your own samples

Besides the data provided in cteam lite, you can add your own data for the genotype analysis. uncompress.pl, ccompress and cmakefilter can help you to prepare your own samples.

## 1. uncompress.pl

uncompress.pl decompresses the `.ccomp.fa.rz` and `.ccompmask.fa.rz` files. The outputs of uncompress.pl can be read by other cTools programs.
```
Usage: uncompress.pl <reference.fa> <sample.ccomp.fa.rz> (compressed hetfa file)

Inputs: Href.fa .ccomp.fa.rz, .ccompmask.fa.rz
Outputs: .fa (uncompressed hetfa file), .fa.fai, .mask.fa (uncompressed mask file), .mask.fa.fai
```

***Notes***: 
Samtools and htsbox are required by uncompress.pl. 
The input files sample.ccomp.fa.rz and sample.ccompmask.fa.rz need to be in the same folder. The output files will be in the same folder of the input files.  

## 2. ccompress

ccompress compresses the `.fa` (hetfa) and `.filter.fa` (mask) files. The outputs of ccompress can be read by other cTools programs.

```
Usage: ccompress <reference.fa> <sample.fa>

Inputs: Href.fa sample.fa, sample.mask.fa
Outputs: .ccomp.fa.rz (compressed hetfa file), .ccompmask.fa.rz (compressed mask file)
```

***Notes***: 
htsbox is required by ccompress.
The input files sample.fa and sample.mask.fa need to be in the same folder. The output files will be in the same folder of the input files.

## 3. cmakefilter

cmakefilter is used to generate .mask.fa files from hetfa files (.fa files of the samples). 

```
Usage: cmakefilter -p <parameter.par>

Inputs:
parameter file: parameter.par 
raw vcfs by chromosome: .vcf.gz 
raw vcf hetfa: rawvcf.hetfa.fa 
cnv filter for the sample: hs37d5_maskCNV_soft.fa
universal filter: x75.fa, hs37d5.fa, pt2__cs-hg19.fa

Outputs: sample.mask.fa, sample.mask.fai
```

***Notes***:
Samtools is required by cmakefilter.
cmakefilter requires 20G memory. For LSF users, this is an example command to run cmakefilter: bsub -W 36:00 -R "rusage[mem=20000]" -o example.out -e example.err cmakefilter -p example.par

Sample parameter_file:

```
vcfdir:   /home/mz128/group/Cteam/cteam_lite6/filter_example/rawvcf
vcfsuffix: short.vcf
hetfa: /home/mz128/group/Cteam/cteam_lite6/filter_example/rawvcf.hetfa.fa
cnv:  /home/mz128/group/Cteam/cteam_lite6/filter_example/hs37d5_maskCNV_soft.fa
gender: M
## gender MUST be specified
href:	/home/mz128/group/Cteam/cteam_lite6/filter_example/hs37d5.fa
chimp:	/home/mz128/group/Cteam/cteam_lite6/filter_example/pt2__cs-hg19.fa
heng75:	/home/mz128/group/Cteam/cteam_lite6/filter_example/x75.fa
maskname:   S_Aleut-1.mask.fa
```

A full parameter list of the parameter_file:

| parameter       | description        |
|-----------------|--------------------|
| cnv      | sample dependent CNV filter (CNVmask.fa, The absolute path of this file should be given here.)    |
| gender    | gender of the sample [M/F]   |
| vcfdir      | directory of the raw vcf files. The raw vcf files should be provided by chromosomes in .gz format. For example 22.vcf.gz.   |
| href | hs37d5.fa (The absolute path of this file should be given here.) |
| chimp | pt2__cs-hg19.fa (The absolute path of this file should be given here.) |
| heng75 | x75.fa (The absolute path of this file should be given here.) |
| vcfsuffix | suffix of the rawvcf files. If the raw vcf files are *.vcf.gz, you should put "vcf" here. |
| hetfa | rawvcf hetfa file (The absolute path of this file should be given here.) |
| fixeddbase | universal filter data base. If this parameter is given, "chimp", "heng75" and "href" will not be used. |
| sampname | sample name for naming the output files (If "mengyao" is given to this parameter, the output mask file name will be "mengyao.mask.fa".) |
| maskname | name of the output file (If this parameter is set, "sampname" will not be used.) |

## 4. filtstats

filtstats calculates the divergence between the input sample and the references based on the input genotype information.  

```
Usage: filtstats -p <parameter.par>

Inputs:
parameter file: parameter.par 
raw vcfs by chromosome: .vcf.gz 
raw vcf hetfa: rawvcf.hetfa.fa 
mask file: sample.mask.fa, sample.mask.fai
human reference genome: Href.fa
chimpanzee genome: chimp.fa

Outputs: statistics to the standard output 
```

The parameter file of filtstats for one sample should be the same as that used for cmakefilter.

A full parameter list of the parameter_file:

| parameter       | description        |
|-----------------|--------------------|
| href | hs37d5.fa (The absolute path of this file should be given here.) |
| chimp | pt2__cs-hg19.fa (The absolute path of this file should be given here.) |
| hetfa | rawvcf hetfa file (The absolute path of this file should be given here.) |
| fixeddbase | universal filter data base. If this parameter is given, "chimp", "heng75" and "href" will not be used. |
| sampname | sample name for naming the output files (If "mengyao" is given to this parameter, the output mask file name will be "mengyao.mask.fa".) |
| maskname | the mask file generated by cmakefilter |
| gender    | gender of the sample [M/F]   |
| iubfile | data base file of the vcf.hetfa.fa files. If "hetfa" is set, this parameter will not be used. | 
| iubmaskfile | data base file of the mask files generated by cmakefilter. If "maskname" is set, this parameter will not be used. |

Sample output:

```
...
chimp_divergence
no mask:  452284128     0.020393
       0  305600824     0.016474
       1  225977068     0.013937
       2  231055892     0.013612
       3  141395946     0.013435
       4  332362214     0.013229
       5  188300238     0.013079
       6  199918126     0.012964
       7  209035070     0.012888
       8  214724740     0.012822
       9 2354360076     0.012396
  Total: 4855014322     0.013705

Href_divergence
no mask:  652633994     0.001859
       0  349764034     0.001430
       1  243539294     0.000956
       2  247185164     0.000932
       3  150639850     0.000918
       4  352794960     0.000907
       5  199219332     0.000891
       6  211105370     0.000895
       7  220339662     0.000884
       8  225948464     0.000881
       9 2458487018     0.000840
  Total: 5311657142     0.001028
```

Description of the output:

The first block shows the divergence between the input sample and the chimpanzee calculated from the input genotype data. The second block shows the divergence between the input sample and the human reference. The three columns are: mask level, base count, and divergence.

------------------------------- 
<!-- Written by Mengyao Zhao on 2015-10-26 -->
*Updated by Mengyao Zhao on 2015-11-13*
