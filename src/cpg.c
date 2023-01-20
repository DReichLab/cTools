/**
* cpulldown.c:	get the genotypes of the given individuls at the given SNP loci from a set of bams
* Author: Nick Patterson
* Revised by: Mengyao Zhao
* Last revise date: 2015-04-27
* Contact: nickp@broadinstitute.org    
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <zlib.h>

#include <nicklib.h>
#include <globals.h>
#include <mcmcpars.h>
#include <getpars.h>
#include <nicksam.h>

#include "faidx.h"
#include "admutils.h"
#include "mcio.h"  
#include "ctools.h"  

#define MAXFL  50   
#define MAXSTR  512

extern int checksizemode ;
char *omode = "eigenstrat" ;
extern int packmode ;            //!< flag - input {is not,is} in packed mode
extern char *packgenos ;         //!< packed genotype data (packit.h)
extern char *packepath ;
extern long packlen;             //!< allocated size of packgenos data space
extern long rlen;                //!< number of bytes in packgenos space that each SNP's data occupies

char *trashdir = "/var/tmp" ;
char *chimpname = "/home/np29/broaddatax/tables//pt2__cs-hg19.fa" ;
faidx_t *chimpfai ;

extern int verbose  ;
int debug = NO ;
int qtmode = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 

char *polarid = NULL ;
int polarindex = -1 ;

char  *genotypename = NULL ;
char  *snpname = NULL ;

char *outname = NULL ;

int checkmode = NO ;


int minchrom = 1 ;
int maxchrom = 23 ;
int xchrom = -1 ;
int db = 1;	// Use .dblist

char *outputname = NULL ;
FILE *ofile ;

double fakespacing = 0.0 ;
char  unknowngender = 'U' ;

char *regname ;
char **iublist ;
char **iubmask ;

char *regstring = NULL ;
char *mask ;

int reglen ; 

void readcommands(int argc, char **argv) ;
char *myfai_fetch(faidx_t *fai, char *reg, int  *plen) ;
char xrefbase(int pos)  ;
int qcpg(int pos, char *alleles)  ;

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   cpg -s snpfile [-o outfile] [-v] [-V]\n" ); 
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-s	snp file (.snp Reich lab format (or plink .map file))\n") ; 
	fprintf(stderr, "\t-o	output file (coding: -1 => No call, 0  => Not CpG,  1 => CpG) \n" );  
	fprintf(stderr, "\t-V	Print more information while the program is running.\n");
	fprintf(stderr, "\t-v	Show version information.\n");
	return 1;
}

int main(int argc, char **argv)
{

  int i, j, k, kk, g, t, tt, ccheck ; 
  int pos, isnp, gval ;
  SNP *cupt ;
  Indiv *indx ;

  int numvind, nignore, numrisks = 1 ;
  int chrom, lastchrom = -1 ;  
  int firstsnp, rlen, mlen  ;
  char *iubarr, *maskarr  ;
  char **dbchar, *dbx ;  
  int *cpgvals ; 

  char ss[100], cbases[2], iub, cm, c1, c2, cval ;

  ofile = stdout; 
  packmode = YES ;
  printf("cpg: version %s\n", version) ; 
  chimpfai = fai_load(chimpname);
  readcommands(argc, argv) ;

  settersemode(YES) ;
  if (outputname != NULL) openit(outputname, &ofile, "w") ;
  if (snpname == NULL) exit(usage()) ;

  numsnps = 
    getsnps(snpname, &snpmarkers, fakespacing,  NULL, &nignore, numrisks) ;

  printf("numsnps: %d\n", numsnps) ;

  ZALLOC(cpgvals, numsnps, int) ; 
  ivclear(cpgvals, -1, numsnps) ;  // -1 unset :: 0 = NoCpG ::  1 > CpG
  lastchrom = - 1 ;

    fprintf(ofile, "%20s ", "ID") ;       
    fprintf(ofile, "%6s %12s ", "Chrom", "Pos") ; 
    fprintf(ofile, "  %4s", "Code") ;
    fprintf(ofile,  "  -1 => No Call;  0 => Not CpG;  1 => CpG" ) ;
    fprintf(ofile, "\n") ; 
  for (isnp = 0 ; isnp < numsnps; ++isnp) { 
   cupt = snpmarkers[isnp]  ; 
   chrom = cupt -> chrom ;
   ccheck = YES ;
   if (chrom < minchrom) ccheck = NO  ;
   if (chrom > maxchrom) ccheck = NO  ;
   if ((xchrom>0) && (chrom != xchrom)) ccheck = NO  ;
   if (ccheck == NO)   { 
    cupt -> ignore = YES ;
    continue ;
   }

  if (lastchrom != chrom) {
   lastchrom = chrom ;
   if (chrom < minchrom) continue ;
   if (chrom > maxchrom) continue ;
   if ((xchrom>0) && (chrom != xchrom)) continue ;
    sprintf(ss, "%d", chrom) ;

   if (chrom == 23) strcpy(ss, "X") ;
   if (chrom == 24) strcpy(ss, "Y") ;
   if (chrom == 25) strcpy(ss, "MT") ;

    freestring(&regname) ;
    regname = strdup(ss) ;

    freestring(&regstring) ; 
    regstring = myfai_fetch(chimpfai, regname, &reglen) ;
    printf("## reg: %s len: %d\n", regname, reglen) ;
   }

    pos = nnint(cupt -> physpos) ;  
    t = qcpg(pos, cupt -> alleles) ;
    fprintf(ofile, "%20s ", cupt -> ID) ; 
    fprintf(ofile, "%6d %12d ", chrom, pos) ; 
    fprintf(ofile, "  %4d", t) ;
    fprintf(ofile, "\n") ; 

 } 

  if (outputname != NULL) fclose(ofile) ;

  printf("##end of run\n") ;
  return 0 ;
}

char xrefbase(int pos) 
{

 if (pos <= 0) return '?' ;
 if (pos >  reglen) return '?' ;
 return toupper(regstring[pos-1]) ;

}

void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "s:o:Vv")) != -1) {

    switch (i)
      {

      case 's':
	snpname = strdup(optarg) ;
	break;

      case 'o':
	outputname = strdup(optarg) ;
	break;


      case 'v':
	printf("version: %s\n", version) ; 
	break; 

      case 'V':
	verbose = YES ;
	break; 

	default:
	exit(usage());

      }
  }
}

char *myfai_fetch(faidx_t *fai, char *reg, int  *plen)
{
  char *treg, *s ;
  treg = strdup(reg) ;

  if (fai==NULL) fatalx("(my_fai_fetch): fai NULL\n") ;

  s = fai_fetch(fai, treg, plen) ;
  if (*plen > 0) {
    free(treg) ;
    return s ;
  }
  free(treg) ;
  return NULL ;
}

int qcpg(int pos, char *alleles) 
{ 
  int  t ;
  char cref, c ;
  static int ncall = 0 ;
  int pos2 ;

  ++ncall ;

   cref = xrefbase(pos) ;
   if ((cref != alleles[0]) && (cref != alleles[1])) return -1 ;   
   
   t = 0 ; 

   if ((alleles[0] == 'C') && (alleles[1] == 'T')) t=1 ;
   if ((alleles[0] == 'T') && (alleles[1] == 'C')) t=1 ;
   if ((alleles[0] == 'G') && (alleles[1] == 'A')) t=2 ;
   if ((alleles[0] == 'A') && (alleles[1] == 'G')) t=2 ;

   if (t==0) return 0 ;
   pos2 = pos + 1 ;
   if (t==2) pos2 = pos-1 ;
   c = xrefbase(pos2) ;
   if (base2num(c) < 0) return -1 ;
   if ((t==1) && (c=='G')) return 1 ;
   if ((t==2) && (c=='C')) return 1 ;
   return 0 ;
}


