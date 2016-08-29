/*
 * ccompress.c: use gzip to compress hetfa files (*.hetfa) and mask files (*.fa)
 * Author: Nick Patterson
 * Revised by: Mengyao Zhao
 * Last revise date: 2014-11-19
 * Contact: mengyao_zhao@hms.harvard.edu 
 */

#include <unistd.h>
#include <nicksam.h>
#include <globals.h> 
#include <mcmcpars.h>
#include <getpars.h>
#include <ctype.h>
#include "bam.h"
#include "faidx.h"
#include "admutils.h"
#include "mcio.h"  

#define MAXSTR 512

int reglen ;
char *regname = NULL ; 
char *regstring = NULL ;
faidx_t *fai;
int chimpmode = NO ;
//char *iubfile = "/home/np29/cteam/release/hetfaplus.dblist" ;
//char *iubmaskfile = "/home/np29/cteam/release/maskplus.dblist" ;
char *iubfile = "/home/mz128/cteam/dblist/hetfa_postmigration.dblist" ;
char *iubmaskfile = "/home/mz128/cteam/dblist/mask_postmigration.dblist" ;
phandle *ph  = NULL ;

SNP **snpmarkers ;
int numsnps ; 

char *trashdir = "/var/tmp" ;
extern int verbose  ;
int qtmode = NO ;

void loadfilebase(char *parname) ;
char *myfai_fetch(faidx_t *fai, char *reg, int  *plen) ;

FILE *fff ; 

char *faname = NULL ;         
char *snpname = NULL ;         
char *wkdir = "./" ;	// default writing directory
char *tempout ;

void readcommands(int argc, char **argv) ;
int setstring(char *iname, unsigned char *ketfa, unsigned char *countfa, int len)  ;
void writefa(FILE *fff, char *regname, char *rrr)  ;

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   ccompress [options] \n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-i <sample name> [default: S_Irula-1]\n");
	fprintf(stderr, "\t-w <working dir / output dir> [default: ./]\n");
	fprintf(stderr, "\t-? Show the instruction.\n\n");
//	fprintf(stderr, "\t-r <chromosome number> \t The chromosome number can be 1-22, X, Y, MT. [default:	all chromosomes]\n\n");
	return 1;
}


int main(int argc, char *argv[])
{
 int x ; 
 char ss[1000] ;
 char *falist[2] ;  
 char *famask[2] ;  
 char *poplist[2] ;  
 char *iublist[2] ;  
 char *iubmask[2] ;  
 char *fasta[2] ;  
 int flen[2] ;
 char c1, c2 ; 
 int npops =2, len ; 
 char outfaname[200]  ;
 char tmpfaname[200]  ;
 char fainame[200] ; 
 char **reglist ;
 int nregs, k, j, t, pos ;
 int numvind, nignore, numrisks = 1 ;
 int loreg[100], hireg[100] ;
 SNP *cupt ;
 char *cval ;

 readcommands(argc, argv);

  if (faname == NULL) fatalx("no -i flag\n") ;
  if (snpname == NULL) fatalx("no -s flag\n") ;
  
  printf("##grabfa -i %s", faname) ;
  printf(" -s %s\n", snpname) ;
  printnl() ;
  iublist[0] = strdup(faname) ;
  strcpy(fainame, faname) ;
  strcat(fainame, ".fai") ;
  nregs = numlines(fainame) ;
  ZALLOC(reglist, nregs, char *) ;
  nregs = getss(reglist, fainame) ;
  nregs = MIN(nregs, 99) ;
  printf("##regs:\n") ;
  printstringsw(reglist, nregs, 5, 10) ;

  numsnps = 
    getsnps(snpname, &snpmarkers, 0.0,  NULL, &nignore, numrisks) ;
   printf("numsnps:  %d\n", numsnps) ;
  ivclear(loreg, 1000*1000*1000, 100) ;
  ivclear(hireg, 0, 100) ;
  ZALLOC(cval, numsnps, char) ;  // dummy
  cclear(cval, CTAB, numsnps) ;
  for (k=0; k<numsnps; k++) { 
   cupt = snpmarkers[k] ; 
   t = cupt -> chrom ;
   if (t>=100) continue ; 
   loreg[t] = MIN(loreg[t], k) ;
   hireg[t] = MAX(loreg[t], k) ;
  }

  for (t=1; t<=24; ++t) { 
   printf("zzlohi:  %d  %d  %d\n", t, loreg[t], hireg[t]) ;
  }

  fasta[0] = fasta[1] = NULL ;  

  for (k=0; k<nregs; ++k) { 

   t = k + 1 ; 
   regname = reglist[k] ;
   printf("zzregname: %d %s\n", k, regname) ;
   if (loreg[t] > hireg[t]) continue ;
   readfa(iublist, fasta, flen, 1) ; 
   len = flen[0] ;
   printf("zz reg: %s %d\n", regname, len) ;
   for (j=loreg[t]; j <= hireg[t]; ++j) {
     cupt = snpmarkers[j] ; 
     pos = nnint(cupt -> physpos) ;
     if (pos<0) continue ;
     if (pos >= len) continue ;
     cval[j] = fasta[0][pos-1] ; 
   }
  }
  for (j=0; j<numsnps; ++j) { 
   cupt = snpmarkers[j] ;
   printf("snp: %s %d %12.6f %12.0f ", cupt -> ID, cupt -> chrom, cupt -> genpos, cupt -> physpos) ;
   printf(" %c %c ", cupt -> alleles[0], cupt -> alleles[1]) ;
   if (cval[j] == CTAB)  printf ("NA") ;
   else printf ("%c", cval[j]) ;
   printnl() ;
  }
  printf("## end of grabfa\n") ;

  return 0 ;

}

void loadfilebase(char *parname) 
{

  ph = openpars(parname) ;
  dostrsub(ph) ;

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
  substring(&treg, "chr", "") ;
  s = fai_fetch(fai, treg, plen) ; 
  if (*plen > 0) { 
    free(treg) ;
    return s ;
  }
  free(treg) ; 
  return NULL ;
}

void readcommands(int argc, char **argv) 
{
  int i;
  phandle *ph ;
  char str[512]  ;

  while ((i = getopt (argc, argv, "i:f:s:r:w:?")) != -1) {

    switch (i)
      {

      case 'i':
	faname = strdup(optarg) ;
	break;

      case 'f':
	faname = strdup(optarg) ;
	break;

      case 'w':
	wkdir = strdup(optarg) ;
	break;

      case 'r':	// Currently this doesn't work.
	regname = strdup(optarg) ;
	break;

      case 's':	
	snpname = strdup(optarg) ;
	break;

      case '?':
	exit(usage());
	break;

	default:
//	printf ("Usage: bad params.... %c\n", i) ;
	fatalx("bad params flag:%c\n, i") ;
      }
  }

}

int getfalist(char **poplist, int npops, char *dbfile, char **iublist) 
{
 char line[MAXSTR+1] ;
 char *spt[MAXFF], *sx ;
 char c ;
 int nsplit ;
 int  t, k, s, nx = 0  ;
 FILE *fff ;
 char *scolon ; 
  
  if (dbfile == NULL) return 0 ;
  openit(dbfile, &fff, "r") ;

  while (fgets(line, MAXSTR, fff) != NULL)  { 	// Read the dblist file.
   nsplit = splitup(line, spt, MAXFF) ; // Store the line into the array spt, devided by \s.
   if (nsplit<1) continue ;
   sx = spt[0] ;
   if (sx[0] == '#') { 
    freeup(spt, nsplit) ;
    continue ;
   }

//	npops is the length of the array poplist. 
//	if spt[0] is the t th element of poplist, return t; spt[0] is not in poplist, return -1.
   t = indxstring(poplist, npops, spt[0]) ;	
   if (t<0) { 
    freeup(spt, nsplit) ; 
    continue ;
   }
    sx = spt[2] ;
    iublist[t] = strdup(sx) ;
    ++nx ;
    freeup(spt, nsplit) ;
  }

   fclose(fff) ;
   return nx ;

}
int readfa(char **falist, char **fasta, int *flen, int n) 
{

 faidx_t *fai = NULL ;
 int k, len ;
 
 ivzero(flen, n) ;
 for (k=0; k<n; ++k) {

  if (falist[k] == NULL) { 
   fasta[k] = NULL ;
   continue ;
  }
  if (fasta[k] != NULL) freestring(&fasta[k]) ;
//fprintf(stderr, "falist[%d]: %s\n", k, falist[k]);
  fai = fai_load(falist[k]) ;
  fasta[k] = myfai_fetch(fai, regname, &len) ;
  flen[k] = len ;
 }

}

