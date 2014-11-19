/*
 * ccompress.c: use gzip to compress hetfa files (*.hetfa) and mask files (*.fa)
 * Author: Nick Patterson
 * Revised by: Mengyao Zhao
 * Last revise date: 2014-11-19
 * Contact: mengyao_zhao@hms.harvard.edu 
 */

#include <unistd.h>
#include <nicksam.h>
#include <getpars.h>
#include "bam.h"
#include "faidx.h"

#define MAXSTR 512
#define MAXFF 20 

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

void loadfilebase(char *parname) ;
char *myfai_fetch(faidx_t *fai, char *reg, int  *plen) ;

FILE *fff ; 

char *iname = "S_Irula-1" ;  
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
//	fprintf(stderr, "\t-r <chromosome number> \t The chromosome number can be 1-22, X, Y, MT. [default:	all chromosomes]\n\n");
	return 1;
}


int main(int argc, char *argv[])
{
 char *faname ;
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
 int nregs, k ;

// regname = strdup("22") ;
 readcommands(argc, argv);

  poplist[0] = strdup("Href") ;

  poplist[1] = strdup(iname) ;
  getfalist(poplist, 2, iubfile, iublist)  ;
  getfalist(poplist, 2, iubmaskfile, iubmask)  ;
  strcpy(fainame, iublist[1]) ;
  strcat(fainame, ".fai") ;
  nregs = numlines(fainame) ;
  ZALLOC(reglist, nregs, char *) ;
  nregs = getss(reglist, fainame) ;
  printf("regs:\n") ;
  printstringsw(reglist, nregs, 5, 10) ;
  

  sprintf(outfaname, "%s/%s.ccomp.fa", wkdir, iname) ;  
  
  sprintf(tmpfaname, "%s/%s.ccomptmp.fa", wkdir, iname) ;  
  openit(tmpfaname, &fff, "w") ;

  fasta[0] = fasta[1] = NULL ;  

  for (k=0; k<nregs; ++k) { 

   regname = reglist[k] ;
//fprintf(stderr, "iublist[0]: %s\niublist[1]: %s\n", iublist[0], iublist[1]);
//fprintf(stderr, "iublist: %s\n", iublist);
   readfa(iublist, fasta, flen, 2) ; 

   len = MIN(flen[0], flen[1]) ;
   for (x=0; x<len; ++x) { 
    c1 = fasta[0][x] ;
    c2 = fasta[1][x] ;

    if (toupper(c1) == c2) { 
     fasta[1][x] = 'Q' ;
    }
   }

  writefa(fff, regname, fasta[1]) ;
  printf("chromosome: %12s done\n", regname) ;
 }
 fclose(fff) ;
 sprintf(ss, "gzip %s", tmpfaname) ;  
 system (ss) ;  
 sprintf(ss, "mv %s.gz %s.gz", tmpfaname, outfaname) ;
 system (ss) ;  
 sprintf(outfaname, "%s/%s.ccompmask.fa", wkdir, iname) ;  
 sprintf(ss, "cp %s %s", iubmask[1], outfaname) ; 
 system(ss) ;
 sprintf(ss, "gzip  %s",  outfaname) ; 
 system(ss) ;

  printf("## end of ccompress\n") ;

	free(iname);
  return 0 ;

}

void writefa(FILE *fff, char *regname, char *rrr) 
{
  int len, tlen, x ; 
  char sss[51], *sx, *pp ;

  fprintf(fff, ">%s\n", regname) ;
  tlen = len = strlen(rrr) ;
  sss[50] = CNULL ;
  sx = rrr ; 
  for (;;) { 
   if (tlen <= 0) break ;
   x = MIN(tlen, 50) ;
   strncpy(sss, sx, x) ;
   sss[x] = CNULL ;
   fprintf(fff, "%s\n", sss) ;
   sx += 50 ;
   tlen -= strlen(sss) ; 
  }
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
  int n, kode ;
  int pops[2] ;

  while ((i = getopt (argc, argv, "i:r:w:?")) != -1) {

    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case 'w':
	wkdir = strdup(optarg) ;
	break;

      case 'r':	// Currently this doesn't work.
	regname = strdup(optarg) ;
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

