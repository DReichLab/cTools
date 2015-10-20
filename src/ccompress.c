/*
 * ccompress.c: use razip to compress hetfa files (*.hetfa) and mask files (*.fa)
 * Author: Nick Patterson
 * Revised by: Mengyao Zhao
 * Last revise date: 2015-10-19
 * Contact: mengyao_zhao@hms.harvard.edu 
 */

#include <unistd.h>
#include <nicksam.h>
#include <getpars.h>
#include <ctype.h>
#include "bam.h"
#include "faidx.h"

#define MAXSTR 512
#define MAXFF 20 

int reglen ;
char *regname = NULL ; 
char *regstring = NULL ;
faidx_t *fai;
int chimpmode = NO ;
//char *iubfile = "/n/data1/hms/genetics/reich/1000Genomes/cteam_remap/B-cteam_lite/v0.2/C-FullyPublic__SignedLetterNoDelay__SignedLetterDelay/C.hetfa.dblist" ;
//char *iubmaskfile = "/n/data1/hms/genetics/reich/1000Genomes/cteam_remap/B-cteam_lite/v0.2/C-FullyPublic__SignedLetterNoDelay__SignedLetterDelay/C.mask.dblist" ;
phandle *ph  = NULL ;

void loadfilebase(char *parname) ;
char *myfai_fetch(faidx_t *fai, char *reg, int  *plen) ;

FILE *fff ; 

char *iname = "S_Irula-1" ;  
//char *wkdir = "./" ;	// default writing directory
char *tempout ;

void readcommands(int argc, char **argv) ;
//int setstring(char *iname, unsigned char *ketfa, unsigned char *countfa, int len)  ;
void writefa(FILE *fff, char *regname, char *rrr)  ;

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   ccompress <reference.fa> <hetfa.fa> \n\n");
//	fprintf(stderr, "Options:\n");
//	fprintf(stderr, "\t-f <hetfa file> Please give absolute path with the file name.\n");
//	fprintf(stderr, "\t-r <Href.fa> Please give absolute path with the file name.\n");
	//fprintf(stderr, "\t-i <sample name> [default: S_Irula-1]\n");
//	fprintf(stderr, "\t-w <output dir> [default: ./]\n");
//	fprintf(stderr, "\t-? Show the instruction.\n\n");
	//fprintf(stderr, "Note: If you are not working on orchestra, -i option is not applicable to you, and you need to use -f and -r options together. If you use -f and -r options, -i option will be ignored.\n\n");
	return 1;
}

int setfalist(char **poplist, int npops, char **iublist) {
	int t;
	for (t = 0; t < npops; ++t) iublist[t] = strdup(poplist[t]);
	return npops;
}  

void readfa(char **falist, char **fasta, int *flen, int n) 
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
  fai = fai_load(falist[k]) ;
  fasta[k] = myfai_fetch(fai, regname, &len) ;
  flen[k] = len ;
 }

}

int main(int argc, char *argv[])
{
 char *faname, *p;
 int x ; 
 char ss[1000] ;
 char *falist[2] ;  
 char *famask[2] ;  
 char *iublist[2] ;  
 char *iubmask[2] ;  
char *poplist[2] ;  
 char *fasta[2] ;  
 int flen[2] ;
 char c1, c2 ; 
 int npops =2, len ; 
 char outfaname[200]  ;
 char tmpfaname[200]  ;
 char fainame[200] ; 
 char **reglist ;
	char *dir, *file;
 int nregs, k;

 readcommands(argc, argv);

//	if (f == 1) {
	poplist[0] = argv[1];
	poplist[1] = argv[2];
	p = strrchr(poplist[1], '.');
	x = strlen(poplist[1]);
	iname = (char*)malloc(x*sizeof(char));
	memset (iname, '\0', sizeof(iname));
//fprintf(stderr, "poplist[1]: %s, x:%d\n", poplist[1], x);
	if (!strcmp(p, ".fa")) strncpy(iname, poplist[1], x - 3); 
	else if (!strcmp(p, ".hetfa")) strncpy(iname, poplist[1], x - 6);
	p = strrchr(iname, '/');
	++p;
	file = (char*)malloc(32 * sizeof(char));
	file = strrchr(poplist[1], '/');
	dir = (char*)malloc(x*sizeof(char));
	strncpy(dir, poplist[1], x - strlen(file));
//fprintf(stderr, "dir: %s\n", dir);

	  setfalist(poplist, 2, iublist)  ;
	  setfalist(poplist, 2, iubmask)  ;
	realloc(iubmask[1], (x + 8) * sizeof(char));
	sprintf(iubmask[1], "%s.filter.fa", iname);

  strcpy(fainame, iublist[1]) ;
  strcat(fainame, ".fai") ;
  nregs = numlines(fainame) ;
  ZALLOC(reglist, nregs, char *) ;
  nregs = getss(reglist, fainame) ;
  printf("regs:\n") ;
  printstringsw(reglist, nregs, 5, 10) ;
  

  //sprintf(outfaname, "%s/%s.ccomp.fa", wkdir, iname) ;  
  
  //sprintf(tmpfaname, "%s/%s.ccomptmp.fa", wkdir, iname) ;  
  sprintf(outfaname, "%s/%s.ccomp.fa", dir, p) ;  
  
  sprintf(tmpfaname, "%s/%s.ccomptmp.fa", dir, p) ;  
  openit(tmpfaname, &fff, "w") ;

  fasta[0] = fasta[1] = NULL ;  

  for (k=0; k<nregs; ++k) { 

   regname = reglist[k] ;
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
 sprintf(ss, "htsbox razip %s", tmpfaname) ;  
 system (ss) ;  
 sprintf(ss, "mv %s.rz %s.rz", tmpfaname, outfaname) ;
 system (ss) ;  
 //sprintf(outfaname, "%s/%s.ccompmask.fa", wkdir, iname) ;  
 sprintf(outfaname, "%s/%s.ccompmask.fa", dir, p) ;  
 sprintf(ss, "cp %s %s", iubmask[1], outfaname) ; 
//fprintf(stderr, "ss: %s\n", ss);
 system(ss) ;
 sprintf(ss, "htsbox razip  %s",  outfaname) ; 
 system(ss) ;

//	sprintf(ss, "rm %s/%s.ccomptmp.fa", dir, p);
//	sprintf(ss, "rm %s.ccomptmp.fa", p);
//	system(ss);
  printf("## end of ccompress\n") ;

	free(dir);
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

/*  while ((i = getopt (argc, argv, "w:?")) != -1) {

    switch (i)
      {

      case 'i':
	iname = strdup(optarg) ;
	break;

      case 'w':
	wkdir = strdup(optarg) ;
	break;

      case 'f':
	//	f = 1;
		poplist[1] = strdup(optarg) ;
	break;

      case 'r':	
	//	f = 1;
	poplist[0] = strdup(optarg) ;
	break;

      case '?':
//	default:
//	printf ("Usage: bad params.... %c\n", i) ;
	exit(usage());
	//fatalx("bad params flag:%c\n, i") ;
      }
  }*/
//fprintf(stderr, "argc: %d\n", argc);
	if (argc < 3) exit(usage());

}

