/*
 * cpoly.c: This program is used to extract heterozygote SNPs from multiple samples
 * Author: Nick Patterson
 * Revised by: Mengyao Zhao
 * Last revise date: 2015-04-30
 * Contact: nickp@broadinstitute.org     
 */

// memory leak. wwill fix soon

#include <sys/wait.h>
#include <stdlib.h>
#include <stdio.h>
#include <nicksam.h>
#include <getpars.h>  
#include <mcmcpars.h>
#include <time.h>  
#include "bam.h"
#include "faidx.h"
#include "globals.h" 
#include "popsubs.h"
#include "mcio.h"
#include "ctools.h"

typedef struct { 
 int *val ; 
 int *maxval ; 
 int ascnum ;
 int np ;
} ASC ;

//char *iubfile = NULL ;
//char *iubmaskfile = NULL ;

char *regname = NULL ; 
//char *parflist = "/home/np29/biology/neander/nickdir/xwdir/may12src/parfxlm" ;

char *parname = NULL ;
int  pagesize = -1 ;  // page size for getiub
int minfilterval = 1 ;
int readfailOK = YES ;

int minderiv = -1 ; 
int maxderiv = 99999 ; 

int minchrom = 1 ;
int xchrom = -1 ;
int maxchrom = 25 ;
char *minch = NULL;
char *maxch = NULL;
int debug = NO ;

char *trashdir = "/var/tmp" ;
int qtmode = NO ;

int doxchrom = YES ; 
int doychrom = YES ; 
int domt = NO ; 


char *polarid = NULL ;
int polarindex = -1 ;

int allowmissing = YES  ;
int allowhets = YES  ;

char *indivname = NULL ;
char *indoutfilename = NULL ;
char *snpoutfilename = NULL ;
char  *genooutfilename = NULL ;

int lopos, hipos  ;
char **fasta ;  // in core bases
char *ascstring = NULL ;
char *noascstring = NULL ;

char *poplistname = NULL ; 
char **poplist ; 
int *hasmask ;
int npops = 0 ;


// bugfix bug when polarize off.   Last pop het didn't work
void readcommands(int argc, char **argv) ;
int getfalist(char **poplist, int npops, char *dbfile, char **iublist) ; 
void clearfainfo(FATYPE *fapt, int mode) ;
//char *myfai_fetch(faidx_t *fai, char *reg, int  *plen) ;
int loadfa(char **poplist, int npops, FATYPE ***pfainfo, char *reg, int lopos, int hipos)  ;
int minfalen(FATYPE **fainfo, int n)  ;


int istransition(char iubc)  ;
int abx(int a, int b)  ; 
int abxok(int abx, int abxmode) ;   
void copyfq(char **names, char **data, int *len, int a, int b) ;
int getiub(char *cc, char *ccmask, FATYPE **fainfo, char *reg, int pos) ; 

int checkpoly(char *cc, char *ccmask, char *pc1, char *pc2)   ;
void printfapt(FATYPE *fapt);
void prints(FILE *fff, int pos, char c1, char c2)  ;
void printgg(FILE *ggg, char *cc, char *ccmask, char c1, char c2, int n) ;

int abxmode = 0 ;
int checkmode = NO ;

ASC **asctable ;
ASC **noasctable ;
int nasc, nonasc ;

static int usage() 
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   cpoly -p <parameter file> [options] \n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "Options: dbhetfa, dbmask now obligatory\n") ;
	fprintf(stderr, "\t-V	Print more information while the program is running.\n");
	fprintf(stderr, "\t-v	Show version information.\n");
	fprintf(stderr, "\t-? 	Show the instruction. (For detailed instruction, please see the document here: https://github.com/mengyao/cTools)\n\n");
	fprintf(stderr, "\t-c 	Checkmode only; don't run\n");
	return 1;
}

int main(int argc, char **argv)
{

 char *spt[MAXFF], *zspt[MAXFF] ;
 int  k, t, t2, q, g, isbad ;
 ASC *ascpt ;
 char *reg ;
 FATYPE **fainfo, *fapt ; 
 int xnpops ;
 int pos, chrom ;
 char *cc, *ccmask ;
 char c1, c2 ; 
 FILE  *fff = NULL, *ggg = NULL ;;
 char ss[1024] ;
 int lo, hi ;

 int numout = 0 ;
 int abxkode ;
 int nmono ;
 int npoly ;
 double ymem ; 

	clock_t start, end;
	float cpu_time;	

	start = clock(); 
 
 hipos = 1000*1000*1000 ;
 lopos = 0 ;
 printf("cpoly: version %s\n", version) ; 

 readcommands(argc, argv) ;

  cputime(0) ;
  calcmem(0) ;


	if (minch != NULL) {
		if (minch[0] == 'X') minchrom = 23;
		else if (minch[0] == 'Y') minchrom = 24;
		else if (!strcmp(minch, "MT")) minchrom = 25; 
		else minchrom = atoi(minch);
        }

         if (maxch != NULL) {
		if (maxch[0] == 'X') maxchrom = 23;
		else if (maxch[0] == 'Y') maxchrom = 24;
		else if (!strcmp(maxch, "MT")) {
                 maxchrom = 25; 
                 domt = YES ; 
                }
		else maxchrom = atoi(maxch);
	}

 if (indivname==NULL) {
   usage() ;
   printf("indivname: omitted\n") ;
   return 1 ;
 }

 if (snpoutfilename != NULL) openit(snpoutfilename, &fff, "w") ;
 else fff = stdout ; 
 if (genooutfilename != NULL) openit(genooutfilename, &ggg, "w") ;

 if (regname != NULL) { 
  if (regname[0] == 'X') xchrom = 23 ; 
  	else if (regname[0] == 'Y') xchrom = 24 ;
	else if (!strcmp(regname, "MT")) xchrom = 25; 
  else xchrom = atoi(regname) ;
 }
 if (xchrom==25) domt = YES ;

 //  fprintf(stderr, "call numlines in [main]\n");
 
  npops = numlines(indivname) ; 
  ZALLOC(poplist, npops, char *) ;
  npops = getss(poplist, indivname) ;
  
 if (npops > 1000) fatalx("too many samples\n") ;

 printf("samplist:\n") ;
 printstrings(poplist, npops) ;

  if (indoutfilename != NULL) {  
    sprintf(ss, "cp %s %s", indivname, indoutfilename) ;
    system (ss) ; 
    printf("%s written\n", indoutfilename) ;
  }

  ZALLOC(hasmask, npops, int) ;
  
  if (polarid != NULL) {
   polarindex = indxstring(poplist, npops, polarid) ;
   if (polarindex<0) fatalx("polarid %s not found\n", polarid) ;
  }

 if (pagesize < 0) { 
  pagesize = (1000*1000*1000)/(2*npops) ;
 }
 xnpops = npops+1  ;

	if (xchrom > 0)  { 
	  chrom = xchrom  ;
	  sprintf(ss, "%d", chrom) ;
	  if (chrom == 23) strcpy(ss, "X") ;
		else if (chrom == 24) strcpy(ss, "Y") ;
		else if (chrom == 25) strcpy(ss, "MT") ;

	 }
	 else sprintf(ss, "%d", 22) ;
	 regname = strdup(ss) ;
	 reg = regname ;

 lo = lopos ; 
 hi = MIN(hipos, lo+pagesize) ;
 loadfa(poplist, npops, &fainfo, reg, lo, hi)  ;
//fprintf(stderr, "main: fapt->fai: %p\n", fainfo[0]->fai);
 printf("npops: %d\n", npops) ;
  for (k=0;  k< npops; ++k) { 
   fapt = fainfo[k] ;
   printfapt(fapt) ;
  }

 ZALLOC(cc, xnpops, char) ; 
 ZALLOC(ccmask, xnpops, char) ; 
 cc[npops] = CNULL ;
 ccmask[npops] = CNULL ; // don't test chimp

 printf("minchrom: %d  maxchrom %d\n", minchrom, maxchrom) ;
 if (checkmode) { 
  printf("checkmode set.  End of Run\n") ;
  return 0 ;  
 }


 nmono = npoly = 0 ;
 for (chrom = minchrom; chrom <= maxchrom; ++chrom) {
  if ((xchrom > 0) && (xchrom != chrom)) continue ;
  if ((chrom == 23) && (doxchrom == NO)) continue ;  
  if ((chrom == 24) && (doychrom == NO)) continue ;  
  if ((chrom == 25) && (domt == NO)) continue ;  

  sprintf(ss, "%d", chrom) ;
  if (chrom == 23) strcpy(ss, "X") ;
  if (chrom == 24) strcpy(ss, "Y") ;
  if (chrom == 25) strcpy(ss, "MT") ;

  freestring(&regname) ;
  regname = strdup(ss) ;
  reg = regname ;
	
  for (pos = lopos ; pos <= hipos; ++pos) { 
    t = getiub(cc, ccmask, fainfo, ss, pos)  ;  
    if (debug) printf("zziub %d\n%s\n%s\n\n", t, cc, ccmask) ;
   if (t==-5) break ;
   if (t<0) continue ;
    
   t = checkpoly(cc, ccmask, &c1, &c2) ; 
   if (t==NO) continue ;
    if (abxmode != 0) {
     abxkode = abx(base2num(c1), base2num(c2)) ;
     if (abxkode < 0) continue ;
     t = abxok(abxkode, abxmode) ;
     if (t==NO) continue ;
    }
   if (c1==c2) { 
// hit
    ++nmono ;
    continue ;
   }
   ++npoly ;

  isbad = NO  ;
  if ((polarindex >= 0) && (allowmissing==NO))  {
   t = 0 ;
   for (k=0; k<npops; ++k) {
    g = gvalm(toupper(cc[k]), ccmask[k], c1, c2, minfilterval) ; // symbol to write (0 1 2 9)
    if (g==9) isbad = YES ;
    t += (2-g) ;
   }
   if (t<minderiv) isbad = YES ;
   if (t>maxderiv) isbad = YES ;
  }
  if (isbad) continue ;


   prints(fff, pos, c1, c2) ;
   printgg(ggg, cc, ccmask, c1, c2, npops) ;

/**
    printf("zzgg: %12d %c %c " , pos, c1, c2)   ;
    printf("%s %s\n", cc, ccmask) ; 
    printnl() ; 
*/
   
 }}

 if (snpoutfilename != NULL) fclose(fff) ;
 if (genooutfilename != NULL) fclose(ggg) ;
 
 printf("## monomorphs: %d  polymorphs %d\n", nmono, npoly) ;

  ymem = calcmem(1)/1.0e6 ;
  printf("##end of cpoly: %12.3f seconds cpu %12.3f Mbytes in use\n", cputime(1), ymem) ;
	
 return 0 ;
}


void printgg(FILE *ggg, char *cc, char *ccmask, char c1, char c2, int n) 
// known not to be triiallelic
{
   int k, g ; 

   if (ggg==NULL) return ;
   for (k=0; k<n; ++k) { 
    g = gvalm(toupper(cc[k]), ccmask[k], c1, c2, minfilterval) ; // symbol to write (0 1 2 9)
    fprintf(ggg,"%d", g) ;
   }
   fprintf(ggg, "\n" ) ;
}

void prints(FILE *fff, int pos, char c1, char c2) 
{
  char sss[100]  ;
  sprintf(sss, "X:%s_%d ", regname, pos) ; 
  fprintf(fff, "%15s ", sss) ;
  fprintf(fff, "%3s 0 %12d ", regname, pos) ;
  fprintf(fff, "%c %c\n", c1, c2) ;

  return ;

}

int checktriallelic(char *pc1, char *pc2, char x1, char x2) 
{
  int aa[4], t ;
  char c1, c2 ; 
  c1 = *pc1; c2 = *pc2 ;
  ivzero(aa, 4) ;

  t = base2num(c1) ; if (t>=0) aa[t] = 1 ;
  t = base2num(c2) ; if (t>=0) aa[t] = 1 ;
  t = base2num(x1) ; if (t>=0) aa[t] = 1 ;
  t = base2num(x2) ; if (t>=0) aa[t] = 1 ;

  t = intsum(aa,4) ;
  if (t>2) return NO ;

  if (base2num(c1) < 0) *pc1 = x1 ; 
  if (t==1) return YES ;
  if (base2num(c2) < 0) {
    if (c1==x1) c2=x2 ;
    if (c1==x2) c2=x1 ;
  }
  *pc2 = c2 ;
  return YES ;
}

int checkpoly(char *cc, char *ccmask, char *pc1, char *pc2)  
{
 int valid[20] ;
 int cnum, t, x, xpol ;
 char cm, cx, cchimp, c1, c2, cpol = '?' ; 
 int isdiff = NO, kasc, k, x1, x2 ;
 char cb[2] ;
 int aa[4] ;

 *pc1 = *pc2 = c1 = c2 = '?' ;
 
  ivzero(aa, 4) ;
  xpol = -1 ;
  if (polarindex>=0) { 
    c1 = cc[polarindex] ;
    cpol = c1 = toupper(c1) ;
    xpol = base2num(cpol) ;
    if (xpol<0) return NO ;
    aa[xpol] = 1 ;
    x1 = xpol ; 
  }
  for (k=0; k<npops; ++k) { 
   cm = ccmask[k] ; 
   cx  = toupper(cc[k]) ;
   if (fvalid(cm, minfilterval) == NO)  cx = '?' ;
   t = isiub2(cx) ; 
   if ((t==NO) && (allowmissing==NO)) return NO ;
   if (t==NO) continue  ; 
   t = iubcbases(cb, cx) ; 
   if ((t==2) & (allowhets==NO))    return NO ;
   x = base2num(cb[0]) ; aa[x] = 1 ;
   x = base2num(cb[1]) ; aa[x] = 1 ;
  }
  t = intsum(aa, 4) ;
  if (t > 2) return NO ;
  if (t ==0 ) return NO ;
  if (xpol<0) { 
   x1 = findfirst(aa, 4, 1) ; 
  }
  aa[x1] = 0 ;
  if (t==1) x2 = x1 ; 
  else {
   x2 = findfirst(aa, 4, 1) ;
  }
   *pc1 = num2base(x1) ;
   *pc2 = num2base(x2) ;
   return YES ;
}
int loadfa(char **poplist, int npops, FATYPE ***pfainfo, char *reg, int lopos, int hipos) 
{
 static int k, numfalist, t;
 static char **falist, **famasklist ;
 static FATYPE **fainfo, *fapt ;
 int *falen ;
 int lo, hi, reflen, len, len_r, len_s, i, tt, l, j, isok ;
 unsigned char c ; 

 static int ncall = 0 ;
	char* region, *ref, *refname = (char*)malloc(256*sizeof(char));
	faidx_t *fai_ref;	// Use this to open the reference sequence. 
  
	++ncall ;

  fflush(stdout) ;

	lo = MAX(1, lopos) ;
	region = (char*)malloc((23+strlen(reg))*sizeof(char));
	sprintf(region, "%s:%d-%d", reg, lo, hipos);

     
	getdbname(iubfile, "Href", &refname);

  if (ncall==1) {
   ZALLOC(falist, npops, char *) ;
   ZALLOC(famasklist, npops, char *) ;

	{
	   numfalist = getfalist(poplist, npops, iubfile, falist) ;	// set falist with the absolute path of hetfa files in .dblist file
	   t = getfalist(poplist, npops, iubmaskfile, famasklist) ; 
	}

   if (numfalist != npops) {
      for (k=0; k<npops; ++k) { 
        if (falist[k] == NULL) printf("no fasta file for: %s\n", poplist[k]) ;
      }
	  fatalx("Do not find the data files. Please check dbhetfa and dbmask in your parameter file.\n") ;
   }
   for (k=0; k<npops; ++k) {
    hasmask[k] = YES ;
    if (famasklist[k] == NULL) {  
     hasmask[k] = NO ;
     continue ;
    }
    t = strcmp(famasklist[k], "NULL") ; 
    if (t==0) {
     hasmask[k] = NO ;
     continue ;
    }
   }

   ZALLOC(fainfo, npops, FATYPE *) ;
   ZALLOC(fasta, npops, char *) ;
   for (k=0; k<npops; ++k) {
    ZALLOC(fainfo[k], 1, FATYPE) ;
    fapt = fainfo[k] ; 
    clearfainfo(fapt, 1) ;

    fapt -> faname = strdup(falist[k]) ;
    fapt -> alias = strdup(poplist[k]) ;
    fapt -> famask = strdup("NULL")   ;
    if (hasmask[k]) fapt -> famask = strdup(famasklist[k]) ;

    verbose = YES;
    if (verbose) {
     printf("%s\n", fapt -> alias) ;
     printf("%s\n", fapt -> faname) ;
     printf("%s\n", fapt -> famask) ;
     printf("loading: %s\n", fapt -> faname) ;
     printnl() ;
   }
   verbose = NO ;
    fapt -> fai = fai_load(fapt -> faname) ;
    fapt -> popnum = k ;
    if (hasmask[k]) fapt -> faimask = fai_load(fapt -> famask) ;
   }
  }

  if (ncall > 1) {
    for (k=0; k<npops; ++k) {
     fapt = fainfo[k] ;
     freestring(&fapt -> rstring) ;
     freestring(&fapt -> mstring) ;
     clearfainfo(fapt, 0) ;
    }
  }

  if (pfainfo != NULL) *pfainfo  = fainfo ;

	fai_ref = fai_load(refname);
	ref = fai_fetch(fai_ref, region, &len_r);
	if (len_r==0) fatalx("bad fetch %s %s\n", refname, region) ; 	// fetch fai

  for (k=0; k<numfalist ; ++k) {
	FILE *fp;
	int byte[2], rz = 0;

     fapt = fainfo[k] ;

        openit(fapt -> faname, &fp, "r") ;
	for (i = 0; i < 2; ++i) byte[i] = getc(fp);
	if (byte[0] == 0x1f && byte[1] == 0x8b) rz = 1;
	fclose(fp);
	
        fapt -> rlen = fai_getlen(fapt->fai, reg) ;
        if (lo>=fapt ->rlen) { 
         reflen = fai_getlen(fai_ref, reg) ; 
         printf("trying to read after chromosome end lo: %d len: %d reflen: %d\n", lo, fapt -> rlen, reflen) ; 
        }
	fapt->rstring = fai_fetch(fapt->fai, region, &len_s);
	if (len_s==0) {
         printf("fetch fails\n") ;
         printfapt(fapt) ; 
         if (readfailOK==NO) { 
          fatalx("bad fetch %s %s\n", fapt->alias, region);	// fetch fai
         }
         else { 
          fapt -> rstring = NULL ; 
         }
        }

	len = len_r < len_s ? len_r : len_s;
//	len = len_s;
	if (rz == 1)	// razipped 
		for (i = 0; i < len; ++i) 
			if (fapt->rstring[i] == 'Q') fapt->rstring[i] = ref[i];
	

      fapt -> regname = strdup(reg) ;
      fapt -> len = len ;
      fapt -> lopos = lo ;
      fapt -> hipos = lo + len - 1 ;
      if (verbose) printf("zzlh %d %d %d %d\n", lopos, hipos, fapt -> lopos, fapt -> hipos) ;

      if (fapt -> faimask == NULL) continue ;
                //sprintf(region, "%s:%d-%d", reg, fapt->lopos, fapt->hipos);
		len = 0;
		fapt->mstring = fai_fetch(fapt->faimask, region, &len);
	  if (len==0) { 
           fapt -> mstring = NULL ; 
           printf("*** warning: bad fetch (mask)  %s %s\n", fapt -> alias, region) ;
          }
      fapt -> mlen = len ;
  }

	free(ref);
	fai_destroy(fai_ref);

  return npops ;

}

char getfacc(FATYPE *fapt, int pos, int xmode) 
// xmode 1:  genotype    xmode 2: mask
{
  int t ;
  t = pos - fapt -> lopos ;

  if (t<0) return '-' ;
  if (xmode == 1) {
   if (t>=fapt -> len) return '-' ;
   return toupper(fapt -> rstring[t]) ;
  }
  if (xmode == 2) {
   if (t>=fapt -> mlen) return '-' ;
   return toupper(fapt -> mstring[t]) ;
  }
  fatalx("badbug\n") ;
}

int getiub(char *cc, char *ccmask, FATYPE **fainfo, char *reg, int pos) 
/** return 
 -2 (pos out of range) 
 -1 polarization fails
 -5 end of chromosome 
  else number of valids 
 
*/
{
  int k, t ;  
  FATYPE *fapt ;
  char *lastreg  ; 
  static int  lastlo = 1000*1000*1000, lasthi=-1 ;
  int newpage = NO, newreg = NO ;
  char regbuff[128] ;
  static long ncnt = 0 ;
  static long ncall = 0 ;

  ++ncall ;
  cclear(cc, '?', npops) ;
  cclear(ccmask, '?', npops) ;
  fapt = fainfo[0] ; 
  lastreg = fapt -> regname ;
  strcpy(regbuff, reg) ;
  if (lastreg == NULL) lastreg = strdup("??") ;

  if (strcmp(lastreg, regbuff) != 0) {
   newpage = YES ;
   newreg =  YES; 
  }
  else { 
   if (pos >= fapt -> rlen) return -5 ;
  }

  lastlo = 1000*1000*1000 ; 
  lasthi = 1 ; 
  for (k=0; k<npops; ++k) { 
   fapt = fainfo[k] ; 
   lastlo = MIN(lastlo, fapt -> lopos) ;
   lasthi = MAX(lasthi, fapt -> hipos) ;
  }

  if (pos < lopos) return -3 ;
  if (pos > hipos) return -3 ;

  if (pos < lastlo) newpage = YES ;
  if (pos > lasthi) newpage = YES ;
  if (ncall == 1) newreg = YES ;

  if (newreg == YES) { 
// printf("zznewrrr %s :: %s %d %d %d\n",lastreg, regbuff,  pos, lastlo, lasthi) ;  fflush(stdout) ;
   fflush(stdout) ;
   freestring(&regname) ;

   regname = strdup(regbuff)  ;  
   lastreg = strdup(regbuff)  ;  

// set falen
   
  }

  if (newpage == YES) { 
   fprintf(stderr, "pos-getiub: %d\n", pos);
   fflush(stdout) ;
   lastlo = pos ;
   lasthi = pos + pagesize ;
   lastlo = MAX(lastlo, lopos) ;
   lasthi = MAX(lasthi, hipos) ;
   lasthi = MIN(lasthi, lastlo+pagesize) ;
   printf("calling loodfa %s %d %d \n", regname, lastlo, lasthi) ;
   fflush(stdout) ;
//fprintf(stderr, "getiub: fapt->fai: %p\n", fainfo[0]->fai);
   loadfa(poplist, npops, &fainfo, regname, lastlo, lasthi)  ;
   printf("newpage: %d %p %d %d\n", pos, topheap(), lastlo, lasthi) ;

   fflush(stdout) ;
  }

  for (k=0; k<npops; ++k) { 
   fapt = fainfo[k] ;
// if (pos<fapt -> lopos) return -2 ;
// if (pos>fapt -> hipos) return -2 ;
   cc[k] = getfacc(fapt, pos, 1) ;
   if (hasmask[k]) ccmask[k] = getfacc(fapt, pos, 2) ;
   else ccmask[k] = '9' ;
  }

  t = 0 ; 
  if ((polarid != NULL) && (base2num(cc[polarindex]) < 0)) return -1 ;
  for (k=0; k<npops; ++k) { 
   if (isiub2(cc[k])) ++t ;
  }

  ++ncnt ;
  if (ncnt == 1) { 
   for (k=0; k<npops; ++k) { 
    fapt = fainfo[k] ;
    printfapt(fapt) ;
   }
  }

  return t ;
}

void readcommands(int argc, char **argv) 

{
  int i, t;
  phandle *ph ;
  char str[512]  ;
  int n, kode ;
  int pops[2] ;

  while ((i = getopt (argc, argv, "p:cvV?")) != -1) {


    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'V':
	verbose = YES ;
	break; 

      case 'v':
	printf("version: %s\n", version) ; 
	break; 

      case 'c':
	checkmode = YES ;
	break; 


      case '?':
		default:
		exit(usage());

      }
  }
         
   if (parname == NULL) exit(usage()) ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

        getstring(ph, "dbhetfa:", &iubfile) ;
        getstring(ph, "dbmask:", &iubmaskfile) ;
	if (! (iubfile && iubmaskfile)) {
		fatalx("Please give values to dbhetfa and dbmask in the parameter file.\n");
        }

   getstring(ph, "regname:", &regname) ;
   getint(ph, "pagesize:", &pagesize) ;
   getint(ph, "minfilterval:", &minfilterval) ;
   getint(ph, "allowmissing:", &allowmissing) ;
   getint(ph, "allowhets:", &allowhets) ;
   getstring(ph, "dbhetfa:", &iubfile) ;
   getstring(ph, "dbmask:", &iubmaskfile) ;
	t = NO;
   getint(ph, "transitions:", &t) ; if (t==YES) abxmode = 3 ;
	t = NO;
   getint(ph, "transversions:", &t) ; if (t==YES) abxmode = 2 ;
   getint(ph, "abxmode:", &abxmode) ; 
	getstring(ph, "minchrom:", &minch) ;
   	getstring(ph, "maxchrom:", &maxch) ;
	getstring(ph, "chrom:", &regname) ;
   getstring(ph, "polarize:", &polarid) ;
   getstring(ph, "indivname:", &indivname) ;
   getstring(ph, "indivoutname:", &indoutfilename) ; /* changed 11/02/06 */
   getstring(ph, "snpoutname:", &snpoutfilename) ; /* changed 11/02/06 */
   getstring(ph, "genooutname:", &genooutfilename) ; 
   getstring(ph, "genotypeoutname:", &genooutfilename) ; 
 
   t = 1 ; getint(ph, "lopos:", &lopos) ; lopos = MAX(lopos, t) ;
   t = BIGINT ; getint(ph, "hipos:", &hipos) ; hipos = MIN(hipos, t) ;
   getint(ph, "doxchrom:",  &doxchrom) ;
   getint(ph, "doychrom:",  &doychrom) ;
   getint(ph, "readfailOK:", &readfailOK) ;
   getint(ph, "domt:",  &domt) ;
   getint(ph, "minderiv:", &minderiv) ;
   getint(ph, "maxderiv:", &maxderiv) ;

   printf("### THE INPUT PARAMETERS\n");
   printf("##PARAMETER NAME: VALUE\n");
   writepars(ph);
   closepars(ph) ;
   fflush(stdout) ;
}

/*
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
*/
void clearfainfo(FATYPE *fapt, int mode)
{

 if (mode==1) {
   freestring(&fapt->faname) ;
   freestring(&fapt->alias) ;
   freestring(&fapt->famask) ;
   fapt -> fai = fapt -> faimask =  NULL ;
 }
 freestring(&fapt -> rstring) ;
 freestring(&fapt -> mstring) ;
 freestring(&fapt -> regname) ;
 fapt -> lopos = fapt -> hipos = -1 ;
 fapt -> len = 0 ;
 fapt -> mlen = 0 ;
 fapt -> rlen = 0 ;
}

void printasc(ASC *ascpt)  
{
 int k ;
 printf("asc: %d\n", ascpt -> ascnum) ;
 for (k=0; k<npops; ++k) { 
  printf("%12s %3d %3d\n", poplist[k], ascpt -> val[k], ascpt -> maxval[k]) ;
 }
 printnl() ;
 fflush(stdout) ;
}

void printfapt(FATYPE *fapt)
{
  int k, np = 0, mlen = fapt -> mlen   ;
  char cc ;

  printf("fapt: %s %s\n",  fapt->faname, fapt->alias) ;
	printf("%p %d %d ", fapt -> fai, fapt -> lopos, fapt -> hipos) ;


  printf("len: %d rlen: %d\n", fapt -> len, fapt -> rlen) ;
  printf("mask: %s %d\n", fapt -> famask, fapt -> mlen) ;
  printnl() ;

  fflush(stdout) ;
}
           

char fixval(char iub, char cm) 
{
  int t ;

  if (cm == CNULL) return iub ;
  if (cm == '-')  return 'N' ; 
  t =(int) (cm) - (int) '0'  ;
  if (t<minfilterval) return '?' ;
  return iub ;
}

int abx(int a, int b) 
{

 if (a<0) return -1 ;
 if (b<0) return -1 ;
 if (a==b) return -1 ;
 if (a>1) return abx(3-a, 3-b) ;

 if (a==0) return b-1 ; 
 if (b==0) return 3 ;  // CA
 if (b==2) return 4 ;  // CG
 if (b==3) return 5 ;  // CT

 fatalx("badbug\n") ;
 
}

int abxok(int abx, int abxmode) { 

 int t ;

 if (abxmode >= 10) { 
  t = abxmode - 10 ;
  if (abx == t) return YES ;
  return NO ;
 }
 
  
   
 switch (abxmode)  { 
 case 0:  
  return YES ;
 case 1:
  if (abx==0) return YES ; // AC
  if (abx==2) return YES ; // AT
  if (abx==3) return YES ; // CA
  return NO ;
 case 2: 
//No AG, CT transversions only
  if (abx==1) return NO ; 
  if (abx==5) return NO ; 
  return YES ; 
 case 3:  // transitions only
  if (abx==1) return YES ; 
  if (abx==5) return YES ; 
  return NO ;
 default:  
  fatalx("abxmode %d not implemented\n", abxmode) ;
 }
}

void ckopen(char *name)  
{
  FILE *fff ; 

  openit(name, &fff, "r") ;
  printf("open OK\n") ; 

  fclose(fff) ;

}
