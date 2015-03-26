/*
* cascertain.c: Pull down the SNPs that match the ascertain criterion.
* Author: Nick Patterson
* Revised by: Mengyao Zhao
* Last revise date: 2015-03-26
* Contact: mengyao_zhao@hms.harvard.edu
*/

#include <sys/wait.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <nicksam.h>
#include <getpars.h>  
#include "bam.h"
#include "faidx.h"
#include "globals.h" 
#include "popsubs.h"
#include "mcio.h"

typedef struct { 
 int *val ; 
 int *maxval ; 
 int ascnum ;
 int np ;
} ASC ;

char *table_path = NULL;
char *regname = NULL ; 
char *snpname = NULL ; 
char *iubfile = NULL;
char *iubmaskfile = NULL; 

char *parname = NULL ;
int  pagesize = 20*1000*1000 ;  // page size for getiub
int minfilterval = 1 ;

int minchrom = 1 ;
int maxchrom = 23 ;
char *minch = NULL;
char *maxch = NULL;
int xchrom = -1 ;

char *monoplistname = NULL ;
char **monosamplist  ;
int  monoval = -1 , nmonosamps = 0 ;

int lopos, hipos  ;
char **fasta ;  // in core bases
char *ascstring = NULL ;
char *noascstring = NULL ;

#define VERSION  "300"    

// monoplistname added  
void readcommands(int argc, char **argv) ;
int getfalist(char **poplist, int npops, char *dbfile, char **iublist) ; 
void clearfainfo(FATYPE *fapt, int mode) ;
char *myfai_fetch(faidx_t *fai, char *reg, int  *plen) ;
int loadfa(char **poplist, int npops, FATYPE ***pfainfo, char *reg, int lopos, int hipos)  ;
int minfalen(FATYPE **fainfo, int n)  ;


int istransition(char iubc)  ;
int abx(int a, int b)  ; 
int abxok(int abx, int abxmode) ;   
void copyfq(char **names, char **data, int *len, int a, int b) ;
int getiub(char *cc, char *ccmask, FATYPE **fainfo, char *reg, int pos) ; 

int setasc(char *ascstring, char **mlist, int nmlist, int mval)  ;
int setnoasc(char *ascstring) ;

int checkasc(ASC **asct, int nasct, char *cc, char *ccmask, char *pc1, char *pc2, char *regname, int pos) ;
void printasc(ASC *ascpt) ;
void printfapt(FATYPE *fapt);
int fvalid(char cm) ;
void prints(FILE *fff, int pos, char c1, char c2)  ;

char **poplist ; 
int *hasmask ;
int npops ;
int seed = 77  ;
int abxmode = 0 ;
int db = 1;	// Use .dblist

ASC **asctable ;
ASC **noasctable ;
int nasc, nonasc ;

static int usage() 
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   cascertain -p <parameter file> [options] \n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-d	directory of the data files (Please set this parameter, if you do not set .dblist files. If this parameter is used to give the data file location, .dblist files will not be used.)\n");
	fprintf(stderr, "\t-V	Print more information while the program is running.\n");
	fprintf(stderr, "\t-v	Show version information.\n");
	fprintf(stderr, "\t-? 	Show the instruction. (For detailed instruction, please see the document here: https://github.com/mengyao/cTools)\n\n");
	return 1;
}

int main(int argc, char **argv)
{

 char *spt[MAXFF], *zspt[MAXFF] ;
 int  k, t, t2 ;
 ASC *ascpt ;
 char *reg ;
 FATYPE **fainfo, *fapt ; 
 int xnpops ;
 int pos, chrom ;
 char *cc, *ccmask ;
 char c1, c2 ; 
 FILE  *fff ;
 char ss[20] ;

 hipos = 1000*1000*1000 ;
 lopos = 0 ;
 int numout = 0 ;
 int abxkode ;
 
 readcommands(argc, argv) ;

	if (minch != NULL) {
		if (minch[0] == 'X') minchrom = 23;
		else minchrom = atoi(minch);
		if (maxch[0] == 'X') maxchrom = 23;
		else maxchrom = atoi(maxch);
	}


 printf ("ascertain: %s reg: %s\n", ascstring, regname) ; 
 if (snpname != NULL) openit(snpname, &fff, "w") ;
 else fff = stdout ; 
 if (regname != NULL) { 
  if (regname[0] == 'X') xchrom = 23 ; 
  //else if (regname[0] == 'Y') xchrom = 24 ;
	//else if (!strcmp(regname, "MT")) xchrom = 25; 
  else xchrom = atoi(regname) ;
 }

 SRAND(seed) ;

 if (monoplistname != NULL) { 
  nmonosamps = numlines(monoplistname) ;
  ZALLOC(monosamplist, nmonosamps, char *) ;
  nmonosamps = loadlist(monosamplist, monoplistname) ;
  if (monoval < 0) { 
   printf("monoval not set: assuming 1 (derived)\n") ;
   monoval = 1 ;
  }
 }

	if (! ascstring) fatalx("No ascertain criterion. Please check your parameter file.\n");
 t = strlen(ascstring) ; 
 if (ascstring[t-1] == ';') ascstring[t] = CNULL ;
 nasc = setasc(ascstring, monosamplist, nmonosamps, monoval) ; // side effect set poplist npops
 nonasc = setnoasc(noascstring) ; // side effect set poplist npops

 xnpops = npops+1  ;

 printf("poplist:\n") ;
 printstrings(poplist, npops) ;
  
 for (k=0; k<nasc; ++k) { 
  ascpt = asctable[k] ;
  printasc(ascpt)  ;
 }
 if (nonasc>0) {
 printf("noasc:\n") ;
 for (k=0; k<nonasc; ++k) { 
  ascpt = noasctable[k] ;
  printasc(ascpt)  ;
 }}
 regname = strdup("22") ; 	// default chromosome 
 reg = regname ;
 ZALLOC(hasmask, npops+1, int) ;

 loadfa(poplist, npops, &fainfo, reg, lopos, hipos)  ;
 printf("npops: %d\n", npops) ;

  for (k=0;  k< npops; ++k) { 
   fapt = fainfo[k] ;
   printfapt(fapt) ;
  }

 ZALLOC(cc, xnpops, char) ; 
 ZALLOC(ccmask, xnpops, char) ; 

 cc[npops] = CNULL ;
 ccmask[npops-1] = CNULL ; // don't test chimp
//fprintf(stderr, "minchrom: %d\nxchrom: %d\nmaxchrom: %d\n", minchrom, xchrom, maxchrom);
 for (chrom = minchrom; chrom <= maxchrom; ++chrom) { 
  if ((xchrom > 0) && (xchrom != chrom)) continue ;
  sprintf(ss, "%d", chrom) ;
  if (chrom == 23) strcpy(ss, "X") ;
  //if (chrom == 24) strcpy(ss, "Y") ;
  //if (chrom == 25) strcpy(ss, "MT") ;
  freestring(&regname) ;
  regname = strdup(ss);
  reg = regname ;
//fprintf(stderr, "reg: %s\n", reg);

  for (pos = lopos ; pos <= hipos; ++pos) { 
  // t = getiub(cc, ccmask, fainfo, reg, pos)  ;
	t = getiub(cc, ccmask, fainfo, ss, pos)  ;  
   if (t==-5) break ;
   if (t<0) continue ;
    
   t = checkasc(asctable, nasc, cc, ccmask, &c1, &c2, regname, pos) ; 
   if (t==YES) {
    if (abxmode != 0) {
     abxkode = abx(base2num(c1), base2num(c2)) ;
     if (abxkode < 0) continue ;
     t = abxok(abxkode, abxmode) ;
     if (t==NO) continue ;
    }
    t2 = checkasc(noasctable, nonasc, cc, ccmask, &c1, &c2, regname, pos) ; 
    if (t2==NO) prints(fff, pos, c1, c2) ;  // printf("snp: %d %s %s\n", pos, cc, ccmask) ;
    if (verbose && (t2==NO)) printf("hit: %d %d %s %c%c\n", chrom, pos, cc, c1 , c2) ;
   }
 }}

 if (snpname != NULL) fclose(fff) ;
 
 printf("## end of cascertain\n") ;

 return 0 ;
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

int checkasc1(ASC *ascpt, int *cnt, int *cnt1) 
{
 int k, mx, x, v ; 
 for (k=0; k<npops; ++k) {
  mx = ascpt -> maxval[k] ;
  x = ascpt -> val[k] ;
  if (mx == 0) continue ;
  v = cnt[k] ;
  if (v == 3) return NO ;
  v = 2 - v ;  // number of derived; dangerous bend 
  if ((mx==2) && (x != v)) return NO ; 
  if ((mx==1) && (x != cnt1[k])) return NO ; 
 }
 return YES ;
}

void mkcnt1(int *cnt1, int *cnt2, int npops, char *cc, char cbase, char *regname, int pos) 
// deterministic hash -> mapdown hets   cb[x]  
{
  int k, v, t, x, z ; 
  char ss[100] ;
  char cb[2], chet, c, *pop ; 

  for (k=0; k<npops; ++k)  {
   v = cnt2[k] ;
   if (v==1) { 
    chet = cc[k] ;
    t = iubcbases(cb, chet) ; 
    if (t != 2) fatalx("logic bug(mkcnt1)\n") ;
    pop = poplist[k] ;
    sprintf(ss, "%s%s%d%c", pop, regname, pos, chet) ;
    x = abs(stringhash(ss)) ;
    x  = x % 2 ; 
    z = 0 ;
    if (cb[x] == cbase) z = 1 ; 
    cnt1[k] = 1-z ; 
    continue ;
   }
   cnt1[k] = (2-v)/2 ;
  }
}

int checkasc(ASC **asct, int nasct, char *cc, char *ccmask, char *pc1, char *pc2, char *regname, int pos) 
// regname, pos used for hash to find single allele of het
{
 int valid[20] ;
 int allelecnt[20], asc1[20] ;
 int cnum, t ;
 char cm, cx, cchimp, c1, c2 ; 
 int isdiff = NO, kasc, k ;
 ASC *ascpt ;

 if (nasct==0) return NO  ;
 cchimp = cc[npops-1] ;
 *pc1 = *pc2 = c1 = c2 = cchimp ;
 
 cnum = base2num(cchimp) ;
 if (cnum<0) return NO ;
 ivclear(valid, 1, npops) ;
 for (k=0; k<npops-1; ++k) { 
  cm = ccmask[k] ; 
  if (fvalid(cm) == NO) valid[k] = 0 ;
  cx  = cc[k] ;
  if (isiub2(cx) == NO)  valid[k] = 0 ;
  if ((valid[k] == YES) && (cx != cchimp)) isdiff = YES ;
 }
 if (isdiff == NO) return NO ;
// now check each ascertainment
 t = mkcnt(allelecnt, cc, ccmask, npops, &c1, &c2) ; 
 *pc1 = c1 ;  *pc2 = c2 ; 
 if (t<0) return NO ;
 if (c1==c2) return NO ;
 mkcnt1(asc1, allelecnt, npops, cc, c1, regname, pos) ;
 for (kasc=0; kasc<nasct; ++kasc) {
  ascpt = asct[kasc] ;
  t = checkasc1(ascpt, allelecnt, asc1) ;       
  if (t==YES) return YES ;
 }
 return NO ;
}

void setasct(ASC *ascpt, char *sx, char **pops, int npops) 
{
 char *spt[MAXFF] ;
 char *spt2[MAXFF] ;
 char *w1 ; 
 int *val, *maxval ;         
 int n1, n2, k, t ;

 if (strlen(sx)<3) fatalx("bad ascertainment string\n") ;
 ZALLOC(ascpt -> maxval, npops, int) ;
 ZALLOC(ascpt -> val, npops, int) ;
 val = ascpt -> val ; maxval = ascpt -> maxval ;

 w1 = strdup(sx) ; 
 substring(&w1, "::", ":") ;
 n1 = splitupx(w1, spt, MAXFF, ',') ;
 freestring(&w1) ;
 for (k=0; k<n1; ++k) { 
  n2 = splitupx(spt[k], spt2, MAXFF, ':') ; 
  if (n2 != 3) fatalx("bad ascertain: %s\n", sx) ;
  t = indxstring(pops, npops, spt2[0]) ;
  if (t<0) fatalx("bad pop. %s %s\n", spt2[0], sx) ;
  val[t] = atoi(spt2[1]) ;
  maxval[t] = atoi(spt2[2]) ;
  freeup(spt2, n2) ;
 }
 freeup(spt, n1) ;
}

void mkmstring(char *w1, char **mlist, int nmlist, int mval)  
{
  char ww[100] ;
  int k, t ;

  w1[0] = CNULL ;
  if (nmlist == 0) return ; 

  for (k=0;  k<nmlist; ++k) { 
   sprintf(ww, "%s::%d:2," ,mlist[k], 2*mval) ;
   strcat(w1, ww) ;
  }
  t = strlen(w1) ; 
  w1[t-1] = CNULL ;
  printf("zzmkmstring:\n") ;  printstring(w1, 50) ;
}

int setasc(char *ascstring, char **mlist, int nmlist, int mval) 
{  
#define MAXP 10
 char *spt[MAXFF] ;
 char *spt2[MAXFF] ;
 char *pops[1000] ; 
 char *w1, *sx ; 
 int np=0, t, n1, n2, k, nasc ;
 char *ww1,  *ww2 ; 

  t = 40*nmlist + strlen(ascstring) + 10 ; ; 
  ZALLOC(ww1, t, char) ; 
  ZALLOC(ww2, t, char) ;
  mkmstring(ww1, mlist, nmlist, mval) ;

// step 1 make list of pops
 w1 = strdup(ascstring) ; 
 substring(&w1, ";", ",") ;
 n1 = splitupx(w1, spt, MAXFF, ',') ;
 freestring(&w1) ;
 for (k=0; k<n1; ++k) { 
  n2 = splitupx(spt[k], spt2, MAXFF, ':') ; 
  if (n2==0) continue ;
  sx = spt2[0] ;
  t = indxstring(pops, np, sx) ; 
  if (t<0) {  
    pops[np] = strdup(sx); ++np ; 
  }
  freeup(spt2,n2) ;
 }
 freeup(spt, n1) ;
 if (np==0) fatalx("no pops!\n") ;

 ZALLOC(poplist, np+nmlist+1, char *) ;
 copystrings(pops, poplist, np) ;
 freeup(pops, np) ;
 copystrings(mlist, poplist+np, nmlist) ;
 np += nmlist ; 
 poplist[np] = strdup("Chimp") ;
 npops = np + 1 ;

 if (checkdup(poplist, npops)) { 
  printf("duplicate sample!\n") ;
  printstrings(poplist, npops) ;
  fatalx("no dups allowed\n") ;
 }

 nasc = splitupx(ascstring, spt, MAXFF, ';') ;
 ZALLOC(asctable, nasc, ASC *) ;
 for (k=0; k<nasc; ++k) { 
  ZALLOC(asctable[k], 1, ASC) ;
  strcpy(ww2,spt[k]) ;  
  if (nmlist>0) { 
   strcat(ww2, ",") ; 
   strcat(ww2, ww1) ;
  }
  setasct(asctable[k], ww2, poplist, npops) ;
  asctable[k] -> ascnum = k+1 ;
 }
 freeup(spt, nasc) ; 
 free(ww1) ;
 free(ww2) ;
 
 return nasc ;
}

int setnoasc(char *ascstring) 
{  
#define MAXP 10
 char *spt[MAXFF] ;
 int k, nonasc ;

 if (ascstring == NULL) return 0 ;

 nonasc = splitupx(ascstring, spt, MAXFF, ';') ; //Split string ascstring with. spt contains the suffix of the splited ascstring. Returns the number of fregments.	
 if (nonasc==0) return 0 ;
 ZALLOC(noasctable, nonasc, ASC *) ;
 for (k=0; k<nonasc; ++k) { 
  ZALLOC(noasctable[k], 1, ASC) ;
  setasct(noasctable[k], spt[k], poplist, npops) ;
  noasctable[k] -> ascnum = k+1 ;
 }
 freeup(spt, nonasc) ; 
 
 return nonasc ;
}

int loadfa(char **poplist, int npops, FATYPE ***pfainfo, char *reg, int lopos, int hipos) 
{
 static int k, numfalist, t, len ;
 static char **falist, **famasklist ;
 static FATYPE **fainfo, *fapt ;
 int *falen ;
 char *ttfasta ;
 int lo, hi ;
 static int ncall = 0 ;
  
  ++ncall ;

  if (ncall==1) {
   ZALLOC(falist, npops, char *) ;
   ZALLOC(famasklist, npops, char *) ;
	if (db == 0) {
	   numfalist = setfalist(poplist, npops, ".fa", falist) ;
	   t = setfalist(poplist, npops, ".filter.fa", famasklist) ;
	} else {
	   numfalist = getfalist(poplist, npops, iubfile, falist) ;	// set falist with the absolute path of hetfa files in .dblist file; falist contains the iubfile names
	   t = getfalist(poplist, npops, iubmaskfile, famasklist) ; 
	}

   if (numfalist != npops) {
      for (k=0; k<npops; ++k) { 
        if (falist[k] == NULL) printf("no fasta file for: %s\n", poplist[k]) ;
      }
      fatalx("Do not find the data files. Please use -d option or set dbhetfa and dbmask in your parameter file.\n") ;
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

    fapt -> faname = strdup(falist[k]) ; 	// faname is the hetfa file name
    fapt -> alias = strdup(poplist[k]) ;
    fapt -> famask = strdup(famasklist[k]) ;

    if (verbose) {
     printf("%s\n", fapt -> alias) ;
     printf("%s\n", fapt -> faname) ;
     printf("%s\n", fapt -> famask) ;
     printf("loading: %s\n", fapt -> faname) ;
     printnl() ;
   }
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

  for (k=0; k<numfalist ; ++k) {
     fapt = fainfo[k] ;
     ttfasta = myfai_fetch(fapt -> fai, reg, &len) ; 
     if (len==0) fatalx("bad fetch %s %s\n", fapt -> faname, reg) ;
      fapt -> rlen = len ;
      lo = MAX(1, lopos) ;
      hi = MIN(len, hipos) ;
      len = hi-lo + 1 ;
      ZALLOC(fapt -> rstring, len+1, char) ;
      strncpy(fapt -> rstring, ttfasta+lo-1, len) ; // indexing is base 1
      fapt -> rstring[len] == CNULL  ;
      freestring(&ttfasta) ;
      fapt -> regname = strdup(reg) ;
      fapt -> len = len ;
      fapt -> lopos = lo ;
      fapt -> hipos = hi ;

      if (hasmask[k] == NO)  continue ; 
      ttfasta = myfai_fetch(fapt -> faimask, reg, &len) ; 
      if (len==0) fatalx("bad fetch (mask)  %s %s\n", fapt -> faimask, reg) ;
      lo = MAX(1, fapt -> lopos) ;
      hi = MIN(len, fapt -> hipos) ;
      len = hi-lo + 1 ;
      ZALLOC(fapt -> mstring, len+1, char) ;
      strncpy(fapt -> mstring, ttfasta+lo-1, len) ; // indexing is base 1
      fapt -> mstring[len] = CNULL  ;
      freestring(&ttfasta) ;
      fapt -> mlen = len ;
  }

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
  fapt = fainfo[0] ; 
  lastreg = fapt -> regname ;
  strcpy(regbuff, reg) ;
  if (lastreg == NULL) lastreg = strdup("??") ;

  if (strcmp(lastreg, regbuff) != 0) {
   newpage = YES ;
   newreg =  YES; 
   fflush(stdout) ;
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
   fflush(stdout) ;
   freestring(&regname) ;

   regname = strdup(regbuff)  ;  
   lastreg = strdup(regbuff)  ;  
  }

  if (newpage == YES) { 
   fflush(stdout) ;
   lastlo = pos ;
   lasthi = pos + pagesize ;
   lastlo = MAX(lastlo, lopos) ;
   lasthi = MAX(lasthi, hipos) ;
   loadfa(poplist, npops, &fainfo, regname, lopos, hipos)  ;
  }

  for (k=0; k<npops; ++k) { 
   fapt = fainfo[k] ;
   if (pos<fapt -> lopos) return -2 ;
   if (pos>fapt -> hipos) return -2 ;
   cc[k] = getfacc(fapt, pos, 1) ;
   if (hasmask[k]) ccmask[k] = getfacc(fapt, pos, 2) ;
   else ccmask[k] = '9' ;
  }

  t = 0 ; 
  if (base2num(cc[npops-1]) < 0) return -1 ;
  for (k=0; k<npops; ++k) { 
   if (isiub2(cc[k])) ++t ;
  }

  ++ncnt ;
  if (ncnt == 1) { 
   printf("zz pos: %s %s\n", cc, ccmask) ;
   for (k=0; k<npops; ++k) { 
    fapt = fainfo[k] ;
    printfapt(fapt) ;
   }
  }

  ccmask[npops-1] = CNULL ;

  return t ;
}

void readcommands(int argc, char **argv) 
{
  int i, t = NO;
  phandle *ph ;
  char str[512]  ;
  int n, kode ;
  int pops[2] ;

  while ((i = getopt (argc, argv, "p:d:vV?")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'd':
	{
		char* p;
		table_path = strdup(optarg) ;
		p = strrchr(table_path, '/');
		if (!p || strcmp(p, "/")) {
			table_path = (char*)realloc(table_path, 256);
			table_path = strcat(table_path, "/");
		}
		db = 0;	// Don't use .dblist
	}
	break;

      case 'V':
	verbose = YES ;
	break; 

      case 'v':
	printf("version: %s\n", VERSION) ; 
	break; 

      case '?':
	default:
	exit(usage());
      }
  }
         
   if (parname == NULL) //return ;
		exit(usage());
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;
	
	if (db == 1) {
	   getstring(ph, "dbhetfa:", &iubfile) ;
	   getstring(ph, "dbmask:", &iubmaskfile) ;
		if (! (iubfile && iubmaskfile))
			fprintf(stderr, "Please use -d option to specify the directory of hetfa and mask files.\nAlternatively, please give values to dbhetfa and dbmask in the parameter file.\n");
	}
//fprintf(stderr, "here, readcommands\n");	
   getstring(ph, "regname:", &regname) ;
   getstring(ph, "snpname:", &snpname) ;
   getint(ph, "pagesize:", &pagesize) ;
   getint(ph, "minfilterval:", &minfilterval) ;
   getint(ph, "seed:", &seed) ;
   getstring(ph, "ascertain:", &ascstring) ;
   getstring(ph, "noascertain:", &noascstring) ;
//printf(stderr, "ascstring: %s\n", ascstring);	
   getint(ph, "transitions:", &t) ; if (t==YES) abxmode = 3 ;
   getint(ph, "transversions:", &t) ; if (t==YES) abxmode = 2 ;
   getint(ph, "abxmode:", &abxmode) ; 
  // getint(ph, "minchrom:", &minchrom) ;
  // getint(ph, "maxchrom:", &maxchrom) ;
   getstring(ph, "minchrom:", &minch) ;
   getstring(ph, "maxchrom:", &maxch) ;
  // getint(ph, "chrom:", &xchrom) ;
   getstring(ph, "chrom:", &regname) ;

   getstring(ph, "monosamples:", &monoplistname) ;
   getint(ph, "monoval:", &monoval) ;
 
   t = 1 ; getint(ph, "lopos:", &lopos) ; lopos = MAX(lopos, t) ;
   t = BIGINT ; getint(ph, "hipos:", &hipos) ; hipos = MIN(hipos, t) ;

   printf("### THE INPUT PARAMETERS\n");
   printf("##PARAMETER NAME: VALUE\n");
   writepars(ph);
   closepars(ph) ;
   fflush(stdout) ;
}

int setfalist(char **poplist, int npops, char *dbfile, char **iublist) {
	int t;
	for (t = 0; t < npops; ++t) {
		iublist[t] = strdup(table_path);
		iublist[t] = (char*) realloc(iublist[t], 64);
		iublist[t] = strcat(iublist[t], poplist[t]);
		if ((!strcmp (poplist[t], "Chimp") || !strcmp (poplist[t], "Href")) && strcmp (dbfile, ".fa")) {
			free (iublist[t]);
			iublist[t] = "NULL";
		} else 
			iublist[t] = strcat(iublist[t], dbfile);
	}
	return npops;
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

  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, MAXFF) ; 
   if (nsplit<1) continue ;
   sx = spt[0] ;
   if (sx[0] == '#') { 
    freeup(spt, nsplit) ;
    continue ;
   }
   t = indxstring(poplist, npops, spt[0]) ; 
   if (t<0) { 
    freeup(spt, nsplit) ; 
    continue ;
   }
    sx = spt[2] ;
    iublist[t] = strdup(sx) ;
    freeup(spt, nsplit) ;
    ++nx ;
  }

   fclose(fff) ;
   return nx ;
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
  printf("%p %d %d n", fapt -> fai, fapt -> lopos, fapt -> hipos) ;
  printf("len: %d rlen: %d\n", fapt -> len, fapt -> rlen) ;
  printf("mask: %s %d\n", fapt -> famask, fapt -> mlen) ;
  printnl() ;

  fflush(stdout) ;
}
           
int fvalid(char cm) 
// is cm indicating valid? 
{
  int t ; 

  if (minfilterval<0) return YES ;
  t = (int) cm - (int) '0' ;

  if (t<0) return NO;
  if (t>=10) return NO ;
  if (t<minfilterval) return NO ;
  return YES ; 
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

int mkcnt(int *cnt, char *iub, char *mask, int npops, char *pc1, char *pc2)  
// number of copies of ref allele
{
  char c1, c2, cm, x1, x2, cbases[2] ; 
  int i, t, x, ok ;
  
  c1 = *pc1; c2 = *pc2 ;

  for (i=0; i<npops; ++i) { 
   iub[i] = fixval(iub[i], mask[i]) ;
   t = iubcbases(cbases, iub[i]) ;
   if (t<0) { 
     cnt[i] = 3  ; iub[i] = '?' ; continue ;
   }
   x1 = cbases[0] ; x2 = cbases[1] ; 
   if ((c1==c2)  && (x1 != c1)) c2 = x1 ;
   if ((c1==c2)  && (x2 != c1)) c2 = x2 ;
   if ((x1 != c1)  && (x1 != c2)) return -1 ; 
   if ((x2 != c1)  && (x2 != c2)) return -1 ; 
   if ((x1==c1) && (t==1)) { 
     cnt[i] = 2; 
     continue ;
   }
   if ((x1==c2) && (t==1)) { 
     cnt[i] = 0; 
     continue ;
   }
   if (t==1) fatalx("logic bug\n") ;
   cnt[i] = 1 ;
  }
  *pc1 = c1 ; *pc2 = c2 ;

  t = 0 ; 
  for (i=0; i< npops-2; ++i) { 
   if (cnt[i] < 3) ++t ;
  }
  return t ;
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
