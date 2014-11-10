#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>

#include <nicklib.h>  
#include <getpars.h>
#include "popsubs.h"
#include "nicksam.h"

#define MAXSTR 2000 
#define MAXFF 100
#define SSIZE 50 

extern int verbose ;
extern char *oldxwstring ; 
extern char **poplist ;
extern char **pophdr ;
extern char *dlist ;
extern int *threshlo, *threshhi, *quallo ;
extern int neanpop, xwpop, npops ;
extern int printtable ;

static int pmdscore = NO ;

static char ***bonetable ;
static int numbline = 0 ; 

static char **bonelist = NULL, **librarylist = NULL ; 
static int numbones, numlibrary, numbadbones= 0  ;
static char **badbonelist ; 

static char ***rdtable ;
static double ***qtable ;
static int numrline = 0 ; 
static char **rdhdrlist, **rdliblist ; 
static char **rdsamplist, **rdpoplist ; 
static int numrhdrs = 0, numrlibs, numrsamps, numrpops ;
static int qhack = NO ;  
int ishack ;
static phandle *ph  = NULL;

extern char **rglist ;
extern int *rgpops, numrg ;

extern char **iublist ;
extern int *iubpops, numiub ;

#define MAXBADRGS 1000 
static  int numbadrgs = 0 ;
static int maxbadrgs = MAXBADRGS ;
static char  *badrgs[MAXBADRGS] ;
static int nopoplist = NO ;

void printbadrgs() 
{
  int k ;

  for (k=0; k<numbadrgs; ++k) { 
   printf("badrgs: %30s\n", badrgs[k]) ;
  }
}

void freebonetable() 
{

   int i, j ;
   for (i=0; i<3; i++) { 
    for (j=0; j<numbline; ++j) { 
     free(bonetable[i][j])  ;
    }
    if (bonetable[i]!= NULL) free(bonetable[i]) ;
   }
   if (bonetable != NULL) free(bonetable) ;
   numbline = 0 ;

   if (numbones>0) freeup(bonelist, numbones) ; 
   if (numlibrary>0) freeup(librarylist, numlibrary) ; 

   numbones = numlibrary = 0 ;
   librarylist = bonelist = NULL ;  // 
}

char *num2h(int k) 
{
// returns hdr  
 if (k>=numbline) return NULL ;
 return bonetable[0][k] ;
}
int loadstr(char **strlist, char **xlist, int numin) 
{
  int k, t, n=0 ;
  char *ss ;

  for (k=0; k<numin; ++k) {  
   ss = strlist[k] ;
   t = indxstring(xlist, n, ss) ; 
   if (t >= 0) continue ;
   xlist[n] = strdup(ss) ; 
   ++n ;
  }
  return n ;
}


int 
getxstr(char ***xx, int maxrow, int numcol, char *fname)
{

  char line[MAXSTR] ;
  char *spt[MAXFF] ;
  char *sx ;
  int nsplit, i, j, num=0, maxff ;
  FILE *fff ;
  static int nbad = 0 ; 

  if (fname == NULL) fff = stdin ; 
  else {
   openit(fname, &fff, "r") ;
  }
  maxff = MAX(MAXFF, numcol) ; 

  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, maxff) ; 
   sx = spt[0] ;
   if (sx[0] == '#') {
    freeup(spt, nsplit) ;
    continue ;
   }
   if (nsplit<numcol) { 
     ++nbad ;
     if (nbad<10) printf("+++ bad line:\n%s\n", line) ;
     continue ;
   }
   if (num>=maxrow) fatalx("too much data\n") ;
   for (i=0; i<numcol; i++)  {
    xx[i][num]  = strdup(spt[i]) ;
   }
   freeup(spt, nsplit) ;
   ++num ;
  }
  if (fname != NULL) fclose(fff) ;
  return num ;
}
int indxsub(char *ss, char **table, int tablen) 
{
  int k ;

  for (k=0; k<tablen; ++k) { 
   if (strstr(ss, table[k]) != NULL) return k ;
  }
  return -1 ;
}


int checkm(char *sa, char *sb) 
{
 int t ;
 if (sa == NULL) return YES ;
 if (sb == NULL) return NO ;
 t = strcmp(sa, sb) ; 
 if (t != 0) return NO ;
 return YES ;
}

/**
char num2base(int k)
{
    char *bases = "ACGT" ;

    if (k<0) return '?' ;
    if (k>3) return '?' ;
    
    return bases[k] ;
}

int base2num(char c)
{
    char cc ; 

    cc = toupper(c) ;

    switch (cc)  {
     case 'A':  return 0;
      break ;
     case 'C':  return 1;
      break ;
     case 'G':  return 2;
      break ;
     case 'T':  return 3;
      break ;
     default:  return -1 ;
    }
}
*/

int kodex (int *aa, int len)  
{
 int i, t ;
 t = 0; 

 for (i=0; i<len; ++i) { 
  t = t*4 + aa[i] ;
 }
 return t ;
}

int kodexb (int *aa, int len, int base)  
{
 int i, t ;
 t = 0; 

 for (i=0; i<len; ++i) { 
  t = t*base + aa[i] ;
 }
 return t ;
}


int kode4 (int *aa)  
{
 return kodex(aa, 4) ;
}
void dekodex (int kode, int *aa, int len)  
{
 int i, t ;

 t = kode ; 
 for (i=len-1; i>=0; --i) { 
  aa[i] = t % 4 ;
  t /= 4 ;
 }
}
void dekodexb (int kode, int *aa, int len, int base)  
{
 int i, t ;

 t = kode ; 
 for (i=len-1; i>=0; --i) { 
  aa[i] = t % base ;
  t /= base ;
 }
}


void dekode4 (int kode, int *aa)  
{
 dekodex(kode, aa, 4) ;
}

void printmat0(double *a, int m, int n) 
{
  int *aa ;  

  ZALLOC(aa, m*n, int) ;  
  fixit(aa, a, m*n) ;

  printimatl(aa, m, n) ;
  free(aa) ;
}

void loadbaseprobs(double **baseprobs, int *basevalid, PILEUP *pilept, int rind)
{

//OOD
   int d, x, k, n, t ;
   BASEP *basept ;

   clear2D(&baseprobs, MAXPOP+1, 4, 0.0) ;  
   ivzero(basevalid, MAXPOP+1) ;

   if (rind<0) return ; 
   baseprobs[MAXPOP][rind] = 1.0 ;
   basevalid[MAXPOP] = YES ;

   for (x=0; x<MAXPOP; ++x) {  
    basept = pilept -> bp[x] ;
    d = basept -> depth ; 
    d = MIN(d, MAXD) ; 
    if (d<=0) continue ; 
    n = 0 ;
    for (k=0; k<d; ++k) {  
     t = base2num(basept -> bases[k]) ;
     if (t<0) continue ;
     baseprobs[x][t] += 1.0 ;
     ++n ;
    }
    if (n==0) continue ;
    bal1(baseprobs[x], 4) ;
    basevalid[x] = 1 ;
   }
}
int loadpprobs(double **baseprobs, int *basevalid, 
  double *pprobs, int *pkode, int *indx, int indlen)  
{

  double pp[256], ppp[256] ; 
  int pk[256], ppk[256] ;

  double pp1[4], y ; 
  int pk1[4] ;
  double *pppt ;
  int *pkpt ;
  int x, s, k, n ;

  if (verbose) printf("yyx %d\n", indlen) ;

  if (indlen == 0) return 0 ;
  if (indlen == 1) { 
   x = indx[0] ; 
   if (verbose) printf("ww %d %d\n", x, basevalid[x]) ;
   if (basevalid[x] == NO) return 0 ;
   s = 0 ;
   for (k=0; k<4; ++k) { 
    y = baseprobs[x][k] ; 
    if (y < 1.0e-10) continue ; 
    pprobs[s] = y ;
    pkode[s] = k ;
    ++s ;
   }
   return s ;
  }
  n = loadpprobs(baseprobs, basevalid, pp, pk, indx, indlen-1) ;
  if (verbose) printf("yy1 %d %d\n", indlen, n) ;
  if (n==0) return 0 ;
  s = loadpprobs(baseprobs, basevalid, pp1, pk1, indx+indlen-1, 1) ;
  if (verbose) { 
    printf("yy2 %d %d %d \n", indlen, 1, s) ;
  }
  if (s==0) return 0 ;
  if ((n*s) == 1) { 
   pprobs[0] = pp[0]*pp1[0] ; 
   pkode[0] = 4*pk[0] + pk1[0] ;
   if (verbose) printf("yy3 %d\n", 1) ;
   return 1 ;
  }
  pppt = ppp ;
  pkpt = ppk ;
  for (k=0; k<s; ++k) { 
   vst(pppt, pp, pp1[k], n) ;
   ivst(pkpt, pk, 4, n) ;
   ivsp(pkpt, pkpt, pk1[k], n) ;
   pppt += n ;
   pkpt += n ;
  }
  copyarr(ppp, pprobs, n*s) ;
  copyiarr(ppk, pkode, n*s) ;
  if (verbose) printf("yy3b %d\n", n*s) ;
  return n*s ;
}

void loadfilebase(char *parname) 
{

  ph = openpars(parname) ;
  dostrsub(ph) ;

}

/**

A1: /groups/reich/nick/broaddata/neander/analysis1 
Yorubaref: A1/NA18507-half-chimp.fa 
Koreanrref: A1/SJK-chimp.fa
Href:  A1/hs2pt.fa // href -> chimp
Chimp: A1/panTro2.fa
Map20:  A1/map20.allchr.fa
Eichler = A1/eichler_safe.fa
*/

void unloadfilebase() 
{

 closepars(ph) ; // sets ph NULL 


}

char *getfaname(char *sname) 
{

static char **fanames = NULL ;
static char **snames = NULL, *sss ;
char bname[256] ;

///char *dir = "/nfs/eva/nickp/data/jun09" ;
char *dir = "/groups/reich/nick/broaddata/neander/analysis1" ;

char *faname ;
int t ;

 if (ph == NULL) { 
  fatalx("(getfaname) loadfilebase not called\n") ;
 }
 
 sss = NULL ; 
 strcpy(bname, sname) ;
 strcat(bname, ":") ;
 getstring(ph, bname, &sss) ;
 if (sss == NULL) fatalx("(getfaname) %s not found\n", sname) ;
 printf ("## getfaname: %s %s\n", sname, sss) ;
 return sss ;

}
int patkode(int a, int b, int c, int d) 
{

 if (a<0) return -1 ;
 if (b<0) return -1 ;
 if (c<0) return -1 ;
 if (d<0) return -1 ;

 if ((a==b) && (a == c) && (a == d)) return 0 ;  // AAAA
 if ((a==b) &&  (a == d)) return 1 ;  // AABA
 if ((a==c) &&  (a == d)) return 2 ;  // ABAA
 if ((a==d) &&  (b == c)) return 3 ;  // ABBA
 if ((a==d) ) return 4 ;   //  ABCA
 if ((b==c) &&  (b == d)) return 5 ;  // BAAA
 if ((a==c) &&  (b == d)) return 6 ;  // BABA
 if ((b==d)) return 7 ; // BACA
 if ((a==b) &&  (c==d)) return 8 ;  // BBAA
 if ((a==b) &&  (a==c)) return 9 ;  // BBBA
 if ((a==b)) return 10 ; // BBCA 
 if ((c==d)) return 11 ; // BCAA
 if ((a==c)) return 12 ; // BCBA 
 if ((b==c)) return 13 ; // CBBA 
 return 14 ; // BCDA

}

char *kodepstring(int kode) 

{ 

 switch (kode) { 

 case 0: return "AAAA" ;
 case 1: return "AABA" ;
 case 2: return "ABAA" ;
 case 3: return "ABBA" ;
 case 4: return "ABCA" ;
 case 5: return "BAAA" ;
 case 6: return "BABA" ;
 case 7: return "BACA" ;
 case 8: return "BBAA" ;
 case 9: return "BBBA" ;
 case 10: return "BBCA" ;
 case 11: return "BCAA" ;
 case 12: return "BCBA" ;
 case 13: return "CBBA" ;
 case 14: return "BCDA" ;
 
 }
return "????" ;

}


int checkag(int *aa, int len) 
// AG or TC forbidden
{
 int cc[4] ;
 int k, x, t ;

 ivzero(cc, 4) ;
 for (k=0; k<len; ++k)  { 
  x = aa[k] ;
  if (x<0) continue ;
  ++cc[x] ;
 }
 t = cc[0]*cc[2] ;  if (t>0) return -1 ;
 t = cc[1]*cc[3] ;  if (t>0) return -1 ;
 return 1 ;
}

int abpat(int *ww, int *aa, int len, int rind) 
{

  int i, var = -1, ret = 0 ;

  ivzero(ww, len) ;
  
  for (i=0; i<len; ++i) { 
   if (aa[i] == rind) { 
    continue ;
   }
   ww[i] = 1 ;
   if (var < 0) var = aa[i] ;
   if (var != aa[i]) ret = -1 ;
  }
  return ret ;
}
int crackhdr(char *hdr, char *lib, char *samp, char *pop) 
{
   char *px ; 
   int t ;
   char ss[1000] ;
   char *scolon ;

    
   px = hdr ;


   strcpy(ss, px)  ;
   scolon = strchr(ss, ':') ;
   if (scolon==NULL) fatalx("bad read: %s\n%s\n", hdr, ss) ;
   scolon = strchr(scolon+1, ':') ;
   if (scolon==NULL) fatalx("bad read: %s\n%s\n", hdr, ss) ;
   *scolon = CNULL ;

   if (rdtable[0] == NULL) fatalx("no readtable\n") ;
   t = indxstring(rdtable[0], numrline, ss) ;  

   if (t<0) {
    fatalx("hdr not found %s\n", hdr) ;
   }

   if (lib != NULL) strcpy(lib,  rdtable[1][t]) ;
   if (samp != NULL) strcpy(samp, rdtable[2][t]) ;
   if (pop != NULL) strcpy(pop,  rdtable[3][t]) ;

   return t ;


}

int setpmdscore(int val) 
{
  pmdscore = val ;
//printf("pmdscore: %d\n", val) ;
}

int setrdtable(char **poplist, int npops, char *rdname) 
{
   char sss[10] ;
   int k, n=0, x, t, i ;
   static int ncall = 0 ;

   ++ncall ;  
   if (rdname == NULL) return 0 ;
   if (ncall==1) loadrdtable(rdname) ;

   x = 0 ;
   for (k=0; k< numrline; ++k) { 
    t = indxstring(poplist, npops, rdtable[3][k]) ;
    if (t<0) {
     freestring(&rdtable[3][k]) ;
     rdtable[3][k] = strdup("Ignore") ;
    }
    else ++x ;
   }  
   for (i=0; i< numrline; i++) { 
    if (!printtable) break ;
    printf("rtable: %4d ", i) ;
    printf("%40s ", rdtable[0][i]) ;
    printf("%10s ", rdtable[1][i]) ;
    printf("%10s ", rdtable[2][i]) ;
    printf("%15s ", rdtable[3][i]) ;
    printnl() ;
   }

   return numrline ;

}


int  loadlist(char **list, char *listname)   
// listname is just a list of names ... 
// dup check made
{
  int n, a, b, t  ;  
  n =  getss(list, listname) ;
  for (a=0; a<n; ++a) {
   for (b=a+1; b<n; ++b) {
    t = strcmp(list[a], list[b]) ;
    if (t==0) {
     fatalx("duplicate in list: %s\n", list[a]) ;
    }
   }
  }
  return n ;
}

char *abxstring(int abx) 
{
  if (abx == 0) return "--" ;
  if (abx == 2) return "tv" ;
  if (abx == 3) return "ti" ;
  if (abx == 10) return "AC" ;
  if (abx == 11) return "AG" ;
  if (abx == 12) return "AT" ;
  if (abx == 13) return "CA" ;
  if (abx == 14) return "CG" ;
  if (abx == 15) return "CT" ;

  return "??" ;

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

int calcvar(PILEUP *pilept, int *vnum, int nean) 
{

 int ccc[4] ; 
 int a, b, ta, tb, rind ; 
 char refbase ;

 *vnum = -1 ; 
 cntallelesnonean(pilept, ccc, nean) ;
 refbase = pilept -> refbase ; 

 rind = base2num(refbase) ; 
 if (rind<0) fatalx("badbug\n") ;
 ccc[rind] = 0 ; 

 if (intsum(ccc,4) == 0) return 0 ;
 ivlmaxmin(ccc, 4, &a, NULL) ;
 ta = ccc[a] ;
 ccc[a] = -99 ;
 
 ivlmaxmin(ccc, 4, &b, NULL) ;
 tb = ccc[b] ;
 if (ta==tb) {
  *vnum = -2 ;
  return 999 ;
 }
 *vnum = a ; 
 
 ccc[rind] = ccc[a] = 0 ;
 return intsum(ccc, 4) ;
}


void cntalleles(PILEUP *pilept, int *ccc) 
{

  cntallelesnonean(pilept, ccc, -99) ;

}

void cntallelesnonean(PILEUP *pilept, int *ccc, int nean) 
{
  BASEP *basept ;
  int x, d, k, t ;

  ivzero(ccc,4) ; 

  for (x=0; x<MAXPOPX; ++x) {  
   if (x==nean) continue ;
   basept = pilept -> bp[x] ;
   if (basept == NULL) return ; 
   d = basept -> depth ; 
   for (t=0; t<d; ++t) { 
    k =  base2num(basept -> bases[t]) ; 
    if (k>=0) ++ccc[k] ;
   }
 }
}

void forcebi(PILEUP *pilept, int xa, int xb) 
{
  BASEP *basept ;
  int j, x, d, k, t ;
  int ccc[4] ;

  ivzero(ccc,4) ; 

  if (xa >=0) ccc[xa] = 1 ;
  if (xb >=0) ccc[xb] = 1 ;

  if (xa>=4) fatalx("badbug\n") ;
  if (xb>=4) fatalx("badbug\n") ;

  for (j=0; j<MAXPOPX; ++j) {  
   basept = pilept -> bp[j] ;
   if (basept == NULL) continue ;
   d = basept -> depth ; 
   x = 0 ;
   for (t=0; t<d; ++t) { 
    k =  base2num(basept -> bases[t]) ; 
    if (k<0) continue ;
    if (ccc[k] == 1) { 
     copybb(basept, x, t) ; 
     ++x ;
    }
   }
   basept -> depth = x ;
 }
 return ;
}
void copybb(BASEP *basept, int a, int b) 
// copy b to a call was buggy (Kasia)
{
 
  if (a==b) return ;
  if (a<0) return ; 
  if (b<0) return ; 
  if (a>=MAXD) return ; 
  if (b>=MAXD) return ;

  basept -> bases[a] = basept -> bases[b] ; 
  
  basept -> qual[a] = basept -> qual[b] ; 
  basept -> diff[a] = basept -> diff[b] ; 
  basept -> strand[a] = basept -> strand[b] ; 
  basept -> bonenumber[a] = basept -> bonenumber[b] ; 
  basept -> libnumber[a] = basept -> libnumber[b] ; 
  basept -> mapq[a] = basept -> mapq[b] ; 

  return ;

}

int checktri(int *ccc, int rind) 
{

  int ddd[4], t, a ; 

  copyiarr(ccc, ddd, 4) ;
  for (a=0; a<4; a++) {  
   if (ddd[a] > 0) ddd[a] = 1 ;
  }
  t = intsum(ddd, 4) ; 
  if (t  > 2) return YES ;
  return NO ;

}

int calcvind(int *ccc, int rind, int *pvind) 
{

  int ddd[4], vind, a ; 

  vind = *pvind = -1 ;

  copyiarr(ccc, ddd, 4) ;
  ddd[rind] = 0  ;

  if (intsum(ddd, 4) == 0) return YES ;

  for (a=0; a<4; ++a) { 
   if (ddd[a] == 0) continue ; 
   if (vind >=0) return NO ;
   vind = a ;
  }
  *pvind = vind ;
  return YES ;
}

int checkicx(int *ccc, int rind, int *pvind) 
{

  int ddd[4], vind, a ; 

  copyiarr(ccc, ddd, 4) ;
  if (ddd[rind] < 2) return NO ; 
  ddd[rind] = 0  ;

  if (intsum(ddd, 4) == 0) return NO ;
  vind = -1 ;
  for (a=0; a<4; ++a) { 
   if (ddd[a] == 0) continue ; 
   if (vind >=0) return NO ;
   vind = a ;
  }
  if (ddd[vind] < 2) return NO ;
  *pvind = vind ;
  return YES ;
}

int loadrdtable(char *rdname)  
{
   int t, i ;

   if (rdname == NULL) fatalx("no read table\n") ;

   t = numlines(rdname) ;
   ZALLOC(rdtable, 4, char **) ;

   for (i=0; i<4; i++) { 
    ZALLOC(rdtable[i], t, char *) ;
   }

   numrline = getxstr(rdtable, t, 4, rdname) ;

   ZALLOC(rdhdrlist, numrline, char *) ;
   ZALLOC(rdliblist, numrline, char *) ;
   ZALLOC(rdsamplist, numrline, char *) ;
   ZALLOC(rdpoplist, numrline, char *) ;

   numrhdrs = loadstr(rdtable[0], rdhdrlist, numrline) ;
   numrlibs = loadstr(rdtable[1], rdliblist, numrline) ;
   numrsamps = loadstr(rdtable[2], rdsamplist, numrline) ;
   numrpops = loadstr(rdtable[3], rdpoplist, numrline) ;

   return numrline ;
}
void freerdtable() 
{

   int i, j ;
   for (i=0; i<4; i++) { 
    for (j=0; j<numrline; ++j) { 
     free(rdtable[i][j])  ;
    }
    if (rdtable[i]!= NULL) free(rdtable[i]) ;
   }
   if (rdtable != NULL) free(rdtable) ;
   numrline = 0 ;

   if (numrhdrs>0) freeup(rdhdrlist, numrhdrs) ; 
   if (numrlibs>0) freeup(rdliblist, numrlibs) ; 
   if (numrsamps>0) freeup(rdsamplist, numrsamps) ; 

   numrhdrs = numrlibs = 0 ;
   rdhdrlist = rdliblist = NULL ;  // 
}
int checkrg(char **rawp) 
{

 int isok = YES, k, t, tt ;
 char *sx ;

 for (k=1; k<numrg; ++k) { 
  sx = rglist[k] ;
  t = indxstring(rglist, k, sx) ;
  if (t>=0) { 
   tt = strcmp(rawp[rgpops[t]], rawp[rgpops[k]]) ; 
   if (tt != 0) isok = NO ; 
   printf("rgdup: %s  samples: %s  %s\n", sx, rawp[rgpops[t]], rawp[rgpops[k]]) ;
  }
 }

 return isok ;

}

int whichpop(char *rgs) 
{
// rgs is rgstring 
   char *s2, *spt[MAXFF], *sx, *sy, *sz ;  
   char *scolon ; 
   char bigb[MAXLEN] ; 
   
   int t, nsplit, k, s ;
   static int ncall = 0 ;

   ++ncall ;
   if (ncall > 100000) ncall = 100000 ;
   if (nopoplist) return 0 ;
   nsplit = splitupxbuff(rgs, spt, MAXFF, ':', bigb, MAXLEN) ;
   if (nsplit != 3) { 
    return -1 ;
   }
   sx = spt[2] ; 
   t = indxstring(rglist, numrg, sx) ;

   if (t<0) { 
    s = indxstring(badrgs, numbadrgs, sx) ; 
    if ((s<0) && (numbadrgs<(maxbadrgs-1))) { 
     badrgs[numbadrgs] = strdup(sx) ;
     ++numbadrgs; 
    }
    return -2 ;
   }
   return rgpops[t] ;
}

   
int whichhdr(char *s1)   
{

   return crackhdr(s1, NULL, NULL, NULL) ;

}

int whichlib(char *s1)   
{
// s1 is header 
   char *s2, sss[6] ;
   char lib[SSIZE] ;
   
   int t ;

   crackhdr(s1, lib, NULL, NULL) ;

   if (strlen(lib) == 0) return -2 ;
  
  if (numrlibs <= 0) return -1 ;
   t = indxstring(rdliblist, numrlibs, lib) ;
   if (t<0) return -1  ; 
   return t ;
}

int getreadx(READ *readpt, char *pline) 
{
 char *spt[MAXFF]  ;
 char line[MAXSTR]  ;
 int nsplit ;
 char *reg, c1, c2, *s1, *s2 ; 
 int pos, len, k ;
 char lib[MAXLEN], samp[MAXLEN], pop[MAXLEN] ;
 static int nbug = 0 ;

  line[MAXSTR] = '\0' ;
  cleanread(readpt) ;
  if (fgets(line, MAXSTR, stdin) != NULL)  {
   nsplit = splitup(line, spt, MAXFF) ; 
   if (nsplit == 0) return -99 ;

   strcpy(readpt -> qname, spt[0]) ;
   strcpy(pline, line) ;

   freeup(spt, nsplit) ;
   return 1 ;
  }
  return -98;
}
void setqhack(int val) 
{
 qhack = val ;
}
int hackit(char *str) 
{
 char *pt ; 
 if (qhack == NO) return NO ; 
 pt = strstr(str, "_p54.0021") ; 
 if (pt == NULL) return NO ; 
 return YES ;
}
int dsindex(char **spt, int nsplit) 
{

  int t, x ;

  if (pmdscore == NO) return -1 ;
  for (x=nsplit-1; x>=0 ; --x) { 
   t = strncmp(spt[x], "DS:Z:", 5) ;
   if (t==0) return x ;
  }

  return -1 ; 

}

int rgindex(char **spt, int nsplit) 
{

  int t, x ;

  for (x=nsplit-1; x>=0 ; --x) { 
   t = strncmp(spt[x], "RG:Z:", 5) ;
   if (t==0) return x ;
  }

  return -1 ; 


}

int readtonl(FILE *fff) 
{
  char c ;

  for (;;) { 
   c = fgetc(fff) ; 
   if (c==EOF) break ;
   if (c==CNL) break ;
  }
  return c ;
}
void printlong(char *line) { 
 static int nbad = 0 ;

 ++nbad ;
 if (nbad>100) return ;

  printf("*** warning:   excessively long read; start is:\n") ;
  printstring(line, 50) ;

}

int getreadfnorg(READ *readpt, FILE *fff) 
{
 char line[MAXSTR+1] ;
 char *spt[MAXFF]  ;
 int nsplit ;
 int  t, ipop, x ;
 char *reg, c1, c2, *s1, *s2 ; 
 int pos, len, k ;
 char  pop[MAXLEN] ;
 char bigb[MAXLEN] ; 
 static int nbug = 0 ;
 static int nbad = 0 ;
 int *pa, *pb, *pc, *pd, *qa, *qb ;
 char xpt[12][MAXLEN] ;

  line[MAXSTR] = '\0' ;
  cleanread(readpt) ;
  if (fgets(line, MAXSTR, fff) != NULL)  {

   t = strlen(line) ; 
   if (t>=(MAXSTR-2)) { 
    readtonl(fff) ; 
    printlong(line) ;
    return -55 ;
   }

   striptrail(line, ' ') ;
   striptrail(line, CNL) ;
   striptrail(line, CTAB) ;

   nsplit = splitupxbuff(line, spt, MAXFF, CTAB, bigb, MAXLEN) ; 
   if (nsplit == 0) return -99 ;
   if (line[0] == '@')  {
     freeup(spt, nsplit) ;
     return -88 ;
   }
   strcpy(readpt -> qname, spt[0]) ;

   for (x=0; x<=10; ++x) { 
    strcpy(xpt[x], spt[x]) ;
   }
   nsplit = 0 ;

   readpt -> ihdr = -1  ; 
   readpt -> ilib = -1  ; 

   ipop = t = whichpop(readpt -> rgstring) ; 
   readpt -> ipop = t = ipop ;


   if (t < 0 )  { 
    return  t ;
   }



   len = strlen(xpt[10]) ;        
   if (len >= MAXLEN) fatalx("bad read\n") ;
   strcpy(readpt -> qual , xpt[10]) ;

   fflush(stdout) ;
   if (quallo != NULL) readpt -> lobasequal = quallo[t] ; 
   if (threshlo != NULL) readpt -> mapthresh =  threshlo[t] ; 

   if (!nopoplist) strcpy(readpt -> pop ,  poplist[t]) ;
   else strcpy(readpt -> pop , "???") ;

   strcpy(readpt -> rname , xpt[2]) ;
   reg = readpt -> rname  ; 
   
   strcpy(readpt -> cigar , xpt[5]) ;

   strcpy(readpt -> seq , xpt[9]) ;
   mkqual(readpt -> iqual, readpt -> qual, len) ;

   readpt -> flag   = atoi(xpt[1]) ;
   readpt -> strand = (readpt -> flag >> 4) & 1 ;
   readpt -> ispair = (readpt -> flag) & 1 ;
   t = readpt -> flag & 0x400 ;  
   if (t>0) readpt -> isdup = YES ;
   pos = readpt -> pos    = atoi(xpt[3]) ;
   readpt -> mapq   = atoi(xpt[4]) ;
   readpt -> isize  = atoi(xpt[8]) ;
   len = readpt -> readlen = strlen(readpt -> seq) ;

   grabreference(reg, pos, len, readpt -> refseq) ;

   readpt -> numposdiff = calcnumposdiff(readpt) ;

   return 1 ;
  }
  return -98;
}

int getreadf(READ *readpt, FILE *fff) 
{
 char line[MAXSTR+1] ;
 char *spt[MAXFF]  ;
 int nsplit ;
 int  t, ipop, x ;
 char *reg, c1, c2, *s1, *s2 ; 
 int pos, len, k ;
 char  pop[MAXLEN] ;
 char bigb[MAXLEN] ; 
 static int nbug = 0 ;
 static int nbad = 0 ;
 int *pa, *pb, *pc, *pd, *qa, *qb ;
 char xpt[12][MAXLEN] ;

  line[MAXSTR] = '\0' ;
  cleanread(readpt) ;
  if (fgets(line, MAXSTR, fff) != NULL)  {

   t = strlen(line) ; 
   if (t>=(MAXSTR-2)) { 
    readtonl(fff) ; 
    printlong(line) ;
    return -55 ;
   }

   striptrail(line, ' ') ;
   striptrail(line, CNL) ;
   striptrail(line, CTAB) ;

   nsplit = splitupxbuff(line, spt, MAXFF, CTAB, bigb, MAXLEN) ; 
   if (nsplit == 0) return -99 ;
   if (line[0] == '@')  {
     freeup(spt, nsplit) ;
     return -88 ;
   }

   strcpy(readpt -> qname, spt[0]) ;
   x = rgindex(spt, nsplit) ;
   if (x < 0) fatalx("bad line (no rg group found):: %s\n", line) ;

   strcpy(readpt -> rgstring, spt[x]) ;
   if (strstr(readpt -> rgstring, "unknown") != NULL)  { 
     return -97 ;
   }

   x = dsindex(spt, nsplit) ;
   if (x<0) {
    readpt -> pmdscore = -98; 
   }
   else { 
    sscanf(spt[x] + 5 , "%lf", &readpt -> pmdscore) ;
   }

   for (x=0; x<=10; ++x) { 
    strcpy(xpt[x], spt[x]) ;
   }
   nsplit = 0 ;

   readpt -> ihdr = -1  ; 
   readpt -> ilib = -1  ; 

   ipop = t = whichpop(readpt -> rgstring) ; 
   readpt -> ipop = t = ipop ;


   if (t < 0 )  { 
    return  t ;
   }



   len = strlen(xpt[10]) ;        
   if (len >= MAXLEN) fatalx("bad read\n") ;
   strcpy(readpt -> qual , xpt[10]) ;

   fflush(stdout) ;
   if (quallo != NULL) readpt -> lobasequal = quallo[t] ; 
   if (threshlo != NULL) readpt -> mapthresh =  threshlo[t] ; 

   if (!nopoplist) strcpy(readpt -> pop ,  poplist[t]) ;
   else strcpy(readpt -> pop , "???") ;

   strcpy(readpt -> rname , xpt[2]) ;
   reg = readpt -> rname  ; 
   
   strcpy(readpt -> cigar , xpt[5]) ;

   strcpy(readpt -> seq , xpt[9]) ;
   mkqual(readpt -> iqual, readpt -> qual, len) ;

   readpt -> flag   = atoi(xpt[1]) ;
   readpt -> strand = (readpt -> flag >> 4) & 1 ;
   readpt -> ispair = (readpt -> flag) & 1 ;
   t = readpt -> flag & 0x400 ;  
   if (t>0) readpt -> isdup = YES ;
   pos = readpt -> pos    = atoi(xpt[3]) ;
   readpt -> mapq   = atoi(xpt[4]) ;
   readpt -> isize  = atoi(xpt[8]) ;
   len = readpt -> readlen = strlen(readpt -> seq) ;

   grabreference(reg, pos, len, readpt -> refseq) ;

   readpt -> numposdiff = calcnumposdiff(readpt) ;

   return 1 ;
  }
  return -98;
}

int calcnumposdiff(READ *readpt) 
{

   char *s1, *s2 ;
   char c1, c2 ;
   int k, t, len ;

   s1 = readpt -> seq ;
   s2 = readpt -> refseq ;
   t = 0 ;

   len = strlen(s1) ;
   for (k=0; k<len; ++k) { 
    c1 = toupper(s1[k]) ;
    c2 = toupper(s2[k]) ;
    if (c1 != c2) { 
     ++t ;
    }
   }
   return t ; 
}

int getread(READ *readpt) 
{

 return getreadf(readpt, stdin) ;

}

void cleanread(READ *readpt) 
{
  readpt -> pop[0] = CNULL ; 
  readpt -> lib[0] = CNULL ; 
  readpt -> samp[0] = CNULL ; 
  readpt -> ipop = readpt -> ilib = readpt -> ihdr = -1 ;
  readpt -> isnean = NO ;
  readpt -> isxw = NO ;
  readpt -> isancient = NO ;

  readpt -> readlen = -1 ;
  readpt -> pos = readpt -> mpos = -9999999 ;
  readpt -> numposdiff = -999 ;
  readpt -> mapthresh = 37  ;
  readpt -> ispair = NO ;
  readpt -> isdup = NO ;
  readpt -> pmdscore = -20000 ;
}
void setancient(READ *readpt) 
{
 char *pop = readpt -> pop ;
 if (strstr(pop, "Nean")) readpt -> isnean = YES ;
 if (strstr(pop, "Altai")) readpt -> isnean = YES ;
 if (strstr(pop, "Vi")) readpt -> isnean = YES ;
 if (strstr(pop, "Mez")) readpt -> isnean = YES ;
 if (strstr(pop, "Feld")) readpt -> isnean = YES ;
 if (strstr(pop, "Sid")) readpt -> isnean = YES ;
 if (strstr(pop, "Denis")) readpt -> isxw = YES ;
 readpt -> isancient = isold(pop) ;
 if (readpt -> isnean) readpt -> isancient = YES ;
 if (readpt -> isxw) readpt -> isancient = YES ;
}

int isold(char *pop) 
{
 if (strstr(pop, "Nean")) return YES ;
 if (strstr(pop, "Altai")) return YES ;
 if (strstr(pop, "Denis")) return YES ;
 if (strstr(pop, "Xwoman")) return YES ;
 return NO ;
}

int isdenisova(char *pop) 
{
 if (strstr(pop, "Denis")) return YES ;
 if (strstr(pop, "Xwoman")) return YES ;
 return NO ;
}

int isneandertal(char *pop) 
{
 if (strstr(pop, "Nean")) return YES ;
 return NO ;
}

void loadqinfo(char *qinfoname, double qnean) 
{
 char line[MAXSTR+1] ;
 char *spt[MAXFF], *sx ;
 char c ;
 int nsplit ;
 int  t, k, s, ind,  a ;
 FILE *fff ;
 char *pname, *pop ;
 double **qq, *q ;
  
  ZALLOC(qtable, numrline, double **) ; 

  for (k=0; k<numrline; ++k) { 
   qq = qtable[k] = initarray_2Ddouble(4, 3, 0.0) ;
   for (a=0; a<4; ++a) { 
    qq[a][0] = 99 ;
    qq[a][1] = 99 ;
    qq[a][2] = 0.0 ;
   }
  }

  if (qinfoname == NULL) fatalx("qinfo compulsory\n") ;
  openit(qinfoname, &fff, "r") ;

  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, MAXFF) ; 
   if (nsplit<1) continue ;
   sx = spt[0] ;
   if (sx[0] == '#') { 
    freeup(spt, nsplit) ;
    continue ;
   }
   
   t = indxstring(rdtable[0], numrline, sx) ; 
   pop = rdtable[3][t] ;
   if (t<0) continue ;
   qq = qtable[t] ; 
   c = spt[2][0] ;
   a = base2num(c) ; 
   q = qq[a] ;
   q[0] = atof(spt[3]) ; // lo
   q[1] = atof(spt[4]) ; // hi
   q[2] = atof(spt[5]) ; // prob
   if (isold(pop)) {
    q[0] = qnean ;
    q[2] = 1.0 ;
   }
   freeup(spt, nsplit) ;
  }

   fclose(fff) ;

}
void printqinfo() 
{
  int a, k ; 
  double *q ;

   for (k=0; k<numrline; ++k) {  
    for (a = 0 ; a < 4; ++a) {
     printf("qtab: %40s ", rdtable[0][k]) ; 
     printf("%c ", num2base(a)) ;
     q = qtable[k][a] ;
     printf("%6.0f ", q[0]) ;
     printf("%6.0f ", q[1]) ;
     printf("%9.3f",  q[2]) ;
     printnl() ;
    }
   }

}

int goodqualp(READ *readpt, int i, THRESH *pthresh)  
{

  int kread, a, q, lo, hi ; 
  double prob, y ;
  double *qt ;
  static long ncall = 0 ;   
  int  skip, t ;
  ++ncall ; 
//  if (ncall > 1000000) ncall = 1000000 ;
  a = base2num(readpt -> seq[i]) ;
  if (a<0) return NO ;
  q = readpt -> iqual[i] ;
  lo = pthresh  -> quallo ;
  hi = pthresh  -> qualhi ;
  if (q<lo) return NO ;
  if (q>hi) return NO ;

 skip = pthresh -> skiplength ;
 if (skip <= 0) return YES ;
 if (i<skip) return NO ;
 t = readpt -> readlen - (i+1) ;
 if (t<skip) return NO ;
 return YES ;
}
int goodqual(READ *readpt, int i)  
{

  int kread, a, q, lo, hi ; 
  double prob, y ;
  double *qt ;
  static long ncall = 0 ;   

  ++ncall ; 
  if (ncall > 1000000) ncall = 1000000 ;
  a = base2num(readpt -> seq[i]) ;
  if (a<0) return NO ;
  if (readpt -> strand == 1) a = 3-a ;
  q = readpt -> iqual[i] ;
  lo = readpt -> lobasequal ;
  if (q<lo) return NO ;
  return YES ;
}

int is_mono(int *aa, int len) 
{
  int  refnum, x ;

  refnum = aa[len-1] ; 

  for (x=0; x< len-1 ; ++x) {
   if (aa[x] < 0) continue ;
   if (aa[x] != refnum) return NO ;
  }
  return YES ; 
}

int is_biallelic(int *aa, int len) 
{
  int  refnum, varnum, x ;

  varnum = -1 ;
  refnum = aa[len-1] ; 

  if (refnum<0) return NO ;

  for (x=0; x< len-1 ; ++x) {
   if (aa[x] < 0) return NO ;
   if (aa[x] == refnum) continue ;
   varnum = aa[x] ;
  }

  for (x=0; x< len-1 ; ++x) {
   if (aa[x] == refnum) continue ;
   if (aa[x] != varnum) return NO ;
  }
  return YES ;

}
void
mkqual(int *iqual, char *qual, int len) 
{
  int k, t ; 

  for (k=0; k<len; ++k) {  
   t = (int) qual[k] ;
   iqual[k] = t -33 ;
   if (ishack) iqual[k] -= 31 ;
   iqual[k] = MAX(iqual[k], 0) ;
  }
}

char *popstring(char **popl, int n) 
{
 static char sss[100] ;
 int k ;

 for (k=0; k<n; ++k) { 
  sss[k] = popl[k][0] ;
 }
 sss[n] = CNULL ;
 return sss  ;
}

int isloshlass(char *pop) 
{
 int t ;
 t = strcmp("Loshbour", pop) ;
 if (t==0) return YES ;
 t = strcmp("Lassithi", pop) ;
 if (t==0) return YES ;
 if (strstr(pop, "Skog")) return YES ;
 return NO ;
}

int getfastalist(char **poplist, int npops, char *dbfile, char **iublist, int *iubpops) 

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
    iublist[nx] = strdup(sx) ;
    if (iubpops != NULL) { 
      iubpops[nx] = t ; 
    }
//  printf("zzziub: %d %d %s\n", nx, t, sx) ;
    ++nx ;
    freeup(spt, nsplit) ;
  }

   fclose(fff) ;
   return nx ;

}

char *getfaiub(int k) 
{

  int x ;

 for (x=0; x<numiub; ++x) { 
  if (iubpops[x] == k)  return iublist[x] ; 
 }
 for (x=0; x<numiub; ++x) { 
   printf("badiub: %d :: %s\n", iubpops[x], iublist[x]) ;
 }

 fatalx("(getfaiub) pop number %d not found\n", k) ;

}
int mkfullpoplist(char **fullpoplist, int *f2pop, char **rawpoplist, int npops) 
{
  int k, nsplit, nfpops = 0  ;
  char *spt[MAXFF], *sx ;

  for (k=0; k<npops; ++k) { 
   nsplit = splitupx(rawpoplist[k], spt, MAXFF, ':') ; 
   copystrings(spt, fullpoplist+nfpops, nsplit) ;
   if (f2pop != NULL) ivclear(f2pop+nfpops, k, nsplit) ;
   nfpops += nsplit ;
   freeup(spt, nsplit) ;
  }
     
  return nfpops ;

}

int getrg(char **rawpoplist, int npops, char *dbfile, char **rglist, int *rgpops) 
{
 char line[MAXSTR+1] ;
 char *spt[MAXFF], *sx ;
 char c ;
 int nsplit ;
 int  t, k, s, nx = 0, na, nb  ;
 FILE *fff ;
 char *scolon ; 
 char *sss ;
 char **trgstring ;
 char **fullpoplist ; 
 int  *f2pop, maxfpops=500, nfpops, popindex ;
  
  if (npops==0) { 
   nopoplist = YES ;
   return 0 ;
  }

  ZALLOC(fullpoplist, maxfpops, char *) ;
  ZALLOC(f2pop, maxfpops, int) ;
  nfpops = mkfullpoplist(fullpoplist, f2pop, rawpoplist, npops) ;

  openit(dbfile, &fff, "r") ;

  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitupx(line, spt, MAXFF, CTAB) ; 
   if (nsplit<1) continue ;
   sx = spt[0] ;
   if ((sx[0] == '#') || (nsplit == 1)) { 
    freeup(spt, nsplit) ;
    continue ;
   }
   if (nsplit < 6) { 
    for (k=0; k<nsplit; ++k) { 
     printf("rg bad db: %d %s\n", k, spt[k]) ;
    }
    fatalx("bad rg line:: %s\n", line) ;
   }
   striptrail(sx, ' ') ;
   t = indxstring(fullpoplist, nfpops, sx) ; 
   if (t<0) { 
    freeup(spt, nsplit) ; 
    continue ;
   }
   popindex = f2pop[t] ;
   for (k=0; k<=5; ++k) { 
    if (!verbose) break ; 
    printf("zzrg: %d %s\n", k, spt[k]) ;
   }
   sx = strdup(spt[5]) ;
   freeup(spt, nsplit) ;

   if (strstr(sx, "[FILE:") != NULL) { 
    substring(&sx, "[FILE:", "") ;
    substring(&sx, "]", "") ;
    substring(&sx, "\n", "") ;
 // printf ("file ::%s::\n", sx) ;
    na = numlines(sx) ; 
    ZALLOC(trgstring, na, char *) ;
    nb = getss(trgstring, sx) ;
    freestring(&sx) ;

    for (k=0; k<nb; ++k) { 
     sx = trgstring[k] ;
     rglist[nx] = strdup(sx) ;
     rgpops[nx] = t ; 
     ++nx ;
    }

    freeup(trgstring, na) ;
    continue ; 
   }
   
   subcolon(sx) ; // : -> space
   nsplit = splitup(sx, spt, MAXFF) ; 
   freestring(&sx) ;
   for (k=0; k<nsplit; ++k) { 
    sx = spt[k] ;
    rglist[nx] = strdup(sx) ;
    rgpops[nx] = popindex ; 
    ++nx ;
   }
   freeup(spt, nsplit) ;
  }
   freeup(fullpoplist, maxfpops) ;
   free(f2pop)  ;

   fclose(fff) ;
   return nx ;

}


int countc(char *uu, int nuu, char ref, char  var) 
{

 int t, j, k, x ; // uu is iub decoded genotypes
 int snpvar, snpref ;

 snpvar = base2num(var) ;
 snpref = base2num(ref) ;
  
 return countg(uu, nuu, snpref, snpvar) ;

}



int countg(char *uu, int nuu, int snpref, int snpvar) 
{

 int t, j, k, x ; // uu is iub decoded genotypes
  
 if (snpref<0) return -1 ;
 if (nuu == 0) return -2 ;
 t = 0 ;
 for (k=0; k<2 ; ++k) { 
  j = MIN(k, nuu-1) ;
  x = base2num(uu[j]) ;
  if ((x != snpvar) && (x != snpref)) return -1 ;
  if (x == snpref) ++t ;
 }

 return t ;

}

void printread(READ *readpt) 
{
 int i, k ;
 double y ;
 char sss[500] ;
 extern char * regstring ;

 printf("qname: %s", readpt -> qname) ;
 printf("%10d %10d %4d\n", readpt -> pos, readpt -> isize, readpt -> numposdiff) ;
 printf("pop: %s  %d\n", readpt -> pop, readpt -> ipop) ;
 printf("rname: %s\n", readpt -> rname) ;
 printf("pos: %12d\n", readpt -> pos) ;
 printf("strand: %d\n", readpt -> strand) ;
 printf("duplicate: %d\n", readpt -> isdup) ;
 printf("cigar: %s\n", readpt -> cigar) ;
 printf("mapq: %d\n", readpt -> mapq) ;
 printf("flag: %0x\n",  readpt -> flag) ;
 printf("isnean: %d\n", readpt -> isnean) ;
 printf("isxw: %d\n", readpt -> isxw) ;
 printf("samp: %s\n", readpt -> samp) ;
 printf("len:  %d\n",  readpt -> readlen) ;
 printf("seq:  %s\n",  readpt -> seq) ;
 printf("ref:  %s\n",  readpt -> refseq) ;
 y = readpt -> pmdscore ;
 if (y > -10000) printf("pmdscore: %9.3f\n", y) ; 

 if (regstring != NULL) {
  cclear(sss, CNULL, 500) ;
  strncpy(sss, regstring + (readpt -> pos) - 1, MIN(readpt -> readlen, 499)) ; 

  for (k=0; k<strlen(sss); ++k) { 
   if (sss[k] == 'N') sss[k] = ' ' ;
   if (sss[k] == 'X') sss[k] = ' ' ;
  }

  printf("reg:  %s\n",  sss) ;
 }
 

/**
 printf("diff: ") ;
 for (i=0; i<readpt -> readlen; ++i) { 
  printf("%1d", readpt -> diff[i]) ;
 }
*/
 for (i=0; i<readpt -> readlen; ++i) { 
  if ((i%10) == 0) printnl()  ;
  printf("qual %d ", readpt -> iqual[i]) ;
 }
 printnl() ;
}




