#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>

#include <nicklib.h>  
#include "pop5subs.h"
#include "nicksam.h"

#define MAXSTR 512 
#define MAXFF 50

extern int verbose ;

static char ***bonetable ;
static int numbline = 0 ; 


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
}

void loadbonetable(char *bfname)  
{
   int t, i ;
   t = numlines(bfname) ;
   ZALLOC(bonetable, 3, char **) ;
   for (i=0; i<3; i++) { 
    ZALLOC(bonetable[i], t, char *) ;
   }
   numbline = getxstr(bonetable, t, 3, bfname) ;
   
   if (verbose == NO) return ;
   for (i=0; i< numbline; i++) { 

    printf("%30s ", bonetable[0][i]) ;
    printf("%30s ", bonetable[1][i]) ;
    printf("%30s ", bonetable[2][i]) ;
    printnl() ;
   }
}

int getbone(char *library, char *bone, char *hdr, char *bonefile) 
{
  
 char *spt[MAXFF] ;
 char *ss[4] ;
 char sss[MAXSTR] ; 
 int nsplit, k ;
   
 if (numbline == 0) loadbonetable(bonefile) ;
 nsplit = split1(hdr, spt,  ':') ;

 ss[0] = spt[0] ;
 ss[1] = ss[3] = ":" ;
 nsplit = split1(spt[1], spt+2, ':')  + 2 ;
//  printf("zz %s %s %s\n", spt[0], spt[1], spt[2]) ;
 ss[2] = spt[2] ;
 
 catx(sss, ss, 4) ;
 if (verbose) printf("hdr: %s\n", sss) ;
 
 k = indxsub(hdr, bonetable[0], numbline) ;

 freeup(spt, nsplit);

 if (k<0) return -1 ;

 strcpy(library, bonetable[1][k]) ;
 strcpy(bone, bonetable[2][k]) ;
 
 return k ;


}

int 
getxstr(char ***xx, int maxrow, int numcol, char *fname)
{

  char line[MAXSTR] ;
  char *spt[MAXFF] ;
  char *sx ;
  int nsplit, i, j, num=0, maxff ;
  FILE *fff ;
  int nbad = 0 ; 

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

/**
char *binary_string(int a, int len) 
{
  static char ss[100] ;
  int t = a, k, i ;
  char *binary = "01" ;

  ss[len] = CNULL ;
  for (i=0; i<len; i++) { 
   k = t % 2 ;
   ss[len-i-1] = binary[k] ;
   t = t/2 ;
  }
  return ss ;
// fragile
}
*/

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

char *getfaname(char *sname) 
{

static char **fanames = NULL ;
static char **snames = NULL, sss[256] ;

char *dir = "/nfs/eva/nickp/data/jun09" ;
char *faname ;
int t ;

 if (fanames == NULL) {
  ZALLOC(fanames, 3, char *) ;
  ZALLOC(snames, 3, char *) ;
  fanames[0] = strdup("NA18507-half-chimp.fa") ;
  fanames[1] = strdup("SJK-chimp.fa") ;
  fanames[2] = strdup("panTro2.fa") ;

  snames[0]  = strdup("Yoruba") ;
  snames[1]  = strdup("Korean") ;
  snames[2]  = strdup("Chimp") ;

 }

 t = indxstring(snames, 3, sname)  ;
 faname = fanames[t] ;
 makedfn (dir, faname, sss, 256) ;
 printf("## %s: %s\n", sname, sss) ;

 return sss ;
}
