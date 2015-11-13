// Written by Nick Patterson
// Revised by Mengyao Zhao on 2015-11-13

#include <sys/wait.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <libgen.h>
#include <nicksam.h>
#include <getpars.h>  
#include "bam.h"
#include "faidx.h"
#include "globals.h" 
#include "popsubs.h"
#include "mcio.h"
#include <nicklib.h>

typedef struct { 
 int depthlo ; 
 int depthhi ; 
 int cat   ; 
 double  ybestw ; 
 double  ybestx ; 
 double bestfrac ;
 double coverage ; 
 double hrefdiverge ; 
 long hrefhit ; 
 long hrefmiss ;
}  BANS ; 

//int debug = NO ;

BANS bans[2000] ;

char *regname = NULL ; 
 char *regnames[30], ss[256] ;
//char *poplistname = NULL ; 
char *iubfile = "/home/mz128/group/Cteam/filter/hetfa.dblist" ;
char *iubmaskfile = "/home/mz128/group/Cteam/filter/mask.dblist" ;
char *parname = NULL ;
char *maskname = NULL ;
char *hetfaname = NULL ; 

char *sampname = NULL ;
//char *dbase = NULL  ;
char *fixeddbase = NULL ;
//char *lov= NULL, *hiv = NULL ;
//char *mapstring = "map90" ;
//char *cnvname = NULL ;

char *chimp = NULL;
char *href = NULL;
//char *heng75 = NULL;

//int hipos, lopos ;
int nfregs ; 
int maxhist = 100 ; 
 int maxqval = 99 ;  // clip depth MQ and MQ0
char *refname, *sampfname ;
char gender = 'U' ; 
int goodmq = 37 ;
int totlen = 0, totvalid = 0 ;  

#define VERSION  "200"      
// first pass at a version for release just containing what is needed
void readcommands(int argc, char **argv) ;
int getfastalist(char **poplist, int npops, char *dbfile, char **iublist, int *iubpops) ;
char *myfai_fetch(faidx_t *fai, char *reg, int  *plen) ;
int minfalen(FATYPE **fainfo, int n)  ;
void  mkhist(char **regnames, long *hh, long *mm, int *blist, 
 int rlo, int rhi, char **fqnames, int nfqnames, int *lovals, int *hivals)   ;
void  mkmask(FILE *fff, char **regnames, BANS **bbans, int *filtval, int rlo, int rhi, char **fqnames, int nfqnames, int *lovals, int *hivals)  ;
void easymask(FILE *ffmask, char **regnames, int rlo, int rhi, char maskval)  ; 
void setfpars(long **hh, BANS **bbans)  ;  
void printfilter(BANS **bans)  ;
void writefaf(FILE *fff, char *regname, char *fastring, int falen) ;

int readfa(char **falist, char **fasta, int *flen, int n) ;
int getdbname(char *dbase, char *name, char **pfqname) ;
int matchit(char ch, char ct)  ;
int hk(char a, char b, char c) ; 
double missfrac(int a, int b) ;
double ymissfrac(double a, double b) ;

void countit(char **fqnames, int rlo, int rhi, long **sschimp, long **sshref)  ;
void printhist(long **ss) ;

static int usage() 
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   filtstats -p <parameter file> [options]\n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-V	Print more information while the program is running.\n");
	fprintf(stderr, "\t-v	Show version information.\n");
	fprintf(stderr, "\t-? 	Show the instruction. (For detailed instruction, please see the document here: https://github.com/mengyao/cTools)\n\n");
	return 1;
}

int main(int argc, char **argv)
{

//Declare a million variables and initialize
 char ca, cb, c1, c2   ;
 int cc;
 int ccc[5];
 int *flen, *fqlen ;
 int i, j, a, b, c, d, w, x, aa, bb ; 
 int t, k, tt ;
 char **fqdata,  **fqnames ; 
   
 int *ismatch, *isrmatch, *isvalid, *hkode, len, fixlen ;
 int lovals[10], hivals[10] ; 
 int nfqnames = 5 ; 
 int  *filtval ;
 FILE *ffmask ;
 char *hengfilt, *cnvfilt ; 
 long **sschimp, **sshref ;

 BANS  *bbans[10], *banspt ;

// step 1 initialize regnames ;
 for (a=1; a<= 22; ++a) { 
  sprintf(ss, "%d", a) ; 
  regnames[a] = strdup(ss) ;
 }
 regnames[23] = strdup("X") ;
 regnames[24] = strdup("Y") ;
 regnames[25] = strdup("MT") ;

 //Read parameter file and initialize
 readcommands(argc, argv) ;
 
	if (hetfaname == NULL) getdbname(iubfile, sampname, &hetfaname) ;

	if (maskname == NULL) getdbname(iubmaskfile, sampname, &maskname) ;


 nfqnames = 4 ;
 ZALLOC(fqnames, nfqnames, char *) ; 
 ZALLOC(flen, nfqnames, int) ; 
 
 
	if (fixeddbase != NULL) {
  		getdbname(fixeddbase, "Chimp", &fqnames[0]) ; //refname = fqnames[0] ;
  		getdbname(fixeddbase, "Href", &fqnames[1]) ; //refname = fqnames[0] ;
 	} else {
 		fqnames[0] = chimp;
 		fqnames[1] = href;
 	}

 fqnames[2] = sampfname = hetfaname ; 
 fqnames[3] = maskname  ; 

 ZALLOC(fqdata, nfqnames, char *) ;

 sschimp = initarray_2Dlong(11, 2, 0) ;
 sshref = initarray_2Dlong(11, 2, 0) ;

 countit(fqnames, 1, 22, sschimp, sshref) ;
 printf("chimp_divergence\n") ; printhist(sschimp) ; 
 printf("Href_divergence\n") ; printhist(sshref) ; 
 printf("## end of filtstats\n") ;
 return 0 ;
}

int vmask(char c1)    
{
  int t ; 
  t = (int) c1 - (int) '0' ; 
  if (t<0) t = 10 ;
  if (t>9) t = 10 ;

  return t ;

}
void countit(char **fqnames, int rlo, int rhi, long **sschimp, long **sshref) 
{
  int *fqlen; 
  char **fqdata ; 
  int nfqnames = 4, a, rnum, i ;
  int lchimp, lhref, t ;
  char cref, chetfa; 
  int kmask ;


  printf("chromosomes: %d-%d\n", rlo, rhi) ;
  ZALLOC(fqlen, nfqnames, int) ;
 for (rnum=rlo; rnum <=rhi; ++rnum) {
  ZALLOC(fqdata, nfqnames, char *) ;
  regname = regnames[rnum] ;
  readfa(fqnames, fqdata, fqlen, nfqnames) ;
  for (a=0; a<nfqnames; ++a) {
   printf("fq: %3d %s %12d\n", a, fqnames[a], fqlen[a]) ;
  }
  t = MIN(fqlen[2], fqlen[3]) ;
  lchimp = MIN(fqlen[0],  t) ; 
  lhref = MIN(fqlen[1],  t) ; 
  for (i=0; i<lchimp ; ++i) { 
   cref    = toupper(fqdata[0][i]) ; 
   chetfa  = toupper(fqdata[2][i]) ; 
   kmask =  vmask(fqdata[3][i]) ; 
   t = matchit(cref, chetfa) ; 
   if (t<0) continue ;
   sschimp[kmask][0] += t ;
   sschimp[kmask][1] += 2-t ;
  }
  for (i=0; i<lchimp ; ++i) { 
   cref = toupper(fqdata[1][i]) ; 
   chetfa = toupper(fqdata[2][i]) ; 
   kmask =  vmask(fqdata[3][i]) ; 
   t = matchit(cref, chetfa) ; 
   if (t<0) continue ;
   sshref[kmask][0] += t ;
   sshref[kmask][1] += 2-t ;
  }

  freeup(fqdata, nfqnames) ;
 }

 free(fqlen) ;
}

void printhist(long **ss)
{
 int k ;
 long a, b, total, tt[2] ; 

   for (k=0; k<11; ++k) {
    total += ss[k][0] ;
    total += ss[k][1] ;
   }
   k = 10 ; 
   printf("%8s ", "no mask:") ;  
   a = ss[k][0] ; b=ss[k][1] ; 
   tt[0] = a ; tt[1] = b ; 
   printf("%10ld ", a + b ) ; 
   printf("%12.6f", ymissfrac(a,b)) ;
   printnl() ;
   for (k=0; k<=9; ++k) { 
    printf("%8d ", k) ;            
    a = ss[k][0] ; b=ss[k][1] ; 
    tt[0] += a ; tt[1] += b ; 
    printf("%10ld ", a + b ) ; 
    printf("%12.6f", ymissfrac(a,b)) ;
    printnl() ;
   }
   printf("%8s ", "Total:") ;  
   a = tt[0]  ; b = tt[1] ; 
   printf("%10ld ", a + b ) ; 
   printf("%12.6f", ymissfrac(a,b)) ;
   printnl() ;
   printnl() ;

}

double ymissfrac(double a, double b)
{
 double t = a + b ;

 if (t<=0.1) return  0 ;
  return  b / t ;

}

void printfilter(BANS **bbans) 
{

  BANS *lastpt, *banspt ; 
  int k ; 
  double yy, yw, yx, ytot ;

  lastpt = NULL ;
  for (k=0; k<=9; ++k) { 
   banspt = bbans[k] ;
   printf("filtercat: %3d", k) ; 
   printf(" %3d %3d ", banspt -> depthlo, banspt -> depthhi) ; 
   yy = banspt -> coverage ; 
   printf("  %9.3f", yy) ;
   printf("  %12.6f", banspt -> bestfrac) ;
   if (lastpt != NULL) { 
     yw = lastpt -> ybestw - banspt -> ybestw ;
     yx = lastpt -> ybestx - banspt -> ybestx ;
     if ((yw + yx) > 0.1)
      printf("   %12.6f", ymissfrac(yw, yx)) ;
     else 
      printf("   %12s", "-" ) ;
   }
   else {
     printf("   %12s", "-" ) ;
   }
   printf(" %10ld", banspt -> hrefhit) ;
   printf(" %10ld", banspt -> hrefmiss) ;
   printf("   %12.6f", banspt -> hrefdiverge) ;
   lastpt = banspt ;
   printnl() ;
  }

}

void setfpars(long **hh, BANS **bbans)  
{
 
// original code .../bestwindow.c 

  long *ttt, *hdata, *mdata ; 
  int j, k, lo, hi, coverpc, a, b ;  // left hand window must be this large 
  double ytot, yw, yx, yy ; 
  long *hbases, *mbases ;
  int u, v,  besta, bestb  ; 
  double bestfrac ; 
  double ybestw, ybestx, y0, y1, xsc, oldbestw, ymaxw, ymaxx, ymax, ysize ;  
  BANS *banspt ; 
  int numbans = 0 ;

  ZALLOC(ttt,  maxhist+1, long) ;
  ZALLOC(hdata, maxhist+1, long) ;
  ZALLOC(mdata, maxhist+1, long) ;
  ZALLOC(hbases,  maxhist+1, long) ;
  ZALLOC(mbases,  maxhist+1, long) ;

   ytot = 0.0 ;
   for (k=0; k<=maxhist; ++k) { 
    ttt[k] = hh[k][0] ; 
    ytot += ttt[k] ;
    hdata[k] = hh[k][1] ; 
    mdata[k] = hh[k][2] ; 
   }

   oldbestw  = -1.0 ; 
   lo =  9999 ; hi = 0 ;

  for (coverpc  = 200; coverpc < 1000; ++coverpc) { 
   bestfrac = 10.0 ;
   besta = bestb = -99 ;

   for (a=1; a<=maxhist; ++a) { 

     u = hdata[a] ; 
     v = mdata[a] ;
   
    if ((u==0) && (v==0)) continue ;

    lvzero(hbases, maxhist+1) ;
    lvzero(mbases, maxhist+1) ;

    if (a>lo) continue ;

    for (b=a; b<=maxhist; ++b) { 

     u=hdata[b] ; 
     v=mdata[b] ;

      yw = hbases[b] = hbases[b-1] + u ; 
      yx = mbases[b] = mbases[b-1] + v ; 

     if ((u==0) && (v==0)) continue ;
     if (b<hi) continue ;
      y1 = (yw+yx)  / ytot ;
      y1 *= 1000.0 ;
      if (y1 < (double) coverpc) continue ;
      xsc = ymissfrac(yw, yx) ; 
      if (xsc < bestfrac) {  
        bestfrac = xsc ; 
        besta = a ; bestb = b ; 
        ybestw = yw ; ybestx = yx ;
      }
     }
  }
// if (bestfrac>1.0) continue ;
   if (besta == -99) continue ;
   if (oldbestw == ybestw) continue ;
   lo = MIN(lo, besta) ;
   hi = MAX(hi, bestb) ;
   oldbestw = ybestw ;
   y1 = (ybestw+ybestx) / ytot ;
   printf("bestlim: %3d %4d ", 0, coverpc) ;  
   printf("%3d %3d ", besta, bestb) ; 
   printf("  %12.0f %12.0f ", ybestw, ybestx) ; 
   printf("%12.6f", bestfrac) ; 
   printf("   %12.6f", y1) ; 
   printnl() ;

   banspt = &bans[numbans] ; 
   banspt  -> depthlo = besta ;
   banspt  -> depthhi = bestb ;
   banspt  -> ybestw = ybestw ;
   banspt  -> ybestx = ybestx ;
   banspt  -> bestfrac = bestfrac ;
   banspt  -> coverage = y1 ;
   banspt  -> hrefdiverge = 0.0 ; 
   banspt  -> hrefhit = banspt -> hrefmiss = 0 ;
   ++numbans ;
 }

 ymaxw = longsum(hdata, maxhist+1) ;
 ymaxx = longsum(mdata, maxhist+1) ;

// add code for no filter 
 ybestw = ymaxw ;
 ybestx = ymaxx ;
 bestfrac = ymissfrac(ybestw, ybestx) ;
 besta = 0 ; bestb = 99 ;
 y1 = (ybestw+ybestx) / ytot ;

  if (ybestw > oldbestw) {
   printf("bestlim: %3d %4d ", 0, 1000) ;  
   printf("%3d %3d ", besta, bestb) ; 
   printf("  %12.0f %12.0f ", ybestw, ybestx) ; 
   printf("%12.6f", bestfrac) ; 
   printf("   %12.6f", y1) ; 
   printnl() ;
   banspt = &bans[numbans] ; 
   banspt  -> depthlo = besta ;
   banspt  -> depthhi = bestb ;
   banspt  -> ybestw = ybestw ;
   banspt  -> ybestx = ybestx ;
   banspt  -> bestfrac = bestfrac ;
   banspt  -> coverage = y1 ;
   banspt  -> cat = -1 ;
   banspt  -> hrefdiverge = 0.0 ; 
   banspt  -> hrefhit = banspt -> hrefmiss = 0 ;
   ++numbans ;
  }

  ymax = ybestw + ybestx ; 

  y0 = ymax/ytot ; y1 = 0.5 ;   // we want max filter at 50% coverage
  if (y0 < y1)  {  // no filtering 
   for (k=0; k<=9; ++k) { 
    bbans[k] = banspt ;
   }
  }

  bbans[0] = banspt ;
  for (k=1; k<=9; ++k) { 
    yy = ((double) (9-k)  * y0 + (double) k * y1 ) / 9.0 ;
     for (j=numbans-1; j>=0; --j) { 
      banspt = &bans[j] ; 
      ysize = banspt -> coverage ;
      if (ysize >= yy) bbans[k] = banspt ;
     }
  }

 
 return ;

}
int setvalid(int *valid, int len, char **fqdata, int nfqnames) 
{
 int i, pos, x ;
 char cref, csamp, ca, cb, cc ;
 int failtype[4] ;

 ivzero(valid, len) ;
 ivzero(failtype, 4) ;
  for (i=0; i<len; ++i) { 
   pos = i+1 ; 
   cref =  toupper(fqdata[0][i]) ;
   csamp = toupper(fqdata[3][i]) ; 

   ca = fqdata[1][i] ; // heng universal filter
   cb = fqdata[2][i] ;  // cnv filter

   x = ranmod(100*1000) ; 
   if (x==0) { 
    printf("zzval %d ", pos);
    printf("%c %c %c %c\n", cref, csamp, ca, cb) ;
   }
  
   x = 99 ;
   if (ca != '1') x = 3 ;    
   if (islower(cb)) x = 2 ;    
   if (base2num(cb) < 0) x = 2 ;    
   if (!isiub2(csamp)) x = 1 ;    
   if (base2num(cref) < 0) x = 0  ;
   if (x<99) { 
    ++failtype[x] ;
    continue ;
   }
   valid[i] = 1 ;
 }
 printf("failtype: %d   ", len) ;
 printimatl(failtype, 1, 4) ;
}

double missfrac(int a, int b) 
{
 int t = a + b ;

 if (t==0) return 0 ;
  return (double) b / (double) t ;

}

int matchit(char ch, char ct) 
{

 char aa[5] , c1, c2 ; 
 int n, t ;

 c1 = toupper(ch) ; 
 c2 = toupper(ct) ; 

 if (verbose) printf("zzmat0: %c %c\n", c1, c2) ; 
 if (!isiub2(c1)) return -1 ;
 if (!isiub2(c2)) return -1 ;

 n = iubdekode(aa, c2) ;  

 if ((n==1) && (c1==c2)) return 2 ;
 if ((n==1) && (c1!=c2)) return 0 ;

// n = 2 
  t = 0 ; 
  if (c1==aa[0]) ++t ;
  if (c1==aa[1]) ++t ;

  if (verbose) printf("zzmat: %d\n", t) ;
  return t ;

}

int getdbname(char *dbase, char *name, char **pfqname) 
{
 char ***names ;  
 int n, k, t, i ; 

 n = numlines(dbase) ;

 ZALLOC(names, 3, char **) ;

 for (i=0; i<=2; ++i) {
  for (k=0; k<n; ++k) { 
   ZALLOC(names[i], n, char *) ;
  }
 }

 n = getnames(&names, n, 3, dbase) ;
 t = indxstring(names[0], n, name) ; 
 if (t<0) fatalx("%s not found in %s\n", name, dbase) ;
 *pfqname = strdup(names[2][t]) ;
 
 for (i=0; i<=2; ++i) { 
  freeup(names[i], n) ;
  free(names[i]) ; 
 }
 free(names) ;
 
 return 1 ; 
}
void mulmat2D(double *x, double **mat, double *w, int a, int b) 
{
  int i, j ;
  double *ww ;

  ZALLOC(ww, a, double) ;
  vzero(ww, a) ; 

  for (i=0; i<a; ++i) { 
   for (j=0; j<b; ++j) { 
    ww[i] += mat[i][j]*w[j] ;
   }
  }
  copyarr(ww, x, a) ;
  free(ww) ;
}
void mulmat2Dt(double *x, double *w, double **mat, int a, int b) 
{
  int i, j ;
  double *ww ;

  ZALLOC(ww, b, double) ;
  vzero(ww, b) ; 

  for (i=0; i<a; ++i) { 
   for (j=0; j<b; ++j) { 
    ww[j] += w[i] * mat[i][j] ;
   }
  }
  copyarr(ww, x, b) ;
  free(ww) ;
}
int readfa(char **falist, char **fasta, int *flen, int n) 
{

 faidx_t *fai ;
 int k, len ;
 
 for (k=0; k<n; ++k) {
  fai = fai_load(falist[k]) ;
  fasta[k] = myfai_fetch(fai, regname, &len) ;
  flen[k] = len ;
 }
}
void readcommands(int argc, char **argv) 

{
  int i, t;
  phandle *ph ;
  char str[512]  ;
  int n, kode ;
  int pops[2] ;
  char *ss = NULL ;

	if (argc < 2) exit(usage());

  while ((i = getopt (argc, argv, "p:vV")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
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

   if (parname == NULL) return ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

 //  getstring(ph, "regname:", &regname) ;
  // getstring(ph, "poplistname:", &poplistname) ;
   getstring(ph, "iubfile:", &iubfile) ;
   getstring(ph, "iubmaskfile:", &iubmaskfile) ;
  // getstring(ph, "mapstring:", &mapstring) ;
   //getstring(ph, "cnv:", &cnvname) ;
   getstring(ph, "gender:", &ss) ;
   if  (ss != NULL)  { 
    gender = ss[0] ;
    freestring(&ss) ;
   }
	//t = 0 ; getint(ph, "lopos:", &lopos) ; lopos = MAX(lopos, t) ;
  // t = BIGINT ; getint(ph, "hipos:", &hipos) ; hipos = MIN(hipos, t) ;

   //getstring(ph, "dbase:", &dbase) ;
   getstring(ph, "hetfa:", &hetfaname) ;
   getstring(ph, "fixeddbase:", &fixeddbase) ;
   getstring(ph, "sampname:", &sampname) ;
  // getstring(ph, "popname:", &sampname) ;
   getstring(ph, "maskname:", &maskname) ;

 //  getstring(ph, "lovals:", &lov) ;
 //  getstring(ph, "hivals:", &hiv) ;
//   getint(ph, "debug:", &debug) ;

	getstring(ph, "href:", &href) ;
    getstring(ph, "chimp:", &chimp) ;
   // getstring(ph, "heng75:", &heng75) ;

   printf("### THE INPUT PARAMETERS\n");
   printf("##PARAMETER NAME: VALUE\n");
   writepars(ph);
   closepars(ph) ;
   fflush(stdout) ;

}

void setlimv(char *vvv, int *vals, int n) 
{
  char *w1 ; 
  char *spt[MAXFF], *spt2[MAXFF], *sx ;
  int nsplit, n2, k ;
  int a, b ;



  if (vvv==NULL) return ;
  w1 = strdup(vvv) ;
  substring(&w1, "::", " ") ;       
  nsplit = splitup(w1, spt, MAXFF) ; 

  for (k=0; k<nsplit; ++k) { 
   sx = spt[k] ;
   n2 = splitupx(sx, spt2, 2, ':') ;
   if (n2!=2) fatalx("bad lim: %s\n", vvv) ;

   a = atoi(spt2[0]) ; 
   if ((a<0) || (a>n)) fatalx("bad limv bounds: %s\n", vvv) ;
   b = atoi(spt2[1]) ; 

   vals[a] = b ;
   freeup(spt2, n2) ;
  }

  freeup(spt, nsplit) ;
  freestring(&w1) ;


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
 fapt -> lopos = fapt -> hipos = -1 ;
 fapt -> len = 0 ;
 fapt -> mlen = 0 ;
 fapt -> rlen = 0 ;

}

int hk(char a, char b, char c)  
{
  int uu[3], t ;
  char aa ;

  ivzero(uu, 3) ;
  aa = a ; if (!isupper(a)) aa = '?' ;
  if (base2num(aa) >=0) uu[0] = 1 ;
  if (b == '3') uu[1] = 1 ;
  aa = c ; if (!isupper(c)) aa = '?' ;
  if (base2num(aa) >=0) uu[2] = 1 ;
  t = 4 * uu[0] + 2 * uu[1] + uu[2]  ;
  return t ;
}

int istransition(char iubc) 
{
  char aa[5]  ;
  int n, abxkode, a, b ; 

  n = iubdekode(aa, toupper(iubc)) ;  
  if (n!=2) return -1 ;
  a = base2num(aa[0]) ; 
  b = base2num(aa[1]) ;
  abxkode = abx(a, b) ;
  return  abxok(abxkode, 3) ;

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

int abxok(int abx, int abxmode) 

{ 

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

void writefaf(FILE *fff, char *regname, char *fastring, int falen) 
{
  char bb[51] ;
  int k ; 

  bb[50] = CNULL ;

  fprintf(fff, ">%s\n", regname) ;
  for (k=0; k<falen; k+=50) { 
   cclear(bb, CNULL, 50) ;
   strncpy(bb, fastring+k, 50) ;
   fprintf(fff, "%s\n", bb) ;
  }
}
