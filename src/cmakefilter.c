// Revised by Mengyao Zhao on 2015-11-09

#include <libgen.h>
#include <nicksam.h>
#include "bam.h"
#include "faidx.h"
#include "globals.h" 
#include "popsubs.h"
#include "mcio.h"
#include <nicklib.h>
#include "opth.h" 

int debug = NO ;


char *regname = NULL ; 
//char *poplistname = NULL ; 
//char *iubfile = "/home/np29/biology/neander/nickdir//xwdir/may12src/altai/myfasta.dblist" ;
char *parname = NULL ;
char *maskname = NULL ;
char *vcfdir = "." ;
char *vcfsuffix = NULL ; 
char *hetfaname = NULL ; 

char *vcfname;// *ref ;
char *sampname = NULL ;
char *dbase = NULL  ;
char *fixeddbase = NULL ;
char *lov= NULL, *hiv = NULL ;
//char *mapstring = "map90" ;
char *cnvname = NULL ;

char *chimp = NULL;
char *href = NULL;
char *heng75 = NULL;

int hipos, lopos ;
int nfregs ; 
int maxhist = 100 ; 
 int maxqval = 99 ;  // clip depth MQ and MQ0
char *refname, *sampfname ;
char gender = 'U' ; 
int goodmq = 37 ;
long totlen = 0, totvalid = 0 ;  

#define VERSION  "300"      
// first pass at a version for release just containing what is needed
void readcommands(int argc, char **argv) ;
int getfastalist(char **poplist, int npops, char *dbfile, char **iublist, int *iubpops) ;
char *myfai_fetch(faidx_t *fai, char *reg, int  *plen) ;
int minfalen(FATYPE **fainfo, int n)  ;
int  mkhist(char **regnames, long *hh, long *mm, int *blist, 
 int rlo, int rhi, char **fqnames, int nfqnames, int *lovals, int *hivals)   ;
void easymask(FILE *ffmask, char **regnames, int rlo, int rhi, char maskval)  ; 
void writefaf(FILE *fff, char *regname, char *fastring, int falen) ;
void  mkmask(FILE *fff, char **regnames, FENTRY **fans,  int rlo, int rhi, char **fqnames, int nfqnames) ;

int readfa(char **falist, char **fasta, int *flen, int n) ;
int readvcf(char *falist, int **fasta, int len, int *flen) ;
int getdbname(char *dbase, char *name, char **pfqname) ;
int matchit(char ch, char ct)  ;
int hk(char a, char b, char c) ; 
double missfrac(int a, int b) ;
double ymissfrac(double a, double b) ;
void printhist(long **ss) ;

void setlimv(char *vvv, int *vals, int n) ;
void calcfiltermask(FILE *ffmask, char **regnames, int rlo, int rhi, 
  int rloout, int rhiout,char **fqnames, int nfqnames, int *lovals, int *hivals)  ;

static int usage() 
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   cmakefilter -p <parameter file> [options]\n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-V	Print more information while the program is running.\n");
	fprintf(stderr, "\t-v	Show version information.\n");
	fprintf(stderr, "\t-? 	Show the instruction. (For detailed instruction, please see the document here: https://github.com/mengyao/cTools)\n\n");
	fprintf(stderr, "Note:\ncmakefilter requires 20G memory. For LSF users, this is an example command to run cmakefilter: bsub -W 24:00 -R \"rusage[mem=20000]\" -o example.out -e example.err cmakefilter -p example.par\n\n");
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
 char *regnames[30], ss[256] ;
 int nfqnames = 5 ; 
 int  *filtval ;
 FILE *ffmask ;
 char *hengfilt, *cnvfilt ; 


// step 1 initialize regnames ;
 for (a=1; a<= 22; ++a) { 
  sprintf(ss, "%d", a) ; 
  regnames[a] = strdup(ss) ;
 }
 regnames[23] = strdup("X") ;
 regnames[24] = strdup("Y") ;
 regnames[25] = strdup("MT") ;

 //Read parameter file and initialize
 vcfname = NULL; // ref = NULL ;
 readcommands(argc, argv) ;
 
if (sampname == NULL) fatalx("no sampname\n") ;
if (hetfaname == NULL) fatalx("no hetfaname\n") ;
if (vcfsuffix == NULL) fatalx("no vcfsuffix\n") ;
if (maskname == NULL) { 
  sprintf(ss, "%s.mask.fa", sampname) ;  
  maskname = strdup(ss) ;
}

 openit(maskname, &ffmask, "w") ;

 nfqnames = 4 ;
 printf("sample name: %s  gender: %c\n", sampname, gender) ; 
 ZALLOC(fqnames, nfqnames, char *) ; 
 
/* 
 getdbname(fixeddbase, "Chimp", &fqnames[0]) ; refname = fqnames[0] ;
 getdbname(fixeddbase, "heng75", &fqnames[1]) ; refname = fqnames[0] ;
*/

	if (fixeddbase != NULL) {
  		getdbname(fixeddbase, "Chimp", &fqnames[0]) ; //refname = fqnames[0] ;
  		getdbname(fixeddbase, "heng75", &fqnames[1]) ; //refname = fqnames[0] ;
 	} else {
 		fqnames[0] = chimp;
 		fqnames[1] = heng75;
 	}
 	refname = fqnames[0];

 fqnames[2] = strdup(cnvname) ;
 fqnames[3] = sampfname = hetfaname ; 

 printstrings(fqnames, nfqnames) ;  fflush(stdout) ;

 nfregs = 3 ; 
 
// set lovals, hivals.  Kill all bases not meeting criterion
 ivclear(lovals,  0, nfregs) ;  
 ivclear(hivals, maxqval, nfregs) ;

 setlimv(lov, lovals, nfregs) ;
 setlimv(hiv, hivals, nfregs) ;

 for (k=0; k<nfregs; ++k) { 
  printf("limv: %2d  %3d %3d\n", k, lovals[k], hivals[k]) ; 
 }

 
if (debug==NO) {
 if (gender == 'M') {
  calcfiltermask(ffmask, regnames, 1, 22, 1, 22, fqnames, nfqnames, lovals, hivals) ; 
  calcfiltermask(ffmask, regnames, 23, 23, 23, 24, fqnames, nfqnames, lovals, hivals) ; 
  easymask(ffmask, regnames, 25, 25, '9') ;
 }
 else {
  calcfiltermask(ffmask, regnames, 1, 23, 1, 23, fqnames, nfqnames, lovals, hivals) ; 
  easymask(ffmask, regnames, 24, 24, '0') ;
  easymask(ffmask, regnames, 25, 25, '9') ;
 }
}
if (debug==YES) {
  calcfiltermask(ffmask, regnames, 21, 22, 21, 22, fqnames, nfqnames, lovals, hivals) ; 
}


 printf("totlen: %ld totvalid: %ld\n", totlen, totvalid) ;

 fclose(ffmask) ;

 sprintf(ss, "samtools faidx %s", maskname) ;
 system(ss) ;

  
 printf("## end of cmakefilter\n") ;
 return 0 ;

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
   printf("%10ld ", a) ; 
   printf("%10ld ", b) ; 
   printf("%12.6f", ymissfrac(a,b)) ;
   printnl() ;
   for (k=0; k<=9; ++k) { 
    printf("%8d ", k) ;            
    a = ss[k][0] ; b=ss[k][1] ; 
    tt[0] += a ; tt[1] += b ; 
    printf("%10ld ", a) ; 
    printf("%10ld ", b) ; 
    printf("%12.6f", ymissfrac(a,b)) ;
    printnl() ;
   }
   printf("%8s ", "Total:") ;  
   a = tt[0]  ; b = tt[1] ; 
   printf("%10ld ", a) ; 
   printf("%10ld ", b) ; 
   printf("%12.6f", ymissfrac(a,b)) ;
   printnl() ;
   printnl() ;

}

void easymask(FILE *ffmask, char **regnames, int rlo, int rhi, char maskval) 
// constant filter
{
  char **fqdata ;  
  char *fqnames[1] ;  
  char *ffout ; 
  int rnum, len ;

 // getdbname(fixeddbase, "Href", &fqnames[0]) ; refname = fqnames[0] ;

	if (fixeddbase != NULL) getdbname(fixeddbase, "Href", &fqnames[0]) ;
 	else fqnames[0] = href;

 	refname = fqnames[0] ;

  for (rnum = rlo; rnum <= rhi; ++rnum) {
   ZALLOC(fqdata, 1, char *) ;
   regname = regnames[rnum] ;
   readfa(fqnames, fqdata, &len, 1) ; 
   ffout = fqdata[0] ; 
   cclear(ffout, maskval, len) ;
   writefaf(ffmask, regname, ffout, len)  ;
   freeup(fqdata, 1) ;
  }

}

long loadtdata(long **tdata, int mq, int mq0, long *hhh, long *mmm, int len, int *bbase) 
{
  int ttt[3], a, kode  ;  
  long total = 0 ;

  lclear2D(&tdata, 100, 2, 0) ;
  for (kode = 0; kode < len; ++kode) {
   dekodeitbb(ttt, kode, 3, bbase) ; 
   total += hhh[kode] ;
   total += mmm[kode] ;
   if (ttt[1] != mq)  continue ;
   if (ttt[2] != mq0) continue ;
   a = ttt[0] ; 
   tdata[a][0] += hhh[kode] ;
   tdata[a][1] += mmm[kode] ;
  }
  return total ;
}

void calcfiltermask(FILE *ffmask, char **regnames, int rlo, int rhi, 
  int rloout, int rhiout, char **fqnames, int nfqnames, int *lovals, int *hivals)  
{

 long *hhh, *mmm ;  
 int nfregs = 3 ;
 int i, j, jbest, fval, lev ; 
 int *filtval ;
 int maxq[3] ; 
 int maxqbase[3] ; 
 int histsize, len ;  
 long **tdata ; 
 int a, b, c, d, x, t, k ;
 long cc, dd, tt, ttbest ;
 int *tlo, *thi ; 
 long *thh, *tmm  ;    
 long total = 0 ;  // total number of hits + misses
 int xlen ;
 double *tcover, *terate, *cover, *erate ;
 int maxlistlen = 1000*1000 ;
 FENTRY **felist, *fept, **febest ;
 int bbans[10] ; 
 double eeold[10], y1, ybest, zbest ;
  
 ivclear(maxq,  maxqval, 3) ; 
 maxq[2] = 1 ; 
 ivsp(maxqbase, maxq, 1, 3) ;

 histsize = iprod(maxqbase, 3) ; 
 printf("histsize: %d\n", histsize) ;
 ZALLOC(hhh, histsize, long) ;
 ZALLOC(mmm, histsize, long) ;

  ZALLOC(cover, maxlistlen, double) ; 
  ZALLOC(erate, maxlistlen, double) ; 
  ZALLOC(felist, maxlistlen, FENTRY *) ; 

  ZALLOC(tlo, maxlistlen, int) ; 
  ZALLOC(thi, maxlistlen, int) ; 
  ZALLOC(thh, maxlistlen, long) ; 
  ZALLOC(tmm, maxlistlen, long) ; 
  
  for (i=0; i<maxlistlen; ++i) {
   ZALLOC(felist[i], 1, FENTRY) ;
   fept = felist[i] ; 
   clearfe(fept) ;
  }


 len = mkhist(regnames, hhh, mmm, maxqbase, rlo, rhi, fqnames, nfqnames, lovals, hivals) ; 
 tdata = initarray_2Dlong(100, 2, 0) ; 
 total  = loadtdata(tdata, 60, 0, hhh, mmm, len, maxqbase) ;
  xlen = mktlist(tdata, 100, tlo, thi, thh, tmm) ;
  for (a=0; a<xlen; ++a) { 
   printf("zzq %3d %3d ", tlo[a], thi[a]) ;
   cc = thh[a] ; dd = tmm[a] ; tt = cc + dd ;
   fept = felist[a] ; 
   fept -> lovals[0] =  tlo[a] ;
   fept -> hivals[0] =  thi[a] ;
   fept -> hit   = cc ;
   fept -> miss  = dd  ;
   fept -> level = 1 ;
   fept -> cover = (double) tt / (double) total ;
   fept -> erate =  (double) dd / (double) tt ; 

   printf("%ld %ld ", fept -> hit, fept -> miss) ;  
   printf("%12.6f ", fept -> cover) ;                   
   printf("%12.6f", fept -> erate) ;
   printnl() ;
  }
// filter design here 
  setfpars(felist, xlen, bbans) ;
  ZALLOC(febest, 10, FENTRY *) ; 
  for (k=0; k<=9; ++k) {
   a = bbans[k] ;
   fept = febest[k] = felist[a] ;
   eeold[k] = fept -> erate  ;  
   ttbest = -1 ; 
  }
  for (b=1; b <=9; ++b) { 
   break; 
// break;  this attempt to use MQ < 60 seems to gain little
   loadtdata(tdata, 60-b, 0, hhh, mmm, len, maxqbase) ;
   xlen = mktlist(tdata, 100, tlo, thi, thh, tmm) ;
   for (k=0; k<=9; ++k) {  
    fept = febest[k] ;
    jbest = -1  ;
    ybest = eeold[k] ;
    for (j=0; j < xlen; ++j) {  // find j with max coverage and slope < eeold[k]  
     cc = thh[j]; dd = tmm[j] ; tt = cc + dd ;
     if (tt<1000) continue ;
     y1 = (double) dd / (double) tt ; 
     if ((y1<ybest) && (tt>ttbest)) {
      jbest = j ;  
      zbest = y1 ;
      ttbest = tt ;  
     }
    }
    if (jbest >= 0)  { 
     j = jbest ;
     printf("zzk: %d lev: %d ", k, b+1) ;
     printf("old erate: %12.6f  new slope: %12.6f\n", eeold[k], zbest) ;
     fept -> level = b+1 ; 
     fept -> lovals[b] = tlo[j] ;
     fept -> hivals[b] = thi[j] ;
     cc = thh[j]; dd = tmm[j] ; 
     fept -> hit  += cc ;
     fept -> miss += dd ;
     cc = fept -> hit; 
     dd = fept -> miss ;
     tt = fept -> hit + fept -> miss ;
     fept -> cover = (double) tt / (double) total ;
     eeold[k] = fept -> erate =  (double) dd / (double) tt ; 
    }
   }
  }
  printf("filter details:\n") ;
  for (k=0; k<=9; ++k) {
   printf("optf_details: %d ", k) ;
   fept = febest[k] ;
   printf("%12.6f ", fept -> cover) ;
   printf("%12.6f ", fept -> erate) ;
   lev = fept -> level ;
   for (b=0; b<lev; ++b) { 
    printf("  %2d ", fept -> lovals[b]) ;
    printf("%2d", fept -> hivals[b]) ;
   }
   printnl() ;
  }
  printf("filter summary:\n") ;
  for (k=0; k<=9; ++k) {
   printf("optfilt: %d ", k) ;
   fept = febest[k] ;
   printf(" %12ld ", fept -> hit) ;
   printf(" %12ld ", fept -> miss) ;
   printf(" %12ld ", fept -> hit + fept -> miss) ;
   printf("%12.6f ", fept -> cover) ;
   printf("%12.6f ", fept -> erate) ;
   printnl() ;
  }
  fflush(stdout) ;

  free(hhh) ;
  free(mmm) ;

  free(cover) ;
  free(erate) ;

  free(tlo) ;
  free(thi) ;
  free(thh) ;
  free(tmm) ;

	if (fixeddbase != NULL) getdbname(fixeddbase, "Href", &fqnames[0]) ;
	else fqnames[0] = href; 

	refname = fqnames[0] ;

 mkmask(ffmask, regnames, febest,  rloout, rhiout, fqnames, nfqnames) ; 

 for (i=0; i<maxlistlen; ++i) {
   free(felist[i]) ;             
 }
 free(felist) ;
}


double ymissfrac(double a, double b)
{
 double t = a + b ;

 if (t<=0.1) return  0 ;
  return  b / t ;

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

/**
   x = ranmod(100*1000) ; 
   if (x==0) { 
    printf("zzval %d ", pos);
    printf("%c %c %c %c\n", cref, csamp, ca, cb) ;
   }
*/
  
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
void mkvcfname(char *vcfname, char *dir, char *reg, char *suff) 
{
  sprintf(vcfname, "%s/%s.%s", dir, reg, suff) ; 
}

void paddata(char **fqd, int pi, int pj, char padval) 
// pad out fqd[pi] with char padval up to length fqd[pj] 
{
  int leni, lenj ; 
  char *fqtemp ;

  leni = strlen(fqd[pi]) ; 
  lenj = strlen(fqd[pj]) ; 

  if (lenj <= leni) return ; 
  ZALLOC(fqtemp, lenj+1, char)  ;  
  cclear(fqtemp, padval, lenj) ;
  fqtemp[lenj] = CNULL ; 
  strncpy(fqtemp, fqd[pi], leni) ;
  freestring(&fqd[pi]) ; 
  fqd[pi] = strdup(fqtemp) ; 
  
}

int  mkhist(char **regnames, long *hh, long *mm, int *blist, 
 int rlo, int rhi, char **fqnames, int nfqnames, int *lovals, int *hivals)  
{

 int **fasta ;  
 int nfregs = 3 ;
 char **fqdata ; 
 int *fqlen, len ; 
 int rnum, a, i, pos ;
 int *valid ; 
 int **vcfvals, vcflen[3], ttt[3] ;
 int *ccc, x, k, cc, t ;
 char cref, csamp ; 
 char vcfname[1024] ;
 int isbad, kode, vv ; 
 long tt, hit, miss ;
 long **mghit, **mgmiss ;

 ZALLOC(fqlen, nfqnames, int) ;

 for (rnum = rlo; rnum <= rhi; ++rnum) { 
  ZALLOC(fqdata, nfqnames, char *) ;

  regname = regnames[rnum] ;
  readfa(fqnames, fqdata, fqlen, nfqnames) ; 
  paddata(fqdata, 1, 3, '1') ;
  for (a=0; a<nfqnames; ++a) { 
   printf("fq: %3d %s %12d\n", a, fqnames[a], fqlen[a]) ;
  }
  ivmaxmin(fqlen, nfqnames, NULL, &len) ; // min fqlen  
// we have now read basic data ; 
  ZALLOC(valid, len, int) ;  
  setvalid(valid, len, fqdata, nfqnames) ;  
  freestring(&fqdata[1]) ; 
  freestring(&fqdata[2]) ; 

  vv = intsum(valid, len) ;
  printf("regname: %s len: %d  valid: %d\n", regname, len, vv) ; 
  totlen  += len ; 
  totvalid += vv ;
  vcfvals = initarray_2Dint(len, 3,  0) ;
  mkvcfname(vcfname, vcfdir, regname, vcfsuffix) ;
  readvcf(vcfname, vcfvals, len, vcflen) ;

  for (i=0; i<len; i++) { 
   if (valid[i] != 1) continue ;
   pos = i+1 ;
   ccc = vcfvals[i] ; 

    cref =  toupper(fqdata[0][i]) ;
    csamp = toupper(fqdata[3][i]) ; 

    if (!isiub2(cref)) continue ;
    if (!isiub2(csamp)) continue ;

    x = matchit(cref, csamp) ;

    isbad = NO ;
    for (k=0; k<3; ++k) { 
     ccc[k] = MIN(ccc[k], blist[k]-1)  ; 
     if (ccc[k] < 0) isbad = YES ; 
    }
    if (isbad) continue ; 
    kode = kodeitbb(ccc, 3, blist) ;

// x is number of hits 
     hh[kode] += x ; 
     mm[kode] += (2-x) ;

  }

  freeup(fqdata, nfqnames) ;
  free2Dint(&vcfvals, len) ;  
  free(valid) ;
 }

 free(fqlen) ;
 mghit  = initarray_2Dlong(3, maxhist, 0) ;
 mgmiss = initarray_2Dlong(3, maxhist, 0) ;

 len = iprod(blist, 3) ;
 for (kode = 0; kode <  len ; ++kode)  {  
  tt = hh[kode] + mm[kode] ;   
  if (tt==0) continue ;
  dekodeitbb(ttt, kode, 3, blist) ; 
  printf("hist: %3d %3d %d ", ttt[0], ttt[1], ttt[2]) ; 
  hit = hh[kode] ; 
  miss = mm[kode] ;
  printf ("%12ld ", hit) ;
  printf ("%12ld ", miss) ;
  printf(" %12.6f\n", ymissfrac(hit, miss)) ; 
  for (a=0; a<3; ++a) { 
   mghit[a][ttt[a]] += hit ;
   mgmiss[a][ttt[a]] += miss ;
  }
 }
 for (a=0; a<3; ++a) { 
  for (k=0; k<maxhist; ++k) {
   hit = mghit[a][k] ;
   miss = mgmiss[a][k] ;
   tt = hit + miss ;
   if (tt==0) continue ;
   printf("marginal: %3d %3d ", a, k) ;
   printf ("%12ld ", hit) ;
   printf ("%12ld ", miss) ;
   printf(" %12.6f\n", ymissfrac(hit, miss)) ; 
  }
 }
 fflush(stdout) ;
 free2Dlong(&mghit, 3) ;
 free2Dlong(&mgmiss, 3) ;
 return len ;
 
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
void *setvcff(char *vcffile, char **vcfn)   
// vcfffile should not exist
{
  char *wkdir = "trash" ;
  char iname[200], tmpname[200], outname[200], ss[1000] ;
  char *bname ;

 sprintf(ss, "mkdir -p %s ", wkdir) ;
 system(ss) ;
 sprintf(iname, "%s.gz", vcffile) ;  
  if (ftest(iname))  {
   printf("got %s\n", iname) ; 
  }
  else { 
   fatalx("gzfetch fails\n") ;
  }
  bname = basename(vcffile) ; 
  sprintf(outname, "%s/%s:%s.vcf", wkdir,sampname,regname ) ;
  sprintf(ss, "cp %s %s.gz", iname, outname) ; 
  printf("zzcopy: %s\n",ss) ;
  system(ss) ;
  sprintf(ss, "gunzip %s.gz", outname) ;
  system(ss) ;
  sprintf(ss, "chmod 664 %s", outname) ; 
  system(ss) ;
  if (ftest(outname) == NO) fatalx("unzip fails\n") ;

  *vcfn = strdup(outname) ;
  
}

int readvcf(char *vcffile, int **fasta, int reflen, int *flen) 
{

        printf("readvcf  file: %s reflen: %d\n", vcffile, reflen) ;
        FILE *vcffp;
        char *vcftmp = NULL, *ctmp ;
        int t ; 
        int depth = -1;
        int MQ = -1;
        int MQ0 = -1;
        int pos = 0;
	int i;
        int *ccc ;
// fasta is preallocated


        if (ftest(vcffile) == NO) { 
          setvcff(vcffile, &vcftmp) ;
          openit(vcftmp, &vcffp, "r") ; // aborts on error
        }
        else {
         openit(vcffile, &vcffp, "r") ; // aborts on error
        }
	char line[4096];

        iclear2D(&fasta, reflen, 3, -1) ;

	while ( fgets ( line, sizeof line, vcffp ) != NULL )  {
        depth = -1;
        MQ = -1;
        MQ0 = -1;
        pos = 0;

        if (line[0] != '#' && line[0] != '-')
        {
         char *fields[MAXFF], *infofields[MAXFF], *dp[MAXFF];
         int nfields = -1, ninfo = -1, ndp = -1;
         nfields = splitupx(line, fields, MAXFF, '\t') ;
         pos = atoi(fields[1]);
         if (pos > reflen) { 
          freeup(fields, nfields) ;
          continue ;
         }
         ccc = fasta[pos-1] ;
         ninfo = splitupx(fields[7], infofields, MAXFF, ';') ;
         int n ;
         for (n=0; n<ninfo; ++n) {
          ndp = splitupx(infofields[n], dp, MAXFF, '=') ;
          if (strcmp (dp[0], "DP") == 0) {
           if (isnumword(dp[1])) depth = atoi(dp[1]);
           ccc[0] = depth;
           //printf("DP: %i\n", depth) ;
          }
          if (strcmp (dp[0], "MQ") == 0) {
           if (isnumword(dp[1])) MQ = round(atof(dp[1]));
           ccc[1] = MQ;
           //printf("MQ: %.5f\n", MQ) ;
          }
          if (strcmp (dp[0], "MQ0") == 0) {
           if (isnumword(dp[1])) MQ0 = round(atof(dp[1]));
           ccc[2] = MQ0;
           //printf("MQ0: %.5f\n", MQ0) ;
          }
          freeup(dp, ndp);
         }
         freeup(fields, nfields);
         freeup(infofields, ninfo);
        }
/**
        t = ranmod(100*1000) ; 
        if (t==0) printf("zzvcf: %d  %d %d %d\n", pos, depth, MQ, MQ0) ;
*/
        
    }
    flen[0] = reflen;
    flen[1] = reflen;
    flen[2] = reflen;
    fclose ( vcffp );
    if (vcftmp != NULL) remove(vcftmp) ;

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

  while ((i = getopt (argc, argv, "p:vV?")) != -1) {

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
	//printf ("Usage: bad params.... \n") ;
	//fatalx("bad params\n") ;
      }
  }

   if (parname == NULL) exit(usage()); //return ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

  // getstring(ph, "regname:", &regname) ;
  // getstring(ph, "poplistname:", &poplistname) ;
   //getstring(ph, "iubfile:", &iubfile) ;
 //  getstring(ph, "mapstring:", &mapstring) ;
   getstring(ph, "cnv:", &cnvname) ;
// getstring(ph, "vcfname:", &vcfname) ;
  // getstring(ph, "ref:", &ref) ;
   getstring(ph, "gender:", &ss) ;
   if  (ss != NULL)  { 
    gender = ss[0] ;
    freestring(&ss) ;
   }
	t = 0 ; getint(ph, "lopos:", &lopos) ; lopos = MAX(lopos, t) ;
   t = BIGINT ; getint(ph, "hipos:", &hipos) ; hipos = MIN(hipos, t) ;

  // getstring(ph, "vcfbase:", &vcfbase) ;
   getstring(ph, "dbase:", &dbase) ;
   getstring(ph, "vcfdir:", &vcfdir) ;
	
	getstring(ph, "href:", &href) ;
    getstring(ph, "chimp:", &chimp) ;
    getstring(ph, "heng75:", &heng75) ;

   getstring(ph, "vcfsuffix:", &vcfsuffix) ;
   getstring(ph, "hetfa:", &hetfaname) ;
   getstring(ph, "fixeddbase:", &fixeddbase) ;
   getstring(ph, "sampname:", &sampname) ;
   getstring(ph, "popname:", &sampname) ;
   getstring(ph, "maskname:", &maskname) ;

   getstring(ph, "lovals:", &lov) ;
   getstring(ph, "hivals:", &hiv) ;
   getint(ph, "debug:", &debug) ;



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

void  mkmask(FILE *fff, char **regnames, FENTRY **fans,  int rlo, int rhi, char **fqnames, int nfqnames)
{

 int **fasta ;  
 int nfregs = 3 ;
 char **fqdata ; 
 int *fqlen, len, reflen ; 
 int rnum, a, i, j, pos ;
 int lo, hi ; 
 int *valid ; 
 int **vcfvals, vcflen[3] ;
 int *ccc, x, k, cc, t ;
 int fval ; 
 long **rstats ;            
 char cref, csamp, maskval ; 
 char *ffout ;
 FENTRY  *fept ;
 char vcfname[1024] ;
 int **ftable ; 
 int lev, dp, mq, mq0 ; 
 int ftok ; // look up ftable

 ftable = initarray_2Dint(100,10,0) ;
 iclear2D(&ftable, 100, 10, 0) ;  

 for (k=1; k<=9; ++k) { 
  fept = fans[k] ; 
  lev = fept -> level ;  
  for (j=0; j<lev; ++j) { 
   lo = fept -> lovals[j] ;
   hi = fept -> hivals[j] ;
   lo = MAX(lo, 0) ; 
   hi = MIN(hi, 99) ; 
   for (i=lo; i<=hi; ++i) { 
    ftable[i][j] = k ; 
   }
  }
 }

 ZALLOC(fqlen, nfqnames, int) ;

 rstats = initarray_2Dlong(11, 2, 0) ; 

 for (rnum = rlo; rnum <= rhi; ++rnum) { 
  ZALLOC(fqdata, nfqnames, char *) ;

  regname = regnames[rnum] ;
  readfa(fqnames, fqdata, fqlen, nfqnames) ; 
  paddata(fqdata, 1, 3, '1') ;
  for (a=0; a<nfqnames; ++a) { 
   printf("fq: %3d %s %12d\n", a, fqnames[a], fqlen[a]) ;
  }
  ivmaxmin(fqlen, nfqnames, NULL, &len) ; // min fqlen  
// we have now read basic data ; 
  ZALLOC(valid, len, int) ;  
  setvalid(valid, len, fqdata, nfqnames) ;  
  vcfvals = initarray_2Dint(len, 3,  0) ;
  mkvcfname(vcfname, vcfdir, regname, vcfsuffix) ;
  readvcf(vcfname, vcfvals, len, vcflen) ;

  reflen = fqlen[0] ;   
  ZALLOC(ffout, reflen+1, char) ;
  ffout[reflen] = CNULL ;
  cclear(ffout, 'N', reflen) ;

  for (i=0; i<len; i++) { 
   pos = i+1 ;
   ccc = vcfvals[i] ; 

    cref =  toupper(fqdata[0][i]) ;
    csamp = toupper(fqdata[3][i]) ; 

    if (!isiub2(cref)) continue ;
    if (!isiub2(csamp)) continue ;

    mq0 = ccc[2] ; 
    mq =  ccc[1] ; j=60-mq ; 
    dp = ccc[0] ;
    dp = MIN(dp, 99) ;
    ftok = YES ; 
    if (valid[i] != 1) ftok = NO ; 
    if (mq0 != 0 ) ftok = NO ; 
    if (j>9) ftok = NO ; 
    if (j<0) continue;  
    if (ftok == NO) fval = 10 ;
    else fval = ftable[dp][j] ;
    if (fval<0) continue ;
    x = matchit(cref, csamp) ;
    rstats[fval][0] += x ;
    rstats[fval][1] += (2-x) ;
    if (fval==10) continue ;
    if (ftok == NO) continue ;
    maskval = (char) (fval + '0') ; 
    ffout[i] = maskval ;
  }


  freeup(fqdata, nfqnames) ; 
  free(valid) ;
  free2Dint(&vcfvals, len) ;
  writefaf(fff, regname, ffout, reflen) ;
 }

 printhist(rstats) ;
 free2Dint(&ftable, 100) ;
 free2Dlong(&rstats, 11) ;

}

