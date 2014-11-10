#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  
#include <nicklib.h>  

#include "kimsubs.h" 

static int verbose = NO ;
static int debug = 0 ;

#define N  9 
void setccc(double *x, int n) ;
void setggg(double *x, int n) ; 
void getcc(double *cc, int n, int a, int b) ;
double *cco, *gco, *nn, *lvec ;      
int nmax = -1 ; 
void m2t(double *ft, double *fm, int n) ;
void t2m(double *gt, double *gm, int n) ;
void dodrift(double *gt, double *ft, int n, double ttau) ;
void setnn(double *lv, int n) ;
void setlv(double *lv, int n) ;

void free_poly(int *poly, double *cpoly)  ; 
int normalize_poly(int **ppoly, double **pcpoly, int n, int numterm) ;
int multiply_poly(int **ppoly, double **pcpoly, int *p1, double *c1, int *p2, double *c2, int n, int n1, int n2) ;
void print_poly(int *poly, double *cpoly, int n, int numterm)  ;
void copy_poly(int *pa, double *ca, int *pb, double *cb, int n, int numterm)  ;
void sort_poly(int **ppoly, double **pcpoly, int n, int numterm, int *order)  ;                      
int liftx(double *wout, double *wval, int *wexp, int lenw, double tau) ;

int kodeitbb(int *xx, int len, int *baselist)  ;
int dekodeitbb(int *xx, int kode, int len, int *baselist) ;

int main  (int argc , char **argv) 
{
  int n=3, maxd=3, numterm=1000, n1, n2 ; 
  int *dmaxa, *poly, *p1,   *p2, *tpp ; 
  double *cpoly, *c1, *c2, *tcc ;
  int i, px ;
  double tau ;

  px = numterm*n ; 
  ZALLOC(poly, px, int) ;
  ZALLOC(cpoly, numterm, double) ;
  ZALLOC(dmaxa, n, int) ;
  ivclear(dmaxa, maxd, n) ;
  for (i=0; i<px; i++) { 
   poly[i] = ranmod(maxd+1) ;
  }
  for (i=0; i<numterm; i++) { 
   cpoly[i] = gauss() ;
  }
  numterm = normalize_poly(&poly, &cpoly, n, numterm) ;
  print_poly(poly, cpoly, n, numterm) ;
  free_poly(poly, cpoly) ;

  n1 = random_poly(&tpp, &tcc, dmaxa, n) ;
  ZALLOC(p1, n1*n, int) ; ZALLOC(c1, n1, double) ;
  copy_poly(tpp, tcc, p1, c1, n, n1) ;

  n2 = random_poly(&tpp, &tcc, dmaxa, n) ;
  ZALLOC(p2, n2*n, int) ; ZALLOC(c2, n1, double) ;
  copy_poly(tpp, tcc, p2, c2, n, n2) ;

  printf("p1 p2:\n") ;
  print_poly(p1, c1, n, n1) ;
  print_poly(p2, c2, n, n1) ;

  numterm = multiply_poly(&tpp, &tcc, p1, c1, p2, c2, n, n1, n2)  ;
  free_poly(p1, c1) ;
  free_poly(p2, c2) ;
  printf("zz3\n") ;
  printmat(tcc, 1, 10) ;
  ZALLOC(poly, numterm*n, int) ; ZALLOC(cpoly, numterm, double) ;
  copy_poly(tpp, tcc, poly, cpoly, n, numterm) ;
  print_poly(poly, cpoly, n, numterm) ;
  tau = 0.1 ;
  numterm = lift_poly(&poly, &cpoly, n, numterm, 1, tau) ;
  return 0 ;
}
int lift_poly(int **ppoly, double **pcpoly, int n, int numterm, int kvar, double tau) 
/* Kimura lift variable kvar with drift tau */
{ 
 int *ww, *w1, *w2, *wexp, d, x, i, t, z, k, j, nexp, xout ; 
 int *wmax, *dmax, maxterm, maxdeg ;
 int *poly,  *tpoly, *pt ; 
 double *cpoly, *tcpoly ;
 int *windex ;
 double  *wval, *wout  ;

 ZALLOC(ww, n, int) ;
 ZALLOC(w1, n, int) ;
 ZALLOC(w2, n, int) ;
 ZALLOC(wmax, n, int)  ;
 ZALLOC(dmax, n, int)  ;

 ivclear(wmax, -1, n) ;
 idperm(ww, n) ; 
 ww[kvar] = n+1000 ; 
 sort_poly(ppoly, pcpoly, n, numterm, ww) ;
 printf("yysort\n") ;
 print_poly(*ppoly, *pcpoly, n, numterm) ;
 poly = *ppoly ; 
 cpoly = *pcpoly ;
 
 for (k=0; k<numterm; ++k) { 
  pt = poly + k*n ;
  for (j=0; j<n; ++j) { 
   wmax[j] = MAX(wmax[j], pt[j]) ;
  }
 }
 ivsp(dmax, wmax, 1, n) ;
 maxterm = iprod(dmax, n) ;
 ZALLOC(tpoly, maxterm*n, int) ;
 ZALLOC(tcpoly, maxterm, double) ;
 ivmaxmin(dmax, n, &maxdeg, NULL) ;

 ZALLOC(wout, maxdeg, double) ;
 ZALLOC(wval, maxdeg, double) ;

 k = 0 ;
 copyiarr(poly+k*n, w1, n) ;
 z = w1[kvar] ;
 w1[kvar] = 0 ;
 x = xout = 0 ;
 wexp[x] = z ;
 wval[x] = cpoly[k] ;
 ++x ;
 
 for (k=1; k<numterm; ++k) { 
  copyiarr(poly+k*n, w2, n) ;
  z = w2[kvar] ;
  w2[kvar] = 0  ;
  t = compiarr(w1, w2, n) ;
  if (t==0) { 
   wexp[x] = z ;
   wval[x] = cpoly[k] ;
   ++x ;
  }
  else { 
   nexp = liftx(wout, wval, wexp, x, tau) ;
   for (d=0; d<=nexp; ++d) {  
    w1[kvar] = d ;
    copyiarr(w1, tpoly+xout*n, n) ;
    tcpoly[xout] = wout[d]  ;
    ++xout  ;
    copyiarr(w2, w1, n) ;
    w1[kvar] = 0 ;
    x = 0 ;
    wexp[x] = z ;
    wval[x] = cpoly[k] ;
    ++x ;
   }
  }
 }
 nexp = liftx(wout, wval, wexp, x, tau) ;
 for (d=0; d<=nexp; ++d) {  
   w1[kvar] = d ;
   copyiarr(w1, tpoly+xout*n, n) ;
   tcpoly[xout] = wout[d]  ;
   ++xout  ;
 }

 x = xout ;
 sort_poly(&tpoly, &tcpoly, n, x, NULL) ;
 free(ww) ;
 free(w1) ;
 free(w2) ;
 free(wmax) ;
 free(dmax) ;
 free(wout) ;
 free(wval) ;

 free_poly(*ppoly, *pcpoly) ;
 ZALLOC(*ppoly, x*n, int) ;
 ZALLOC(*pcpoly, x, double) ;
 copy_poly(tpoly, tcpoly, *ppoly, *pcpoly, n, x) ;
 free_poly(tpoly, tcpoly) ;

 return x ;

}
int liftx(double *wout, double *wval, int *wexp, int lenw, double tau) 
{

 int nexp, i, k ; 

 ivmaxmin(wexp, lenw, &nexp, NULL) ;

 vzero(wout, nexp+1) ;
 for (i=0; i<lenw; i++) {   
  k = wexp[k] ; 
  wout[k] = wval[i] ;
 }
 liftit(wout, wout, nexp, tau) ;
 return nexp ;


}
void sort_poly(int **ppoly, double **pcpoly, int n, int numterm, int *order)                         
/* order = NULL os OK */
{
  int **ipt, *ind, *pt ; 
  int *poly, *wpoly ; 
  double *cpoly, *wcpoly ; 
  int j, t, k ; 
  double y ;

  static int *tpoly ;
  static double *tcpoly ;
 
 poly = *ppoly  ;
 cpoly = *pcpoly  ;

 ZALLOC(ipt, numterm, int *) ;
 ZALLOC(ind, numterm, int) ;
 pt = poly ;
 for (k=0; k<numterm; ++k) {
  ipt[k] = pt ;
  pt += n ;
 }

 ipsortitp(ipt, ind, numterm, n, order) ;

 ZALLOC(wpoly, n*numterm, int) ;
 ZALLOC(wcpoly, numterm, double) ;

 for (k=0; k<numterm; ++k) { 
  pt = ipt[k] ;  
  j = ind[k] ;
  y = cpoly[j] ;
  copyiarr(pt, wpoly+k*n, n) ;
  wcpoly[k] = y ;
 }

 free_poly(*ppoly, *pcpoly) ; 
 free(ipt) ; 
 free(ind) ; 


 ZALLOC(*ppoly, numterm*n, int) ;
 ZALLOC(*pcpoly, numterm, double) ;
 copy_poly(wpoly, wcpoly, *ppoly, *pcpoly, n, numterm) ;
 free_poly(wpoly, wcpoly) ;

}




void copy_poly(int *pa, double *ca, int *pb, double *cb, int n, int numterm) 
// pb cb preallocated
{
 copyiarr(pa, pb, n*numterm) ;
 copyarr(ca, cb, numterm) ;


}

int  random_poly(int **ppoly, double **pcpoly, int *dmax, int n) 
{
 int *wmax  ; 
 static int *tpoly = NULL ;
 static double *tcpoly = NULL;
 int numterm, x ;

 ZALLOC(wmax, n, int) ;
 ivsp(wmax, dmax, 1, n) ;
 numterm = iprod(wmax, n) ;

 free_poly(tpoly, tcpoly) ;
 ZALLOC(tpoly, numterm*n, int) ;
 ZALLOC(tcpoly, numterm, double) ;

 for (x=0; x<numterm; ++x) { 
  dekodeitbb(tpoly+x*n, x, n, wmax) ;
  tcpoly[x] = gauss() ;
 }


 *ppoly  = tpoly ;
 *pcpoly = tcpoly ;

 free(wmax) ;
 return numterm ;

}

void free_poly(int *poly, double *cpoly)  
{
 if (poly != NULL) free(poly) ;
 if (cpoly != NULL) free(cpoly) ;
 poly = NULL ;
 cpoly = NULL ;
}
void print_poly(int *poly, double *cpoly, int n, int numterm) 
{
 int i, j, t ;


 for (i=0; i<numterm; i++) { 
  printf("%5d: ", i) ;
  for (j=0; j<n; j++) { 
   t = poly[i*n+j] ; 
   printf("%3d ", t) ;
  }
  printf ("  %12.6f", cpoly[i]) ;
  printnl() ;
 }

}

int normalize_poly(int **ppoly, double **pcpoly, int n, int numterm) 
{
/** 
 takes polynomial withe duplicated terms and returns sorted normalized poly
 n number of variables, numtern number of terms
*/
  int **ipt, *ind, *pt ; 
  int *poly, *wpoly ; 
  double *cpoly, *wcpoly ; 
  int j, t, k, x ; 
  double y ;

  static int *tpoly ;
  static double *tcpoly ;
 
 poly = *ppoly  ;
 cpoly = *pcpoly  ;
 ZALLOC(ipt, numterm, int *) ;
 ZALLOC(ind, numterm, int) ;
 pt = poly ;
 for (k=0; k<numterm; ++k) {
  ipt[k] = pt ;
  pt += n ;
 }
/**
 printf("zznosort %d %d\n", numterm, n) ;
 printimat(ipt[0], 1, n) ;
 printimat(ipt[1], 1, n) ;
 printnl() ;
*/

 ipsortit(ipt, ind, numterm, n) ;

/**
 printimat(ind, 1, 10) ;
 printf("zzsort\n") ;
 printimat(ipt[0], 1, n) ;
 printimat(ipt[1], 1, n) ;
 printimat(ipt[numterm-1], 1, n) ;
 printnl() ;
*/

 ZALLOC(wpoly, n*numterm, int) ;
 ZALLOC(wcpoly, numterm, double) ;

 pt = ipt[0] ;
 x = 0 ;
 copyiarr(pt, wpoly+x*n, n) ;
 for (k=0; k<numterm; ++k) { 
  j = ind[k] ;
  y = cpoly[j] ;
  t = compiarr(ipt[k], pt, n) ; 
  if (t!=0) { 
   ++x ; 
   pt = ipt[k] ;
   copyiarr(pt, wpoly+x*n, n) ;
  }
  wcpoly[x] += y ;
 }

 free_poly(*ppoly, *pcpoly) ; 
 free(ipt) ; 
 free(ind) ; 

 ++x ;  // numterm

 ZALLOC(*ppoly, x*n, int) ;
 ZALLOC(*pcpoly, x, double) ;
 copy_poly(wpoly, wcpoly, *ppoly, *pcpoly, n, x) ;
 free_poly(wpoly, wcpoly) ;

 return x  ;
}
int multiply_poly(int **ppoly, double **pcpoly, int *p1, double *c1, int *p2, double *c2, int n, int n1, int n2) 
{
  int **ipt, *ind, *pt ; 
  int i, j, t, k, x, numpterm ; 
  double y ;

 static int *tpoly = NULL ;
 static double *tcpoly = NULL;
 
 numpterm = n1*n2 ;

 free_poly(tpoly, tcpoly) ;
 ZALLOC(tpoly, numpterm*n, int) ;
 ZALLOC(tcpoly, numpterm, double) ;


 x = 0 ; 
 for (i=0; i<n1; i++) { 
  for (j=0; j<n2; j++) { 
   ivvp(tpoly+x*n, p1+i*n, p2+j*n, n) ;
   tcpoly[x] = c1[i]*c2[j] ;
   ++x ;
  }
 }

 
 printf("zz1 %d %d\n", x, n) ;
 printimat(tpoly, 10, n) ;
 printmat(tcpoly, 1, 10) ;
 x  = normalize_poly(&tpoly, &tcpoly, n, x) ;
 printf("zz2 %d\n", x) ;
 printimat(tpoly, 10, n) ;
 printmat(tcpoly, 1, 10) ;

 ZALLOC(*ppoly, x*n, int) ;
 ZALLOC(*pcpoly, x, double) ;
 copy_poly(tpoly, tcpoly, *ppoly, *pcpoly, n, x) ;
 free_poly(tpoly, tcpoly) ;
 return x ;

}



int kodeitbb(int *xx, int len, int *baselist)  
{
  int t = 0 , i, base ;

  for (i=0; i<len; ++i) {
   base = baselist[i] ;
   t *= base ;
   t += xx[i] ;
  }
  return t ;
}

int dekodeitbb(int *xx, int kode, int len, int *baselist)
{
// return weight

  int i, t, base ;

  t = kode ;
  for (i=len-1; i>=0; --i) {
   base = baselist[i] ;
   xx[i] = t % base ;
   t /= base ;
  }
  return intsum(xx, len) ;

}


