#include <stdio.h>
#include <math.h> 

#include "statsubs.h" 
#include "vsubs.h" 

#define EPS1 .001 
#define EPS2 1.0e-12 
#define ZLIM 20 
#define QSIZE 10 


static double *bern ; /* bernouilli numbers */
static int bernmax = 0 ;
static double *ztable = NULL, *ptable = NULL ;
static double  ptiny  ;
static int numbox = QSIZE*ZLIM ;

static double zzprob(double zval) ;
static double znewt(double z, double ptail) ;

static double ltlg1(double a, double x)  ;
static double ltlg2(double a, double x)  ;
static double rtlg1(double a, double x)  ;
static double rtlg2(double a, double x)  ;
static double pochisq (double x, int df) ;
static double pof (double F, int df1, int df2) ;

void cinterp(double val, double x0, double x1, 
  double f0, double f0p, double f1, double f1p, double *fv, double *fvp)  ;
int firstgtx(double val, double *tab, int n) ;
static int gtx(double *tab, int lo, int hi, double val)  ;
void gettw(double x,  double *tailp, double *densp)   ;

static int twtabsize = -1 ;
static double *twxval, *twxpdf, *twxtail ;

double twdensqq(double twstat) 
// Tracy-Widom prob density
{
  double dens, tail ;

  gettw(twstat, &dens, &tail) ; 
  return dens ;

}

double twtailqq(double twstat) 
// Tracy-Widom right tail
{

  double dens, tail ;
  static int ncall = 0 ;  
    

  ++ncall ;

  gettw(twstat, &tail, &dens) ; 
/**
  printf("zz %9.3f %9.3f\n", twstat, tail) ;
  if (ncall==10) abort() ;
*/
  return tail ;

}

static int fgtx(double *tab, int lo, int hi, double val) 
{

    int k ;
    
    if (val >= tab[hi]) return hi+1 ;  
    if (val < tab[lo])  return lo ;
    k = (lo+hi)/2 ;
    if (val <= tab[k]) return fgtx(tab, lo+1, k, val) ; 
    return fgtx(tab, k, hi-1, val) ;
}

int firstgtx(double val, double *tab, int n)
// tab sorted in ascending order 
{
 return fgtx(tab, 0, n-1, val) ;
}

void
gettw(double x, double *tailp, double *densp)   

{
     int k, n  ; 
     double x0, x1, f0, f1, f0p, f1p ;  
     double *xx[3] ;


  if (twtabsize = -1)  {
    k = numlines("TWTAB") ;
    ZALLOC(twxval, k, double) ;
    ZALLOC(twxpdf, k, double) ;
    ZALLOC(twxtail, k, double) ;
    xx[0] = twxval ;
    xx[1] = twxtail ;
    xx[2] = twxpdf ;
    twtabsize = getxx(xx, k, 3, "TWTAB") ;
  }
  n = twtabsize ;

     k = firstgtx(x, twxval, n) ;    
     
     if (k<=0) {  
       *tailp = 1.0 ; 
       *densp = 0.0 ; 
       return ;
     }

     if (k>=n) { 
       *tailp = twdensx(x)  ;
       *densp  = twtailx(x) ;
       return ; 
     }

     x0 = twxval[k-1] ; 
     x1 = twxval[k] ;  
     f0  = twxtail[k-1] ;
     f0p = twxpdf[k-1] ;
     f1 =  twxtail[k] ;
     f1p = twxpdf[k] ;

// now do cubic interpolation
     cinterp(x, x0, x1, 
      f0, -f0p, f1, -f1p, tailp, densp) ;
      *densp = - *densp ;

/**
     printf("zzz %9.3f %9.3f %9.3f\n", x0, x1, x) ;
     printf("zz1 %9.3f %9.3f %9.3f\n", f0, f1, *tailp) ;
     printf("zz2 %9.3f %9.3f %9.3f\n", f0p, f1p, *densp) ;
*/

}

void cinterp(double val, double x0, double x1, 
  double f0, double f0p, double f1, double f1p, double *fv, double *fvp) 
// cubic interpolation val should be between x0 and x1
// fv is function at x 
// fvp is derivative

{

    double inc, yval, f, fp, a0, b0, a1, b1, a2, a3 ; 
    double c0, c1, cc0, cc1 ;

    inc = x1-x0 ;
    yval = (val-x0)/inc ;
    b0 = f0 ;
    b1 = f0p*inc ;  
    c0 = f1 ; 
    c1 = f1p*inc ;
 
    a0 = b0 ; 
    a1 = b1 ; 
    cc0 = c0 - (a0+a1) ; 
    cc1 = c1 - a1 ;  
    a2 = 3*cc0-cc1 ; 
    a3 = cc1 - 2*cc0 ; 

    f = a3 ; 
    f *= yval ; 
    f += a2 ; 
    f *= yval ; 
    f += a1 ; 
    f *= yval ; 
    f += a0 ;
    *fv = f ;                        
    fp = 3*a3 ; 
    fp *= yval ; 
    fp += 2*a2 ; 
    fp *= yval ; 
    fp += a1 ;

    *fv = f ; 
    *fvp = fp/inc ;

}

