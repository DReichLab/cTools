#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h>
#include <getpars.h> 


static int verbose = NO ; 
#define  VERSION  "100"  

double loghprob(int n, int a, int m, int k)  ;
int ranhprob(int n, int a, int m)  ;

int main() 
{
  int n=50,  m=16, a=5  ;
  int k, t, iter ;
  int hist[100] ;
  double y  ;

  ivzero(hist, 100) ;
  for (iter =1 ; iter <= 10000; ++iter) { 
   t = ranhprob(n, a, m) ;
   ++hist[t] ;
  }
  for (k=0; k<=a ; ++k) { 
   y = exp(loghprob(n, a, m, k)) ;
   printf("%3d %3d %9.3f\n", k, hist[k], y) ;
  }

}

double loghprob(int n, int a, int m, int k) 
// http://www.math.uah.edu/stat/urn/Hypergeometric.xhtml
{
/** 
 n balls a black.  Pick m without replacement  
 return log prob (k black)
*/

double ytop, ybot ;

if (k<0) return -1.0e30 ;
if (k>a) return -1.0e30 ;
if (k>m) return -1.0e30 ;
if ((m-k)>(n-a)) return -1.0e30 ;

 ytop = logbino(a, k)  + logbino(n-a, m-k) ;
 ybot = logbino(n, m) ;
 return ytop - ybot ;

}

int ranhprob(int n, int a, int m) 
// rejection sampling.  Devroye
{
 double v, y ;
 double pm, logpm, w, ru, rw, rat ;
 int mode, k, x, zans ;

 v = (double) (a+1)*(m+1) / (double) (n+1) ; 
 mode = (int) v ;

/**
 for (k=-5; k<=5; ++k) {  
  x = mode+k ;
  y = exp(loghprob(n, a, m, x)) ;
  printf("%4d %4d %12.6f\n", mode, x, y)  ;
 }
*/

 logpm = loghprob(n, a, m, mode) ;
 pm = exp(logpm) ;              
 w = 1 + pm ; 
 for (;;) { 
  ru = DRAND() ;
  rw = DRAND() ;
  if (ru <= w/(1+w)) y = DRAND()*w/pm ;
  else y = (w+ranexp())/pm ;
  x = nnint(y) ; 
  if (ranmod(2)==0) x = -x ;
  zans = mode+x ;
  if (zans<0) continue ;
  if (zans>a) continue ;
  rat = exp(loghprob(n, a, m, zans)-logpm) ; 
  rw *= MIN(1, exp(1.0-pm*y)) ;
  if (rw <= rat) break ;
 }
 return zans ;
 
}
