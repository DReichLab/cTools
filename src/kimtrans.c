#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  
#include <nicklib.h>  

#include "kimsubs.h" 

int verbose = NO ;
static int debug = 0 ;

void setwynn(int val) ; 
void evalkim00(double x, double p) ;
double chenstrook(double x, double y, double t) ;
#define N 1000 

/** 
 simply prinbt Kimura function for given x, t 
*/
int main  (int argc , char **argv) 
{
   double x, t, y, ktrans, cs ;
   int j, wynnval = YES ;
  
/**
  evalkim00(0.01, 0.02) ;
  return 0 ;
*/
  x = atof(argv[1]) ;
  t = atof(argv[2]) ;
  if (argc>3) wynnval = atoi(argv[3]) ;
  setwynn(wynnval) ;
  setmaxpol(800) ;
   
  printf("## x:   %10.4f  t:  %10.4f\n", x, t) ;
  for (j=0; j<=N; j++) { 
   y = (double) j ;  y /= (double) N ; 
   verbose = NO ;
// if (j==0) verbose = YES ;
   ktrans = evalkim(1-y, 1-x, t, 1.0) ;
   cs = chenstrook(1-x, 1-y, t) ;
   printf("%12.6f  %12.6f %12.6f\n", y, ktrans, cs) ;
  }
  return 0 ;
}


/**
   Copyright (c) 1998 
   Kapteyn Institute Groningen
   All Rights Reserved.
*/

double logqfun(double z) 
{

// q[z_] := Sqrt[z] BesselI[1, 2*Sqrt[z]];

  double logbessi1( double x), logbessi1(double x) ;
  double y  ;

  if (z==0.0) return -1.0e30 ;
  y = sqrt(z) ;
  return log(y) + logbessi1(2*y) ;
}

double qfun(double z) 
{

// q[z_] := Sqrt[z] BesselI[1, 2*Sqrt[z]];
// derivative at 0 is 1

  double bessi1( double x), logbessi1(double x) ;
  double y  ;



  if (z==0.0) return  0 ;     
  y = sqrt(z) ;
  return y*bessi1(2*y) ;


}

double logbessi1( double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double ax,ans,logans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
      logans = log(ans) ;
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      logans = log(ans) ; 
      logans += ax ; logans -=  log(ax)/2 ;
   }
   return logans ; 
}




double bessi1( double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
   }
   return x < 0.0 ? -ans : ans;
}
double ps(double x) 
{
 double y  ; 
 if (x==0.0) return 0 ;

 y = asin(sqrt(x)) ;
 return y*y ;

}
double dps(double x) 
// deriv(ps(x))
{
 double y  ; 
 if (x==0.0) return 1 ;

 y = asin(sqrt(x)) ;
 y /= sqrt(x*(1-x)) ;
 return y ;  

}


double chenstrook(double x, double y, double t) 
{

 double tt = t/2 ;
 double px, py, dpx, dpy ;
 double qarg, earg, top1, top2, bot1, bot2 ;
 double ans, z ;

/**
YS: Mathematica 
psi[x_] := (ArcSin[Sqrt[x]])^2;
dpsi[x_] := ArcSin[Sqrt[x]]/(Sqrt[1 - x] Sqrt[x]);
q[z_] := Sqrt[z] BesselI[1, 2*Sqrt[z]];
ChenStroock[x_, y_, t_] := 1/(y (1 - y))*Exp[-(psi[x] + psi[y])/(t/2)]/Sqrt[dpsi[x]*dpsi[y]]*q[psi[x]*psi[y]/(t/2)^2];
*/

 if (x>0.5) return chenstrook(1-x, 1-y, t) ;

 if (y==1.0) return chenstrook(x, .999999, t) ;
 px = ps(x) ; py = ps(y) ;
 dpx = dps(x) ; dpy = dps(y) ;

 qarg = (px*py)/(tt*tt) ;
 earg = -(px+py)/tt ; 
 if (qarg > 10) { 
  z = logqfun(qarg) ;
  top1 = 1.0 ; 
  top2 = exp(z+earg) ;
 }
 else {
  top1 = exp(earg) ;
  top2 = qfun(qarg) ;
 }
 bot1 = y*(1-y) ;
 bot2 = sqrt(dpx*dpy) ;

 if (y==0.0) { 
   top2 = px/(tt*tt) ;  // deriv wrt y 
   bot1 = 1.0  ;  // deriv wrt y 
 }
 if (y==1.0) { 
  printf("zzz %12.6f %12.6f %12.6f %12.6f\n", top1, top2, bot1, bot2) ;
  printf("zz2 %12.6f %12.6f %12.6f %12.6f\n", px, py, qarg, qfun(qarg)) ;
 }

 ans = (top1*top2)/(bot1*bot2) ;
 return ans ;
  



}


