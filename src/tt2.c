#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  
#include <strsubs.h> 
#include <vsubs.h>
#include <sortit.h>

#include <ranmath.h> 
#include "gsubs.h" 
#include "statsubs.h" 

double z2prob(double p)  ;

main() 
{
    double z1,z2,z3,p ;
    int i ;
    for (i=0; i<15; i++) {  
      z1 = (double) i + DRAND() ;
      p = ntail(z1) ;
      z2 = zprob(p) ;
      z3 = z2prob(p) ;
      printf("%g %g %g %g\n",z1,z2,z3,p-ntail(z3));  
    }
}
double z2prob(double p) 
{
   double za, zb, h, err, eabs, ltop, lbot ;
   double pi, piq, piql ;
   int iter ;

   pi = acos(0.0) * 2.0 ;
   piq = sqrt(2.0*pi) ;
   piql = log(piq) ;
   
   za = zprob(p) ;
   for (iter=0; iter<10; ++iter)  { 
    err= p-ntail(za) ;
    if (err==0.0) return za ;
    eabs = fabs(err) ; 
    lbot = (-0.5*za*za ) -piql ;  
    ltop = log(eabs) ;
    h = exp(ltop-lbot) ;
    if (err>0) h = -h ;
    printf("%3d %20.10f  %20.10f  %20.10f\n", iter, za, za+h, err) ;
    za = za+h ;
  }
  return za ;
}




