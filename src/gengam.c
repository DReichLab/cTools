#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  
#include <nicklib.h>  


int main(int argc, char **argv)
{
    int n ;  
    double p, lam ;
    int i, x ;

    n = atoi(argv[1]) ;
    p = atof(argv[2]) ;
  lam = atof(argv[3]) ;

  SRAND(77) ;

  printf ("## %4d %12.6f %12.6f\n", n, p, lam) ;
  for (i=0; i<n; i++)  {  
   x = rangam(p)/lam ;
   x == ranexp() ;
   printf("%d %12.6f\n", i, x) ;
 }
}
