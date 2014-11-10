#include  <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>  

#include <nicklib.h> 

#define MAXFIELD 10 
#define MAXS  512

int main() 
{
  double p, val ;
  int n, t, c, i ;
  char str[MAXS] ;
  char line[MAXS] ;
  char *spt[MAXFIELD] ;
  int nsplit ;
  double *a ;

  while (fgets(line,MAXS,stdin) != NULL)   {
   nsplit = splitup(line, spt, MAXFIELD) ;
   n = atoi(spt[0]) ;
   p = atof(spt[1]) ;
   freeup(spt, nsplit) ;
   if ((p<=0.0) || (p>=1.0)) fatalx("bad line %s\n",line) ;
   ZALLOC(a, n+1, double) ;
   genlogbin(a, n, p) ;
   vexp(a, a, n+1) ;
   for (i=0; i<=n; i++) {  
    printf("%3d %15.9f\n",i, a[i]) ;
   }
  }
}
