#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h>
#include <getpars.h> 

static int verbose = NO ; 
#define  VERSION  "100"  

#define N 4 

int main() 
{
   int x = 7, xx ;   
   char c ;

   c = int2char(x) ; 
   xx = char2int(c) ;

   printf("%c %d\n", c, xx) ;
}
