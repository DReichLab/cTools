#include <nicklib.h>    
#include "getpars.h" 


int   
main  (int argc , char **argv) 
{
   phandle *ph ;
   char *parname ; 
   parname = strdup("actpar") ;
   ph = openpars(parname) ;
   writepars(ph) ;
   printf("\n\n") ;
   dostrsub(ph) ;
   writepars(ph) ;
   closepars(ph) ;
}
