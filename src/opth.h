#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h> 
#include <getpars.h> 

#define MAXSTR  512


typedef struct {
  int lovals[10], hivals[10] ;
  double cover, erate ;
  long hit, miss, level ;
}  FENTRY ;  

void readcommands(int argc, char **argv) ;
void clearfe(FENTRY *fept)   ;
int mktlist(long **dd, int dlen, int *tlo, int *thi, long *thh, long *tmm) ;
void setfpars(FENTRY **felist, int felen,  int *bbans)  ;
