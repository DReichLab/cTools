#include <nicklib.h>
#include "faidx.h" 

#ifndef NJPSAM 

#define MAXD    100
#define MAXDD   10
#define MAXPOP  299   
#define MAXPOPX 300
#define MAXLEN  2000 

typedef struct  {
 char *faname ; 
 char *alias ;
 char *famask ;
 int lopos ; 
 int hipos ; 
 int len ;
 int rlen ;
 int mlen ; 
 int popnum ;
 int isref ;
 char *regname ;
 char *rstring ;
 char *mstring ;
 faidx_t *fai ;
 faidx_t *faimask ;
} FATYPE ;

typedef struct {
  char qname[MAXLEN] ; 
  char pop[MAXLEN] ; 
  int ipop ;
  int ilib ;
  int ihdr ;
  int isamp ;
  int flag ;
  int isnean ;
  int isxw ;
  int isancient ;
  int strand ;
  char rname[MAXLEN] ; 
  int pos ;
  int mapq ; 
  char cigar[MAXLEN]; 
  char mrname[MAXLEN] ; // null => rname "=" => rname
  int mpos ;
  int isize ;
  int bonenumber ; 
  char seq[MAXLEN] ; 
  char qual[MAXLEN] ; 
  int  iqual[MAXLEN] ;
  int  diff[MAXLEN] ;
  int  readlen ;  // strlen(seq) = strlen(qual)  
  char refseq[MAXLEN] ;    
  int  posdiff[MAXLEN] ; 
  int numposdiff ;
  char lib[MAXLEN] ;
  char samp[MAXLEN] ;
  char rgstring[MAXLEN] ;
  int mapthresh ;
  int lobasequal ;
  int ispair ;
  int isdup ;
  double pmdscore ;
} READ  ;

typedef struct {
 char *chrom ; 
 int  cpos ;
 int  len  ; 
 char *seq ; 
 int  *aqual ; 
}  AREC ; 

typedef struct {
 char refbase ;
 int ipop ;
 int pos ;
 char bases[MAXD] ;
 char base ;
 char last ; 
 char next ;
 int sbase ;  // strand for base
 int  qual[MAXD] ;
 int  diff[MAXD] ; 
 int strand[MAXD] ;
 double pmd[MAXD] ;  
 int  depth ;
 int  bonenumber[MAXD] ;
 int  libnumber[MAXD] ;
 char strandbase[2] ;
 int  mapq[MAXD] ;
 int pickmapq ;
 int pickbonenumber ;
 int pickindex ;
} BASEP ;

typedef struct { 
 char refbase ;
 char refnext ;
 char reflast ;
 int  pos ;
 BASEP *bp[MAXPOPX] ;
 int isprocessed ; 
} PILEUP ;

typedef struct { 
 int lomapthresh ; 
 int himapthresh ; 
 int quallo ; 
 int qualhi ; 
 int lobasedepth ;
 int hibasedepth ;
 int skiplength ;
 int ispair ;
 int setthresh ; 
 double pmdlo ;
 double pmdhi ;
}  THRESH ;

#endif 
#define NJPSAM  

