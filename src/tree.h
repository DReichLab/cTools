
#define MAXMUT 1024
typedef struct {
  int pos ; 
  char alleles[2] ;
  int  mtype ;
  char pat[40] ;
 struct  NODE *qtree ;
}  MUTATE ;

typedef struct { 
  int pnum ; 
  int num ;
  int lev ;
  char name[40] ;
  int child[MAXMUT] ;
  int nchild ;
  MUTATE *muts ; 
  int nmuts ;
} NODE ;


