#include "nicksam.h" 


void printbadrgs()  ;
void freebonetable()  ;
int loadbonetable(char *bfname)  ;
int loadstr(char **strlist, char **xlist, int numin) ;
int bone2num(char *ss) ;
int library2num(char *ss) ;
int crackhdr(char *hdr, char *lib, char *samp, char *pop) ;
int getbone(char *library, char *bone, char *hdr, char *bonefile) ;
int getxstr(char ***xx, int maxrow, int numcol, char *fname) ;
int indxsub(char *ss, char **table, int tablen) ;
int check(char *sa, char *sb) ;
char *popstring(char **popl, int n) ; // first character

void loadfilebase(char *parname) ;
void unloadfilebase() ;

void freerdtable()  ;
int loadrdtable(char *rdname)  ;
int whichhdr(char *s1) ;
int whichlib(char *s1) ;
int whichpop(char *s1) ;
void fixrgstring(char *rgs) ;
int checkrg(char **rawp) ;
int oldwhichpop(char *s1, char *lib, char *samp, char *pop) ;


int base2num(char c) ;
char num2base(int k) ;
int kodex (int *aa, int len)  ;
int kode4 (int *aa)  ;
void setqhack(int val) ;
void dekodex (int kode, int *aa, int len)  ;
void dekode4 (int kode, int *aa)  ;
int kodexb (int *aa, int len, int base)  ;
int setpmdscore(int val) ;
void dekodexb (int kode, int *aa, int len, int base)  ;
char *binary_string(int a, int len) ;
void printmat0(double *a, int m, int n) ;
void loadbaseprobs(double **basep, int *basevalid, PILEUP *pilept, int rind) ;
char *getfaname(char *sname) ;
int patkode(int a, int b, int c, int d) ;
char *kodepstring(int kode) ; 
int checkag(int *aa, int len) ;
int abpat(int *ww, int *aa, int len, int rind) ;
int setrdtable(char **plist, int nplist, char *rdname) ;
int loadbadbonelist(char *badlist)   ;
int  loadlist(char **list, char *listname)   ;
char *num2h(int k) ;
int abx(int a, int b) ;
int abxok(int abx, int abxmode) ; 
char *abxstring(int abx) ;

int isold(char *pop) ;
int isdenisova(char *pop) ;
int isneandertal(char *pop) ;
void setancient(READ *readpt) ;

int calcvar(PILEUP *pilept, int *vnum, int nean) ;
void cntalleles(PILEUP *pilept, int *ccc) ; 
void cntallelesnonean(PILEUP *pilept, int *ccc, int nean) ;
void forcebi(PILEUP *pilept, int xa, int xb) ; 
void copybb(BASEP *basept, int a, int b) ;
int calcvind(int *ccc, int rind, int *pvind) ;
int checkicx(int *ccc, int rind, int *pvind) ;
int checktri(int *ccc, int rind) ;

int getreadx(READ *readpt, char *line) ;
int getreadf(READ *readpt, FILE *fff) ;
int getreadfnorg(READ *readpt, FILE *fff) ;
int getread(READ *readpt) ;
int calcnumposdiff(READ *readpt) ;
void cleanread(READ *readpt) ;
void loadqinfo(char *qinfoname, double qnean) ;
void printqinfo() ;
int goodqualp(READ *readpt, int i, THRESH *pthresh)  ;
int goodqual(READ *readpt, int i)  ;
int is_mono(int *aa, int len)  ;
int is_biallelic(int *aa, int len)  ;
void mkqual(int *iqual, char *qual, int len) ; 
void getpophdr(char **pophdr, char **poplist, int npops) ;
int getfastalist(char **poplist, int npops, char *dbfile, char **rglist, int *rgpops)  ;
int getrg(char **poplist, int npops, char *dbfile, char **rglist, int *rgpops)  ;
int getrgj(char **rawpoplist, int npops, char *dbfile, char **rglist, int *rgpops) ;
void printrghit()  ;
int mkfullpoplist(char **fullpoplist, int *f2pops, char **rawpoplist, int npops) ;
int countg(char *uu, int nuu, int snpref, int snpvar)  ;
int countc(char *uu, int nuu, char var, char ref)  ;
void printread(READ *readpt) ;
int grabreference(char *reg, int pos, int len, char *refstring)  ;
void getnumreads(long *plong, long *ptotal)  ;
void makergcats(char *rrr, char *bamname)  ;

char *getfaiub(int k) ;
int isloshlass(char *pop) ;

#ifndef NJPQINFO
typedef struct {
 char popname [15] ; 
 int popind ;
 int chimpbasenum ;
 int strand ;
 int loqual ; 
 int hiqual ;
 double loprob ;
} QINFO ;

// specifies loqual and hiqual for pop given chimp base  
// refinement if qual = loqual base only accepted with probability loprob

#endif 
#define NJPQINFO

#ifndef _POPDEF_ 
#define _POPDEF_


#endif

