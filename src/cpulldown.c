/*
* cpulldown.c:	get the genotypes of the given individuls at the given SNP loci
* Author: Nick Patterson
* Revised by: Mengyao Zhao
* Last revise date: 2014-11-26
* Contact: mengyao_zhao@hms.harvard.edu
*/

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <nicklib.h>
#include <globals.h>
#include <mcmcpars.h>
#include <getpars.h>
#include <nicksam.h>

#include "faidx.h"
#include "admutils.h"
#include "mcio.h"  

#define WVERSION   "131" 

// fai_destroy called
#define MAXFL  50   
#define MAXSTR  512

//char *iubfile = "/home/np29/cteam/release/hetfaplus.dblist" ;
//char *iubmaskfile = "/home/np29/cteam/release/maskplus.dblist" ;

char *iubfile = "/home/mz128/cteam/dblist/hetfa_postmigration.dblist" ;
char *iubmaskfile = "/home/mz128/cteam/dblist/mask_postmigration.dblist" ;

extern enum outputmodetype outputmode  ;
extern int checksizemode ;
char *omode = "eigenstrat" ;
extern int packmode ;            //!< flag - input {is not,is} in packed mode
extern char *packgenos ;         //!< packed genotype data (packit.h)
extern char *packepath ;
extern long packlen;             //!< allocated size of packgenos data space
extern long rlen;                //!< number of bytes in packgenos space that each SNP's data occupies

char *trashdir = "/var/tmp" ;
extern int verbose  ;
int qtmode = NO ;
Indiv **indivmarkers;
SNP **snpmarkers ;
int numsnps, numindivs ; 

char  *genotypename = NULL ;
char  *snpname = NULL ;

char *indoutfilename = NULL ;
char *snpoutfilename = NULL ;
char  *genooutfilename = NULL ;

char  *indivname = NULL ;

int minfilterval = 0 ;  // note default
int checkmode = NO ;

char **samplist ; 
int nsamps ;
int *hasmask ;

int minchrom = 1 ;
int maxchrom = 97 ;
int xchrom = -1 ;

char *outputname = NULL ;
FILE *ofile ;

double fakespacing = 0.0 ;
char  unknowngender = 'U' ;

char *regname ;
char **iublist ;
char **iubmask ;

char *fasta ;
char *mask ;

void readcommands(int argc, char **argv) ;
long setgenoblank (SNP **snpmarkers, int numsnps, int numindivs)   ;
int readfa1(char *faname, char **pfasta, int *flen) ;
char *myfai_fetch(faidx_t *fai, char *reg, int  *plen) ;
int  mksamplist(char **samplist, Indiv **indivmarkers, int numindivs) ;
int getfalist(char **poplist, int npops, char *dbfile, char **iublist)  ;
int checkr(char **samplist, int nsamps, char **iublist, char **iubmask) ;
char fixval(char iub, char cm) ; 
void getfasta(char **pfasta, char **pmask, int *rlen, int *mlen, int kk) ;

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   cpulldown -p <parameter file> [options] \n\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "\t-V	verbose\n");
	fprintf(stderr, "\t-c	checkmode\n");
	fprintf(stderr, "\t-v	Show version information.\n");
	fprintf(stderr, "\t-? 	Show the instruction. (For detailed instruction, please see the document here: https://github.com/mengyao/cTools)\n\n");
	return 1;
}

int main(int argc, char **argv)
{

  int i, j, k, kk, g, t, tt, ccheck ; 
  int pos, isnp, gval ;
  SNP *cupt ;
  Indiv *indx ;

  int numvind, nignore, numrisks = 1 ;
  int chrom, lastchrom = -1 ;  
  int firstsnp, rlen, mlen  ;
  char *iubarr, *maskarr  ;

  char ss[100], cbases[2], iub, cm, c1, c2, cval ;

  ofile = stdout; 
  packmode = YES ;
  readcommands(argc, argv) ;
  settersemode(YES) ;
  if (outputname != NULL) openit(outputname, &ofile, "w") ;

  setomode(&outputmode, omode) ;

  numsnps = 
    getsnps(snpname, &snpmarkers, fakespacing,  NULL, &nignore, numrisks) ;

// fakespacing 0.0 (default)

  numindivs = getindivs(indivname, &indivmarkers) ;
  ZALLOC(samplist, numindivs, char *) ;
  nsamps = mksamplist(samplist, indivmarkers, numindivs) ;

  ZALLOC(hasmask, nsamps, int) ;

  printf("numsnps: %d  numindivs:  %d\n", numsnps, numindivs) ;

  ZALLOC(iublist, nsamps, char *) ;
  ZALLOC(iubmask, nsamps, char *) ;
  getfalist(samplist,nsamps, iubfile, iublist)  ;
  getfalist(samplist,nsamps, iubmaskfile, iubmask)  ;

   for (k=0; k<nsamps; ++k) {
    hasmask[k] = YES ;
    if (iubmask[k] == NULL) {  
     hasmask[k] = NO ;
     continue ;
    }
    t = strcmp(iubmask[k], "NULL") ; 
    if (t==0) {
     hasmask[k] = NO ;
     continue ;
    }
   }


  checkr(samplist, nsamps, iublist, iubmask) ;

  for (k=0; k<nsamps; ++k) { 
   printf("%15s\n%s\n%s\n", samplist[k], iublist[k], iubmask[k]) ;
   printnl() ;
  }
  if (checkmode == YES) return 0 ;
  setgenoblank (snpmarkers,  numsnps,  numindivs)   ;

  ZALLOC(iubarr, numindivs+1, char) ;
  ZALLOC(maskarr, numindivs+1, char) ;
  maskarr[numindivs] = iubarr[numindivs] = CNULL ;
  for (isnp = 0 ; isnp < numsnps; ++isnp) { 
   cupt = snpmarkers[isnp]  ; 
   chrom = cupt -> chrom ;
   ccheck = YES ;
   if (chrom < minchrom) ccheck = NO  ;
   if (chrom > maxchrom) ccheck = NO  ;
   if ((xchrom>0) && (chrom != xchrom)) ccheck = NO  ;
   if (ccheck == NO)   { 
    cupt -> ignore = YES ;
    continue ;
   }
   if (chrom == lastchrom) continue ;
   lastchrom = chrom ;
   if (chrom < minchrom) continue ;
   if (chrom > maxchrom) continue ;
   if ((xchrom>0) && (chrom != xchrom)) continue ;
    sprintf(ss, "%d", chrom) ;
    freestring(&regname) ;
    regname = strdup(ss) ;
    firstsnp = isnp ;
    for (k=0; k<numindivs; ++k) { 
     indx = indivmarkers[k] ;
     kk = indxindex(samplist, nsamps, indx -> ID) ;
     getfasta(&fasta, &mask, &rlen, &mlen, kk) ;
// snp array must have ref allele on forward strans
     printf("fasta retrieved: %d %s %d %d\n", chrom, indx -> ID, rlen, mlen)  ;  fflush(stdout) ;
     if (rlen==0) continue ; 
     for (j=firstsnp; j<numsnps; ++j) { 
      cupt = snpmarkers[j] ;
      if (cupt -> chrom != chrom) break ;
      pos = cupt -> physpos ;
      t = pos-1 ;
      if (t>=rlen) continue ;
      gval = -1 ;
      iub = fasta[t] ;
      if (hasmask[kk]) cm = mask[t] ;
      else cm = '9' ;
      iub = toupper(iub) ;
      cval = fixval(iub, cm) ;
      tt = iubcbases(cbases,  cval) ;   
      if (tt<0) continue ;  // invalid default
      c1 = cupt -> alleles[0] ;
      c2 = cupt -> alleles[1] ;

      if ((cbases[0] != c1) && (cbases[0] != c2)) { 
       if (verbose) printf("triallelic: %s %s %c %c %c %c\n",  cupt -> ID, indx -> ID, c1, c2, cbases[0], cbases[1]) ;
       continue ;
      }
      if ((cbases[1] != c1) && (cbases[1] != c2)) { 
       if (verbose) printf("triallelic: %s %s %c %c %c %c\n",  cupt -> ID, indx -> ID, c1, c2, cbases[0], cbases[1]) ;
       continue ;
      }
      t = 0 ;
      if (cbases[0] == c1) ++t ;
      if (cbases[1] == c1) ++t ;
      putgtypes(cupt, k, t) ;
//    printf("zz %d %d %d %c %c %c\n", j, k, t, iub, cbases[0], cbases[1]) ;
     }
    
   }
  }

  outfiles(snpoutfilename, indoutfilename, genooutfilename, 
   snpmarkers, indivmarkers, numsnps, numindivs, -1, NO) ;


  printf("##end of run\n") ;
  return 0 ;
}
void readcommands(int argc, char **argv) 

{
  int i,haploid=0;
  char *parname = NULL ;
  phandle *ph ;
  char str[5000]  ;
  char *tempname ;
  int n ;

  while ((i = getopt (argc, argv, "p:cvV?")) != -1) {

    switch (i)
      {

      case 'p':
	parname = strdup(optarg) ;
	break;

      case 'v':
	printf("version: %s\n", WVERSION) ; 
	break; 

      case 'c':
	checkmode = YES ;
	break; 


      case 'V':
	verbose = YES ;
	break; 

      case '?':
	default:
	exit(usage());

//	printf ("Usage: bad params.... \n") ;
//	fatalx("bad params\n") ;
      }
  }

         
   if (parname == NULL) //return ;
		exit(usage());
  // pcheck(parname,'p') ;
   printf("parameter file: %s\n", parname) ;
   ph = openpars(parname) ;
   dostrsub(ph) ;

   getint(ph, "minfilterval:", &minfilterval) ;
   getstring(ph, "genotypename:", &genotypename) ;
   getstring(ph, "snpname:", &snpname) ;
   getstring(ph, "indivname:", &indivname) ;
 //  getstring(ph, "indoutfilename:", &indoutfilename) ;
   getstring(ph, "indivoutname:", &indoutfilename) ; /* changed 11/02/06 */
 //  getstring(ph, "snpoutfilename:", &snpoutfilename) ;
   getstring(ph, "snpoutname:", &snpoutfilename) ; /* changed 11/02/06 */
  // getstring(ph, "genooutfilename:", &genooutfilename) ; 
   getstring(ph, "genotypeoutname:", &genooutfilename) ; /* changed 11/02/06 */
   getstring(ph, "outputformat:", &omode) ;  
   getint(ph, "minchrom:", &minchrom) ;
   getint(ph, "maxchrom:", &maxchrom) ;
   getint(ph, "chrom:", &xchrom) ;
 //  getstring(ph, "dbhetfa:", &iubfile) ;
 //  getstring(ph, "dbmask:", &iubmaskfile) ;
   writepars(ph) ;
   closepars(ph) ;

}
long setgenoblank (SNP **snpmarkers, int numsnps, int numindivs)   
{
  double y ;
  char *pbuff ;
  long ngenos = 0 ; 
  int j, k ; 
  SNP *cupt ;
  // rlen is number of bytes needed to store each SNP's genotype data

  y = (double) (numindivs * 2) / (8 * (double) sizeof (char)) ;   
  rlen = nnint(ceil(y)) ;
  packlen = rlen*numsnps ;
  packmode = YES ;
  ZALLOC(packgenos, packlen, char) ;
  pbuff = packgenos ;
  clearepath(packgenos) ;

  for (j=0; j<numsnps; ++j) { 

         cupt = snpmarkers[j] ;
         ZALLOC(cupt -> gtypes, 1, int) ;
         cupt -> pbuff = pbuff ;  
         pbuff += rlen ;
        cupt -> ngtypes = numindivs ;
        for (k=0; k<numindivs; ++k) {
          putgtypes(cupt, k, -1) ;      // initialize all individuals to "missing data"
          ++ngenos ;
        }
    }
    return ngenos ;
}

int readfa1(char *faname, char **pfasta, int *flen) 
{

 faidx_t *fai ;
 int k, len, t ;
 char *ttfasta ;
 char ssreg[20] ;
 int ntry = 0, itry ;
 
 if (pfasta != NULL) freestring(pfasta) ;

 *flen = 0 ;
 *pfasta = NULL ;

  if (faname  == NULL) { 
   return  0;
  }
  fai = fai_load(faname) ;
  strcpy(ssreg, regname) ;
  t = strcmp(regname, "23") ; if (t==0) strcpy(ssreg, "X") ;
  t = strcmp(regname, "24") ; if (t==0) strcpy(ssreg, "Y") ;
  t = strcmp(regname, "90") ; if (t==0) strcpy(ssreg, "MT") ;
  ttfasta = myfai_fetch(fai, ssreg, &len) ;
  fai_destroy(fai) ; // close files
  *flen = len ;
  if (ttfasta == NULL) { 
   printf ("***warning: no hetfa for %s :: %s\n", faname, ssreg) ;
   return len ;
  }
  *pfasta = strdup(ttfasta) ;
  freestring(&ttfasta) ;
  return len ;

}

int getfalist(char **poplist, int npops, char *dbfile, char **iublist)  

{
 char line[MAXSTR+1] ;
 char *spt[MAXFF], *sx ;
 char c ;
 int nsplit ;
 int  tt, t, k, s, nx = 0  ;
 FILE *fff ;
 char *scolon ; 
  
  if (dbfile == NULL) return 0 ;
  openit(dbfile, &fff, "r") ;

  while (fgets(line, MAXSTR, fff) != NULL)  {
   nsplit = splitup(line, spt, MAXFF) ; 
   if (nsplit<1) continue ;
   sx = spt[0] ;
   if (sx[0] == '#') { 
    freeup(spt, nsplit) ;
    continue ;
   }
   t = indxstring(poplist, npops, spt[0]) ; 
   if (t<0) { 
    freeup(spt, nsplit) ; 
    continue ;
   }
    sx = spt[2] ;
    tt = strcmp(sx, "NULL") ; 
    if (tt = 0) {  
     iublist[t] = NULL ;
    }
    else iublist[t] = strdup(sx) ;
    freeup(spt, nsplit) ;
    ++nx ;
  }

   fclose(fff) ;
   return nx ;

}
char *myfai_fetch(faidx_t *fai, char *reg, int  *plen)
{
  char *treg, *s ;
  treg = strdup(reg) ;

  if (fai==NULL) fatalx("(my_fai_fetch): fai NULL\n") ;

  s = fai_fetch(fai, treg, plen) ;
  if (*plen > 0) {
    free(treg) ;
    return s ;
  }
  free(treg) ;
  return NULL ;
}

void clearfainfo(FATYPE *fapt, int mode)
{

 if (mode==1) {
   freestring(&fapt->faname) ;
   freestring(&fapt->alias) ;
   freestring(&fapt->famask) ;
   fapt -> fai = fapt -> faimask =  NULL ;
 }
 freestring(&fapt -> rstring) ;
 freestring(&fapt -> mstring) ;
 freestring(&fapt -> regname) ;
 fapt -> lopos = fapt -> hipos = -1 ;
 fapt -> len = 0 ;
 fapt -> mlen = 0 ;
 fapt -> rlen = 0 ;

}

void printfapt(FATYPE *fapt)
{
  int k, np = 0, mlen = fapt -> mlen   ;
  char cc ;

  printf("fapt: %s %s\n",  fapt->faname, fapt->alias) ;
  printf("%x %d %d n", fapt -> fai, fapt -> lopos, fapt -> hipos) ;
  printf("len: %d rlen: %d\n", fapt -> len, fapt -> rlen) ;
  printf("mask: %s %d\n", fapt -> famask, fapt -> mlen) ;
  printnl() ;

  fflush(stdout) ;
}
           
int fvalid(char cm) 
// is cm indicating valid? 
{
  int t ; 

  if (minfilterval < 0) return YES ;
  t = (int) cm - (int) '0' ;

  if (t<0) return NO;
  if (t>=10) return NO ;
  if (t<minfilterval) return NO ;
  return YES ; 
}


char fixval(char iub, char cm) 
{
  int t ;

  if (cm == CNULL) return iub ;
  if (cm == '-')  return 'N' ; 
  t =(int) (cm) - (int) '0'  ;
  if (t<minfilterval) return '?' ;
  return iub ;
}

int mkcnt(int *cnt, char *iub, char *mask, int npops, char *pc1, char *pc2)  
// number of copies of ref allele
{
  char c1, c2, cm, x1, x2, cbases[2] ; 
  int i, t, x, ok ;
  
  c1 = *pc1; c2 = *pc2 ;

  for (i=0; i<npops; ++i) { 
   iub[i] = fixval(iub[i], mask[i]) ;
   t = iubcbases(cbases, iub[i]) ;
   if (t<0) { 
     cnt[i] = 3  ; iub[i] = '?' ; continue ;
   }
   x1 = cbases[0] ; x2 = cbases[1] ; 
   if ((c1==c2)  && (x1 != c1)) c2 = x1 ;
   if ((c1==c2)  && (x2 != c1)) c2 = x2 ;
   if ((x1 != c1)  && (x1 != c2)) return -1 ; 
   if ((x2 != c1)  && (x2 != c2)) return -1 ; 
   if ((x1==c1) && (t==1)) { 
     cnt[i] = 2; 
     continue ;
   }
   if ((x1==c2) && (t==1)) { 
     cnt[i] = 0; 
     continue ;
   }
   if (t==1) fatalx("logic bug\n") ;
   cnt[i] = 1 ;
  }
  *pc1 = c1 ; *pc2 = c2 ;

  t = 0 ; 
  for (i=0; i< npops-2; ++i) { 
   if (cnt[i] < 3) ++t ;
  }
  return t ;
}
int  mksamplist(char **samplist, Indiv **indivmarkers, int numindivs) 
{
  int k, n=0 ;
  Indiv *indx ;

  for (k=0; k<numindivs; ++k) { 
    indx = indivmarkers[k] ;
    if (indx -> ignore) continue ;
    samplist[n] = strdup(indx -> ID) ;
    ++n ;
  }
  return n ;

}

int checkr(char **samplist, int nsamps, char **iublist, char **iubmask) 
// semi debug routine; check hetfa and mask files can be read
{
  int nfail=0, t, k, ret, x ;
  char *sname ;

  for (k=0; k<nsamps; ++k) {
   sname = samplist[k] ;
// printf("checkr: %d %s\n", k,sname) ;
   t = 0 ;
   if (ftest(iublist[k]) == NO) { 
    printf("%s hetfa fail to open\n", sname) ;
    ++t ;
   }
   if (iubmask[k] == NULL) { 
    printf("%s mask fail to open\n", sname) ;
    ++t ;
   }
   else {
    x = strcmp(iubmask[k], "NULL") ;
    if ((x != 0) && ftest(iubmask[k]) == NO) { 
     printf("%s mask fail to open\n", sname) ;
     ++t ;
    }
   }
   if (t==0) printf("%s OK to read\n", sname) ; 
   nfail += t ;
  }
  return nfail ;
}

void getfasta(char **pfasta, char **pmask, int *rlen, int *mlen, int kk) 
{
  char *faname, *famask ;
  int tmlen=0, len, t  ;

  faname = iublist[kk] ;
  famask = iubmask[kk] ;

  freestring(pfasta) ;
  freestring(pmask) ;

  readfa1(faname, pfasta, &len) ;
  if (famask == NULL) t = 0 ; 
  else t = strcmp(famask, "NULL") ;
  if (t != 0)  {
   readfa1(famask, pmask, &tmlen) ;
   len = MIN(len, tmlen) ;
  }
  else { 
   tmlen = len ;
   ZALLOC(*pmask, len, char) ;
   cclear(*pmask, '9', len) ;
  }
  *rlen = len ;
  *mlen = tmlen ;
 
}
