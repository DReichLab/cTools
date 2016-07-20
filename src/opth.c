#include "opth.h" 

#define MAXSTR  512

int bbans[10] ;

extern int verbose ;
char *iname = "hist.out"  ;

int maxlistlen = 100*1000 ;
int maxdepth = 99 ; 
int minmq = 60 ;


void setfpars(FENTRY **felist, int felen,  int *bbans)  
{
  int *bestindex, j, k, t ; 
  int ansindex[10] ;
  int coverpc, maxcindex ; 
  double ycover, erate, olderate  ; 
  double y0, y1, yy ;
  FENTRY *fept ; 

  ZALLOC(bestindex, 1000, int) ;
  ivclear(bestindex, -1, 1000) ;
  for (coverpc  = 200; coverpc < 1000; ++coverpc) { 
   ycover = (double) coverpc / 1000.0 ;
   for (k=0; k<felen; ++k) {   
    fept = felist[k] ; 
    if (fept -> cover < ycover) continue ; 
    erate = fept -> erate ; 
    t = bestindex[coverpc] ;
    if (t<0) { 
     bestindex[coverpc] = k ;  
     continue ; 
    }
    olderate =  felist[t] -> erate ;
    if (erate < olderate) bestindex[coverpc] = k  ; 
   }
  }
  y0 = -2 ;
  for (k=0; k<felen; ++k) {   
   fept = felist[k] ; 
   if (fept -> cover > y0) { 
    maxcindex = k ;
    y0 = fept -> cover ;
   }
   y0 = MAX(y0, fept -> cover) ;
  }
  y1 = 0.5 ;
  ivclear(ansindex, maxcindex, 10) ;
  if (y0 >= y1 )  {  // actual filtering 
    for (k=1; k<=9; ++k) { 
     yy = ((double) (9-k)  * y0 + (double) k * y1 ) / 9.0 ;
     erate = 99 ; 
     for (j=0; j<felen; ++j) { 
      fept = felist[j] ; 
      if (fept -> cover < yy) continue ; 
      t = ansindex[k] ;  
      olderate = felist[t] -> erate ;
      erate = fept -> erate ;
      if (erate < olderate) ansindex[k] = j ; 
     }
   }
  }

  copyiarr(ansindex, bbans, 10) ;
  free(bestindex) ;

}
int mktlist(long **dd, int dlen, int *tlo, int *thi, long *thh, long *tmm) 
{
 long h, m, hh, mm ;
 long a, b, x, t, v ; 
 int k, lo, hi  ; 
 int *zlo, *zhi, *indx ;
 long *zhh, *zmm ; 
 double *cover, *erate ; 
 double yerate, besterate ;
 long total = 0 ; 
 long tlim = 1000 ; // must have this much
 int a1, a2, a3, b1, b2, b3 ; 
 double cov1, cov2, cov3, err1, err2, err3 ;
 double aa, ee2 ;

 t = dlen*dlen ; 
 ZALLOC(zlo, t, int) ;
 ZALLOC(zhi, t, int) ;
 ZALLOC(zhh, t, long int) ;
 ZALLOC(zmm, t, long int) ;
 ZALLOC(indx, t, int) ;



 x = 0 ;
 for (a=0; a<dlen; ++a) { 
  h = dd[a][0] ; 
  m = dd[a][1] ; 
  total += (h+m) ;
  t = (h+m) ; if (t==0) continue ; 
  hh = h ; mm = m ;
  t = hh + mm ; 
  if (t<tlim) continue ;
  zlo[x] = zhi[x] = a ; 
  zhh[x] = hh ;  zmm[x] = mm ;
  ++x ;
 for (b=a+1; b<dlen; ++b) { 
  h = dd[b][0] ; 
  m = dd[b][1] ; 
  t = (h+m) ; if (t==0) continue ; 
  hh += h ; mm += m ;
  if (hh<0) fatalx("badhh\n") ;
  t = hh + mm ; 
  if (t<tlim) continue ;
  zlo[x] = a ; 
  zhi[x] = b ; 
  zhh[x] = hh ;  zmm[x] = mm ;
  if (zhh[x]<0) fatalx("badzhh\n") ;
  ++x ;
 }}
// printf("zzhx %d %ld\n", x, total) ;
 ZALLOC(cover, x, double) ; 
 ZALLOC(erate, x, double) ; 
 for (k=0; k<x; ++k) { 
  t = zhh[k] + zmm[k] ; 
  cover[k] = (double) t / (double) total ;
  erate[k] = (double) zmm[k] / (double) t ;
//  printf("zzerate::%d  %d %d %ld %ld  %9.3f %9.3f\n", k, zlo[k], zhi[k], zhh[k], zmm[k], cover[k], erate[k]) ;
 }
 if (x>0) sortit(cover, indx, x) ; 
 besterate = 100 ; 
 v = 0 ;
 lo = -1; hi = 999 ;
 a1 = b1  = 0; 
 for (a=1; a<x-2; ++a) { 
  break ; 
  b2 = indx[a] ; 
  b3 = indx[a+1] ; 
  a2 = a ; a3 = a+1 ;
  cov1 = cover[a1] ;
  err1 = erate[b1] ;
  cov2 = cover[a2] ;
  err2 = erate[b2] ;
  cov3 = cover[a3] ;
  err3 = erate[b3] ;
  aa = (cov2-cov1)/(cov3-cov1) ; 
  ee2 = err1 + aa*(err3-err2) ; 
  if (ee2<err2) { 
   cover[a2] -1 ; continue ;
  }
  a1 = a2 ; 
  b1 = b2 ;
 }



 for (a=x-1; a>=0; --a) { 
  b = indx[a] ;  
  yerate = erate[b] ;
/**
  if (a==0) printf("zzww %d %d %d  %d %d %12.6f %12.6f %12.6f\n", 
 b, zlo[b], zhi[b], zhh[b], zmm[b], cover[a], erate[b], besterate) ;
*/
  if (yerate>=besterate) cover[a] = -1 ; 
  if (zlo[b] < lo) cover[a] = -1 ;
  if (zhi[b] > hi) cover[a] = -1 ;
  if (cover[a] >= 0) {
   lo = zlo[b] ; 
   hi = zhi[b] ;
   besterate = yerate ; 
  }
 }
 for (a=0; a<x; ++a) { 
  if (cover[a]<0) continue ;
  b = indx[a] ;  
  tlo[v] = zlo[b] ; 
  thi[v] = zhi[b] ; 
  thh[v] = zhh[b] ; 
  tmm[v] = zmm[b] ; 
  ++v ; 
 }
       
 free(cover) ; 
 free(erate) ; 

 free(zlo) ; 
 free(zhi) ; 
 free(zhh) ; 
 free(zmm) ; 
 free(indx) ; 

 return v ;

}

void clearfe(FENTRY *fept)  
{
 ivclear(fept -> lovals, 0, 10) ; 
 ivclear(fept -> hivals, 0, 10) ; 
 fept -> cover = 0 ; 
 fept -> erate = 1000 ;
 fept -> level = 0 ;
}

