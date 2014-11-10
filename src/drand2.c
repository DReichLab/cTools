
double drand2()  
{
  double x, y ;
  double maxran, maxran1 ;
  static double eps = -1.0 ;
/** 
 DRAND is quantized 1/2^31 
 call it twice and get max precision 
*/

  if (eps < 0.0) {
   maxran = 1.0-DBL_EPSILON  ;
   maxran1 = (double) (BIGINT-1) / (double) BIGINT ;
   eps = maxran - maxran1 ;
  }

  x = DRAND() ;
  y = DRAND() ;
  return x + y * eps ;
}

