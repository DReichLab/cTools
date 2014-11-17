#include <nicklib.h>

/** 
 Dawson's Integral 
 exp(-t*t) \int ( exp(x^2), x = 0..t) 
 loosely based on mcerror.for 
*/

double dawson(double t) 
{
        
        double  z1,  cs, cr, cl ;
        double z1sq, cer ;
        int k ;

        z1 =  fabs(t) ;
        if (z1 <= 1.0e-12) return 0 ; 
        z1sq = -t*t ;
        if (z1 < 4.5) {  
         cs = cr = z1 ;
         for (k= 1; k <= 120; ++k) {  
            cr *= z1sq/((double) k + 0.5) ;
            cs += cr ;
            if (fabs(cr/cs) < 1.0e-15) break ;
         }
         cer = cs ;
        }
        else {
         cl = 1/z1 ;
         cr = cl ;
         for (k=1; k<=13; ++k) {
          cr *= -((double) k-0.5) / z1sq ;
          cl += cr ;
          if (fabs(cr/cl) < 1.0e-15) break ;
         }
         cer = 0.5*cl ;
        }
        if (t<0) cer = -cer ; 
        return cer ;
}

