int qcpg(char *snpid, int pos, char *alleles) 
{ 
  int  t ;
  char cref, c ;
  static int ncall = 0 ;
  int pos2 ;

  ++ncall ;

   cref = xrefbase(pos) ;
   if (cref != alleles[0]) printf("funny: %20s %c %c\n", snpid, alleles[0], cref) ;
   
   t = 0 ; 

   if ((alleles[0] == 'C') && (alleles[1] == 'T')) t=1 ;
   if ((alleles[0] == 'T') && (alleles[1] == 'C')) t=1 ;
   if ((alleles[0] == 'G') && (alleles[1] == 'A')) t=2 ;
   if ((alleles[0] == 'A') && (alleles[1] == 'G')) t=2 ;

   if (t==0) return NO ;
   pos2 = pos + 1 ;
   if (t==2) pos2 = pos-1 ;
   c = xrefbase(pos2) ;
   if ((t==1) && (c=='G')) return YES ;
   if ((t==2) && (c=='C')) return YES ;
   return NO ;
}


