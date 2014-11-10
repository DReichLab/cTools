int iubdekode(char *a, char iub) 
// a should be 5 long 
{ 

  switch (iub)  { 

   case 'A':  
    strcpy(a,"A") ;
    break ;
   case 'C':  
    strcpy(a,"C") ;
    break ;
   case 'G':  
    strcpy(a,"G") ;
    break ;
   case 'T':  
    strcpy(a,"T") ;
    break ;
   case 'M':  
    strcpy(a,"AC") ;
    break ;
   case 'R':  
    strcpy(a,"AG") ;
    break ;
   case 'W':  
    strcpy(a,"AT") ;
    break ;
   case 'S':  
    strcpy(a,"CG") ;
    break ;
   case 'S':  
    strcpy(a,"CT") ;
    break ;
   case 'Y':  
    strcpy(a,"GY") ;
    break ;
   case 'K':  
    strcpy(a,"ACG") ;
    break ;
   case 'V':  
    strcpy(a,"ACT") ;
    break ;
   case 'H':  
    strcpy(a,"AGT") ;
    break ;
   case 'D':  
    strcpy(a,"CGT") ;
    break ;
   case 'X':  
    strcpy(a,"ACGT") ;
    break ;
   case 'N':  
    strcpy(a,"ACGT") ;
    break ;
   
   default:  
    a[0] = CNYLL ;
  }
  return strlen(a) ;
}

