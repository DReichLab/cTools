/** 
 Appropriate settings for Reich Lab.  
 Users must either specify hetfa and mask files in parameter file, or modify this include file and recompile
*/

#define CTOOLDIR  "/n/data1/hms/genetics/reich/1000Genomes/cteam_remap/B-cteam_lite/v0.2/C-FullyPublic__SignedLetterNoDelay__SignedLetterDelay/"  

char *version = "1100" ;
char *iubfile =  CTOOLDIR"/C.hetfa.dblist"   ;
char *iubmaskfile =  CTOOLDIR"/C.mask.dblist"   ;

