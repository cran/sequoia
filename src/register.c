#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h> 


static R_NativePrimitiveArgType dupType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 SpecsInt
//  STRSXP,  //  GenoFileFR
  INTSXP,  // 3 GenoFR
  INTSXP,  // 4 nDupGenos
  INTSXP,  // 5 DupGenosFR
  INTSXP,  // 6 nMisMFR
};

static R_NativePrimitiveArgType psType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 SpecsInt
  REALSXP, // 3 SpecsDbl
//  STRSXP,  // 4 GenoFileFR  char **GenoFileFR,
  INTSXP,  // 5 GenoFR
  INTSXP,  // 6 SexRF
  INTSXP,  // 7 BYRF
  REALSXP,  // 8 APRF
  INTSXP,  // 9 parentsRF
  REALSXP, // 10 LrRF
  INTSXP,  // 11 OhRF
  INTSXP,  // 12 nAmb
  INTSXP,  // 13 AmbigID
  INTSXP,  // 14 AmbigSex
  INTSXP,  // 15 AmbigAgeDif
  INTSXP,  // 16 AmbigRel
  REALSXP, // 17 AmbigLR
  INTSXP,  // 18 AmbigOH
  INTSXP,  // 19 Nd
  INTSXP,  // 20 DumParRF
  REALSXP, // 21 DumLrRF
  INTSXP,  // 22 DumBYRF
  INTSXP,  // 23 DumNoff
  INTSXP,  // 24 DumOff  
  REALSXP, // 25 TotLL
};


void F77_NAME(duplicates)(int *Ng, int *SpecsInt, int *GenoFR, 
  int *nDupGenos, int *DupGenosFR, int *nMisMFR);

void F77_NAME(makeped)(int *Ng, int *SpecsInt, float *SpecsDbl,
  int *GenoFR, int *SexRF, int *BYRF, int *APRF, int *parentsRF, 
  float *LrRF, int *OhRF, int *nAmb, int *AmbigID, int *AmbigSex, 
  int *AmbigAgeDif, int *AmbigRel, float *AmbigLR, int *AmbigOH,
  int *Nd, int *DumParRF, float *DumLrRF, int *DumBYRF, 
  int *DumNoff, int *DumOff,  float *TotLL);  

void F77_NAME(deallocall)();


static const R_FortranMethodDef fortranMethod[] = {
  {"duplicates", (DL_FUNC) &F77_NAME(duplicates), 6, dupType},
  {"makeped", (DL_FUNC) &F77_NAME(makeped), 24, psType},
  {"deallocall", (DL_FUNC) &F77_NAME(deallocall), 0},
  {NULL, NULL, 0}
};


void attribute_visible R_init_sequoia(DllInfo *info)
{
  R_registerRoutines(info,
                     NULL,          // .C
                     NULL,          // .Call
                     fortranMethod, // .Fortran
                     NULL);         // .External
  R_useDynamicSymbols(info, FALSE);
}
