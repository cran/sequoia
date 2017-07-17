#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

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
  INTSXP,  // 4 GenoFR
  INTSXP,  // 5 SexRF
  INTSXP,  // 6 BYRF
  REALSXP,  // 7 APRF
  INTSXP,  // 8 parentsRF
  REALSXP, // 9 LrRF
  INTSXP,  // 10 OhRF
  INTSXP,  // 11 nAmb
  INTSXP,  // 12 AmbigID
  INTSXP,  // 13 AmbigSex
  INTSXP,  // 14 AmbigAgeDif
  INTSXP,  // 15 AmbigRel
  REALSXP, // 16 AmbigLR
  INTSXP,  // 17 AmbigOH
  INTSXP,  // 18 Nd
  INTSXP,  // 19 DumParRF
  REALSXP, // 20 DumLrRF
  INTSXP,  // 21 DumBYRF
  INTSXP,  // 22 DumNoff
  INTSXP,  // 23 DumOff  
  REALSXP, // 24 TotLL
};

extern void F77_NAME(duplicates)(int *Ng, int *SpecsInt, int *GenoFR, 
  int *nDupGenos, int *DupGenosFR, int *nMisMFR);

extern void F77_NAME(makeped)(int *Ng, int *SpecsInt, float *SpecsDbl,
  int *GenoFR, int *SexRF, int *BYRF, int *APRF, int *parentsRF, 
  float *LrRF, int *OhRF, int *nAmb, int *AmbigID, int *AmbigSex, 
  int *AmbigAgeDif, int *AmbigRel, float *AmbigLR, int *AmbigOH,
  int *Nd, int *DumParRF, float *DumLrRF, int *DumBYRF, 
  int *DumNoff, int *DumOff,  float *TotLL);  

extern void F77_NAME(deallocall)();


static const R_FortranMethodDef FortranEntries[] = {
  {"duplicates", (DL_FUNC) &F77_NAME(duplicates), 6, dupType},
  {"makeped", (DL_FUNC) &F77_NAME(makeped), 24, psType},
  {"deallocall", (DL_FUNC) &F77_NAME(deallocall), 0},
  {NULL, NULL, 0}
};


void R_init_sequoia(DllInfo *info)
{
  R_registerRoutines(info,
                     NULL,          // .C
                     NULL,          // .Call
                     FortranEntries, // .Fortran
                     NULL);         // .External
  R_useDynamicSymbols(info, FALSE);
}
