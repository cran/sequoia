#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static R_NativePrimitiveArgType dupType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 SpecsInt
  STRSXP,  // 3 GenoFileFR
  INTSXP,  // 4 GenoFR
  INTSXP,  // 5 nDupGenos
  INTSXP,  // 6 DupGenosFR
  INTSXP,  // 7 nMisMFR
};

static R_NativePrimitiveArgType psType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 SpecsInt
  REALSXP, // 3 SpecsDbl
  STRSXP,  // 4 GenoFileFR
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


void F77_NAME(duplicates)(int *nMax, int *nDupGenoID, int *nDupLhID, 
  int *nDupGenos, int *nSexless, int *DupGenoIdFR, int *DupLhIdFR, 
  int *DupGenosFR, int *SexlessFR);

void F77_NAME(makeped)(int *Ng, int *SpecsInt, float *SpecsDbl,
  char **GenoFileFR, int *GenoFR, int *SexRF, int *BYRF, int *parentsRF, 
  float *LrRF, int *OhRF, int *nAmb, int *AmbigID, int *AmbigSex, 
  int *AmbigAgeDif, int *AmbigRel, float *AmbigLR, int *AmbigOH,
  int *Nd, int *DumParRF, float *DumLrRF, int *DumBYRF, 
  int *DumNoff, int *DumOff,  float *TotLL);  

void F77_NAME(deallocall)();


static const R_FortranMethodDef fortranMethod[] = {
  {"duplicates", (DL_FUNC) &F77_NAME(duplicates), 7, dupType},
  {"makeped", (DL_FUNC) &F77_NAME(makeped), 25, psType},
  {"deallocall", (DL_FUNC) &F77_NAME(deallocall), 0},
  {NULL, NULL, 0}
};


void R_init_sequoia(DllInfo *info)
{

  R_registerRoutines(info,
                     NULL,          // .C
                     NULL,          // .Call
                     fortranMethod, // .Fortran
                     NULL);         // .External
}
