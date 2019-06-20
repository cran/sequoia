#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

static R_NativePrimitiveArgType psType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 SpecsInt
  REALSXP, // 3 SpecsDbl
  INTSXP,  // 4 GenoFR
  INTSXP,  // 5 SexRF
  INTSXP,  // 6 BYRF
  REALSXP, // 7 APRF
  INTSXP,  // 8 parentsRF
  REALSXP, // 9 LrRF
  INTSXP,  // 10 OhRF
  INTSXP,  // 11 Nd
  INTSXP,  // 12 DumParRF
  REALSXP, // 13 DumLrRF
  INTSXP,  // 14 DumBYRF
  REALSXP, // 15 TotLL
	INTSXP,  // 16 nDupGenos
  INTSXP,  // 17 DupGenosFR
  INTSXP,  // 18 nMisMFR
  REALSXP, // 19 DupGenoLR
};

static R_NativePrimitiveArgType ambigType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 SpecsInt
  REALSXP, // 3 SpecsDbl
  INTSXP,  // 4 GenoFR
	INTSXP,  // 5 SexRF
  INTSXP,  // 6 BYRF
  REALSXP, // 7 APRF
  INTSXP,  // 8 parentsRF
//	INTSXP,  // 9 Nd      
  INTSXP,  // 10 DumParRF
	INTSXP,  // 11 nAmb
  INTSXP,  // 12 AmbigID
  INTSXP,  // 13 AmbigRel
  REALSXP, // 14 AmbigLR
  INTSXP,  // 15 AmbigOH
	INTSXP,  // 16 nTrio
  INTSXP,  // 17 trioID
  REALSXP, // 18 trioLR
	INTSXP,  // 19 trioOH
};

static R_NativePrimitiveArgType eType[] = {
	INTSXP,  
  INTSXP, 
	INTSXP,
  REALSXP,
};

extern void F77_NAME(makeped)(int *Ng, int *SpecsInt, float *SpecsDbl,
  int *GenoFR, int *SexRF, int *BYRF, float *APRF, int *parentsRF, 
	float *LrRF, int *OhRF, int *Nd, int *DumParRF, float *DumLrRF, 
	int *DumBYRF, float *TotLL,
	int *nDupGenos, int *DupGenosFR, int *nMisMFR, float *DupGenoLR);  

extern void F77_NAME(findambig)(int *Ng, int *SpecsInt, float *SpecsDbl,
  int *GenoFR, int *SexRF, int *BYRF, float *APRF, int *parentsRF, 
  int *DumParRF, int *nAmb, int *AmbigID, int *AmbigRel, float *AmbigLR, 
	int *AmbigOH, int *nTrio, int *trioID, float *trioLR, int *trioOH);

extern void F77_NAME(deallocall)();

extern void F77_NAME(mkerrors)(int *Nind, int *nSnp, int *GenoFR,
	float *EProbFR);

static const R_FortranMethodDef FortranEntries[] = {
	{"makeped", (DL_FUNC) &F77_NAME(makeped), 19, psType},
  {"findambig", (DL_FUNC) &F77_NAME(findambig), 18, ambigType},
  {"deallocall", (DL_FUNC) &F77_NAME(deallocall), 0},
	{"mkerrors", (DL_FUNC) &F77_NAME(mkerrors), 4, eType},
  {NULL, NULL, 0, NULL}
};


void R_init_sequoia(DllInfo *info)  // attribute_visible -> error
{
  R_registerRoutines(info,
                     NULL,          // .C
                     NULL,          // .Call
                     FortranEntries, // .Fortran
                     NULL);         // .External
  R_useDynamicSymbols(info, FALSE);
//	R_forceSymbols(info, TRUE);  available from R 3.0.0 
}
