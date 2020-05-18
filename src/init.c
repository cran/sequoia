#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

static R_NativePrimitiveArgType psType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 SpecsIntGlb
  INTSXP,  // 3 SpecsIntMkPed
  REALSXP, // 4 SpecsDbl
  REALSXP, // 5 ErrV
  INTSXP,  // 6 GenoFR
  INTSXP,  // 7 SexRF
  INTSXP,  // 8 BYRF
  REALSXP, // 9 APRF
  INTSXP,  // 10 parentsRF
  REALSXP, // 11 LrRF
  INTSXP,  // 12 OhRF
  INTSXP,  // 13 Nd
  INTSXP,  // 14 DumParRF
  REALSXP, // 15 DumLrRF
  INTSXP,  // 16 DumBYRF
  REALSXP, // 17 TotLL
	REALSXP, // 18 AP_OUT
};

static R_NativePrimitiveArgType dupType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 SpecsInt
  REALSXP, // 3 SpecsDbl
  REALSXP, // 4 ErrV
  INTSXP,  // 5 GenoFR
  INTSXP,  // 6 SexRF
  INTSXP,  // 7 BYRF
  REALSXP, // 8 APRF
  INTSXP,  // 9 nDupGenos
  INTSXP,  // 10 DupGenos
  INTSXP,  // 11 nMisMatch
  INTSXP,  // 12 SnpdBoth
  REALSXP, // 13 DupLR
};

static R_NativePrimitiveArgType ambigType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 SpecsInt
  INTSXP,  // 3 SpecsIntAmb
  REALSXP, // 4 SpecsDbl
  REALSXP, // 5 ErrV
  INTSXP,  // 6 GenoFR
  INTSXP,  // 7 SexRF
  INTSXP,  // 8 BYRF
  REALSXP, // 9 APRF
  INTSXP,  // 10 parentsRF
  INTSXP,  // 11 DumParRF
  INTSXP,  // 12 nAmb
  INTSXP,  // 13 AmbigID
  INTSXP,  // 14 AmbigRel
  REALSXP, // 15 AmbigLR
  INTSXP,  // 16 AmbigOH
  INTSXP,  // 17 nTrio
  INTSXP,  // 18 trioID
  REALSXP, // 19 trioLR
  INTSXP,  // 20 trioOH
};


static R_NativePrimitiveArgType pedLLRType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 SpecsInt
  REALSXP, // 3 SpecsDbl
  REALSXP, // 4 ErrV
  INTSXP,  // 5 CalcLLR
  INTSXP,  // 6 GenoFR
  INTSXP,  // 7 Sex
	INTSXP,  // 8 BY
	REALSXP, // 9 AP
	INTSXP,  // 10 parentsRF
	INTSXP,  // 11 OHRF
	REALSXP, // 12 LRRF
	INTSXP,  // 13 SnpdBoth
	INTSXP,  // 14 dumparrf
	REALSXP, // 15 dumLRrf
	INTSXP,  // 16 dumbyrf
};


static R_NativePrimitiveArgType eType[] = {
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
};


extern void F77_NAME(makeped)(int *ng, int *specsintglb, int *specsintmkped, 
  double *specsdbl, double *errv, int *genofr, int *sexrf, int *byrf, 
	double *aprf, int *parentsrf, double *lrrf, int *ohrf, 
	int *nd, int *dumparrf, double *dumlrrf, int *dumbyrf, double *totll, double *apout);
	
extern void F77_NAME(duplicates)(int *ng, int *specsint, double *specsdbl, 
	double *errv, int *genofr, int *sexrf, int *byrf, double *aprf, 
	int *ndupgenos, int *dupgenos, int *nmismatch, int *snpdboth, double *duplr);

extern void F77_NAME(findambig)(int *ng, int *specsint, int *specsintamb, double *specsdbl,
  double *errv, int *genofr, int *sexrf, int *byrf, double *aprf, int *parentsrf, 
  int *dumparrf, int *namb, int *ambigid, int *ambigrel, double *ambiglr, int *ambigoh, 
	int *ntrio, int *trioids, double *triolr, int *triooh);

extern void F77_NAME(getpedllr)(int *ng, int *specsint, double *specsdbl,
  double *errv, int *calcllr, int *genofr, int *sexrf, int *byrf, double *aprf, int *parentsrf, 
	int *ohrf, double *lrrf, int *snpdboth, int *dumparrf, double *dumlrrf, int *dumbyrf);
	
extern void F77_NAME(deallocall)();

extern void F77_NAME(mkerrors)(int *nind, int *nsnp, int *genofr, double *eprobfr);



static const R_FortranMethodDef FortranEntries[] = {
	{"makeped", (DL_FUNC) &F77_NAME(makeped), 18, psType},
	{"duplicates", (DL_FUNC) &F77_NAME(duplicates), 13, dupType},
  {"findambig", (DL_FUNC) &F77_NAME(findambig), 20, ambigType},
	{"getpedllr", (DL_FUNC) &F77_NAME(getpedllr), 16, pedLLRType},
  {"deallocall", (DL_FUNC) &F77_NAME(deallocall), 0},
	{"mkerrors", (DL_FUNC) &F77_NAME(mkerrors), 4, eType},
  {NULL, NULL, 0, NULL}
};


void attribute_visible R_init_sequoia(DllInfo *info)  // attribute_visible -> error
{
  R_registerRoutines(info,
                     NULL,          // .C
                     NULL,          // .Call
                     FortranEntries, // .Fortran
                     NULL);         // .External
  R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);  //available from R 3.0.0
}
