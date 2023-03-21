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
  INTSXP,  // 9 LYRF
  REALSXP, // 10 APRF
  INTSXP,  // 11 parentsRF
  REALSXP, // 12 LrRF
  INTSXP,  // 13 OhRF
  INTSXP,  // 14 Nd
  INTSXP,  // 15 DumParRF
  REALSXP, // 16 DumLrRF
  INTSXP,  // 17 DumBYRF
  REALSXP, // 19 TotLL
	REALSXP, // 19 AP_OUT
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
  INTSXP,  // 3 SpecsIntMkPed
  REALSXP, // 4 SpecsDbl
  REALSXP, // 5 ErrV
  INTSXP,  // 7 GenoFR
  INTSXP,  // 8 Sex
	INTSXP,  // 9 BY
	REALSXP, // 10 AP
	INTSXP,  // 11 parentsRF
	INTSXP,  // 12 OHRF
	REALSXP, // 13 LRRF
	INTSXP,  // 14 SnpdBoth
	INTSXP,  // 15 dumparrf
	REALSXP, // 16 dumLRrf
	INTSXP,  // 17 dumbyrf
};

static R_NativePrimitiveArgType pairLLType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 Np
  INTSXP,  // 3 SpecsInt
  REALSXP, // 4 SpecsDbl
  REALSXP, // 5 ErrV
  INTSXP,  // 6 GenoFR
  INTSXP,  // 7 BYRF
  REALSXP, // 8 AP
  INTSXP,  // 9 pairIDs
	INTSXP,  // 10 pairSex
  INTSXP,  // 11 pairAgeDiff
  INTSXP,  // 12 pairFocal
  INTSXP,  // 13 pairk
  INTSXP,  // 14 dropP
	INTSXP,  // 15 parentsRF
  INTSXP,  // 16 dumparRF
	REALSXP, // 17 LLRF
	INTSXP,  // 18 TopRF
	REALSXP, // 19 dLRF
};

static R_NativePrimitiveArgType BYprobType[] = {
  INTSXP,  // 1 Ng
  INTSXP,  // 2 Nx
  INTSXP,  // 3 nAP
  INTSXP,  // 4 nYearsIn
	INTSXP,  // 5 BY
	INTSXP,  // 6 LYRF
	REALSXP, // 7 AP
	INTSXP,  // 8 parentsRF
	REALSXP, // 9 byprobv
};

static R_NativePrimitiveArgType eType[] = {
  INTSXP,
  INTSXP,
  INTSXP,
  REALSXP,
};


static R_NativePrimitiveArgType relType[] = {
  INTSXP,
  INTSXP,
  INTSXP,
  INTSXP,
};


extern void F77_NAME(makeped)(int *ng, int *specsintglb, int *specsintmkped,
  double *specsdbl, double *errv, int *genofr, int *sexrf, int *byrf, int *lyrf,
	double *aprf, int *parentsrf, double *lrrf, int *ohrf,
	int *nd, int *dumparrf, double *dumlrrf, int *dumbyrf, double *totll, double *apout);

extern void F77_NAME(duplicates)(int *ng, int *specsint, double *specsdbl,
	double *errv, int *genofr, int *sexrf, int *byrf, double *aprf,
	int *ndupgenos, int *dupgenos, int *nmismatch, int *snpdboth, double *duplr);

extern void F77_NAME(findambig)(int *ng, int *specsint, int *specsintamb, double *specsdbl,
  double *errv, int *genofr, int *sexrf, int *byrf, double *aprf, int *parentsrf,
  int *dumparrf, int *namb, int *ambigid, int *ambigrel, double *ambiglr, int *ambigoh,
	int *ntrio, int *trioids, double *triolr, int *triooh);

extern void F77_NAME(getpedllr)(int *ng, int *specsint, int *specsintmkped,
  double *specsdbl, double *errv, int *genofr, int *sexrf, int *byrf,
  double *aprf, int *parentsrf, int *ohrf, double *lrrf, int *snpdboth, int *dumparrf,
  double *dumlrrf, int *dumbyrf);

extern void F77_NAME(getpairll)(int *ng, int *np, int *specsint, double *specsdbl,
  double *errv, int *genofr, int *byrf, double *aprf,
  int *pairids, int *pairsex, int *pairagediff, int *pairfocal, int *pairk,
  int *dropp, int *parentsrf, int *dumparrf, double *llrf, int *toprf, double *dlrf);

extern void F77_NAME(getbyprobs)(int *ng, int *nx, int *nap, int *nyearsin, int *byrf,
  int *lyrf, double *aprf, int *parentsrf, double *byprobv);

extern void F77_NAME(deallocall)(void);

extern void F77_NAME(mkerrors)(int *nind, int *nsnp, int *genofr, double *eprobfr);

extern void F77_NAME(getrel)(int *nind, int *pedrf, int *nrel, int *relv);


static const R_FortranMethodDef FortranEntries[] = {
	{"makeped", (DL_FUNC) &F77_NAME(makeped), 19, psType},
	{"duplicates", (DL_FUNC) &F77_NAME(duplicates), 13, dupType},
  {"findambig", (DL_FUNC) &F77_NAME(findambig), 20, ambigType},
	{"getpedllr", (DL_FUNC) &F77_NAME(getpedllr), 16, pedLLRType},
  {"getpairll", (DL_FUNC) &F77_NAME(getpairll), 19, pairLLType},
  {"getbyprobs", (DL_FUNC) &F77_NAME(getbyprobs), 9, BYprobType},
  {"deallocall", (DL_FUNC) &F77_NAME(deallocall), 0},
	{"mkerrors", (DL_FUNC) &F77_NAME(mkerrors), 4, eType},
  {"getrel", (DL_FUNC) &F77_NAME(getrel), 4, relType},
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
