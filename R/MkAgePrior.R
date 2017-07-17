#
# see also AgeDiff_priors.R

# this is P(age|relationship) / P(Age),

# to be used in P(relationship|age) = P(age|relationship) * P(relationship) / P(age)

# where P(Age) is the emperical age distribution in the sample.


# setwd("E:/Manuscripts/Inbreeding depression/AllData/Sequoia")

#======================================================================
#======================================================================

#' Age priors
#'
#' Calculate age-based prior probabilities for various categories of pairwise
#' relatives.
#'
#' if UseParents = TRUE, Retrieve age distributions of maternal & paternal
#' parents, siblings and grandparents from assigned parents, to use as input
#' for sibship clustering and grandparent assignment.
#' If the lifehistory file indicates a single age class, \eqn{MS = PS = 1} and
#'  \eqn{MGM = PGF = MGF = UA = 0}.
#'
#' @param UseParents use the age distribution of assigned parents. Otherwise,
#'   equal probabilities across all age differences are assumed.
#' @param nAgeClasses number of age classes; age prior matrix will have
#'   nAgeClasses + 1 rows.
#' @param Parents  dataframe with scaffold pedigree of assigned parents;
#'   columns ID - Dam - Sire
#' @param LifeHistData dataframe with 3 columns:
#'  \itemize{
#'  \item{ID: }{max. 30 characters long,}
#'  \item{Sex: }{1 = females, 2 = males, other numbers = unkown,}
#'  \item{Birth Year: }{(or hatching year) Zero and negative numbers are
#'    interpreted as missing values.}}
#'
#' @return A matrix with the probability ratio of two individuals being a
#'   certain type of relative versus being a random draw, conditional on the
#'   (absolute) age difference between those individuals.
#'
#'   One row per absolute age difference (0 - nAgeClasses), and one column for
#'   each relationship type, with abbreviations:
#'      \item{MS}{Maternal siblings}
#'      \item{PS}{Paternal siblings}
#'      \item{MGM}{Maternal grandmother}
#'      \item{PGF}{Paternal grandfather}
#'      \item{MGF}{Maternal grandfathers and paternal grandmothers}
#'      \item{UA}{Avuncular}
#'      \item{M}{Mothers}
#'      \item{P}{Fathers}
#'
MakeAgeprior <- function(UseParents = FALSE,
                         nAgeClasses = 0,
                         Parents = NULL,
                         LifeHistData = NULL)
{
  if (UseParents & !is.null(Parents)) {
    PropAss <- with(Parents, sum(!is.na(dam)) + sum(!is.na(sire))) / (2*nrow(Parents))
  } else {
    PropAss <- 0
  }

  if (nAgeClasses > 1) {
    if (!UseParents | PropAss < 0.25) {   # rather arbitrary threshold.
        AP <- matrix(1, nAgeClasses, 8,
                     dimnames=list(1:nAgeClasses,
                            c("MS", "PS", "MGM", "PGF", "MGF", "UA", "M", "P")))
        AP[1, c("M", "P")] <- 0
        AP[1:2, c("MGM", "PGF", "MGF")] <- 0
    } else {
      AP <- MkAgePrior(LHData = LifeHistData,
                       Par = Parents,
                       MaxT = nAgeClasses-1)
    }
  } else {
    AP <- as.matrix(t(c(MS = 1, PS = 1, MGM = 0, PGF = 0, MGF = 0, UA = 0,
                        M=0, P=0)))
  }
  AP
}


#=====================================================================
#=====================================================================

MkAgePrior <- function(LHData = NULL,
                       Par = NULL,
                       MaxT = NULL) {
names(LHData) <- c("ID", "Sex", "BY")
names(Par) <- c("ID", "Mother", "Father")

Ped <- merge(Par, LHData[,c("ID","BY")], by.x="ID", by.y="ID", all.x=T)
Ped <- merge(Ped, LHData[,c("ID","BY")],
                  by.x="Mother", by.y="ID", suffixes=c("",".mum"),all.x=T)
Ped <- merge(Ped, LHData[,c("ID","BY")],
                  by.x="Father", by.y="ID", suffixes=c("",".dad"),all.x=T)


################

BYdifs <- function(x) {
  dM <- abs(outer(x, x, "-"))
  dM[upper.tri(dM, diag = TRUE)] <- NA
  c(stats::na.exclude(c(dM)))
}

tbl <- function(x, n) table(factor(x, levels = 0:n))


################
# all

BYdiff.all <- BYdifs(Ped$BY)

if (is.null(MaxT)) MaxT <- diff(range(LHData$BY[LHData$BY>=0], na.rm=TRUE))
x <- c(0:MaxT)  # age difference

A.all <- tbl(BYdiff.all, n=MaxT)
P.all <- A.all/sum(A.all)

################
# maternal age

MumAge <- with(Ped, BY-BY.mum)
A.mum <- tbl(MumAge, n=MaxT)
LR.mum <- (A.mum/sum(A.mum))/P.all


################
# paternal age

DadAge <- with(Ped, BY-BY.dad)
A.dad <- tbl(DadAge, n=MaxT)
LR.dad <- (A.dad/sum(A.dad))/P.all


####################
# maternal sibships

tblMum <- table(Ped$Mother)
multiMum <- names(tblMum[tblMum>1])

MSA.L <- plyr::dlply(Ped[which(Ped$Mother %in% multiMum & !is.na(Ped$BY)),],
               "Mother", function(x) x$BY)
nMS <- sapply(MSA.L, length)

MSAD <- unlist(lapply(MSA.L, BYdifs))
A.MS <- tbl(MSAD, n=MaxT)
LR.MS <- (A.MS/sum(A.MS))/P.all



######################
# paternal sibships

tblDad <- table(Ped$Father)
multiDad <- names(tblDad[tblDad>1])

PSA.L <- plyr::dlply(Ped[which(Ped$Father %in% multiDad & !is.na(Ped$BY)),],
               "Father", function(x) x$BY)
PSAD <- unlist(lapply(PSA.L, BYdifs))
A.PS <- tbl(PSAD, n=MaxT)
LR.PS <- (A.PS/sum(A.PS))/P.all

################
### age of grandparents ###

Age.M <- outer(x,x,"+")

AgeL <- list()
for (i in 0:MaxT) {
  AgeL[[i+1]] <- which(Age.M==i,arr.ind=T)
}

# paternal grandfather
PGF.pr.M <- outer(A.dad/sum(A.dad), A.dad/sum(A.dad), "*")
PGF.pr <- numeric(MaxT)
for (i in 0:MaxT) {
  PGF.pr[i+1] <- sum(PGF.pr.M[AgeL[[i+1]]])  # sum over age combinations probabilities
}

# maternal grandmother
MGM.pr.M <- outer(A.mum/sum(A.mum), A.mum/sum(A.mum), "*")
MGM.pr <- numeric(MaxT)
for (i in 0:MaxT) {
  MGM.pr[i+1] <- sum(MGM.pr.M[AgeL[[i+1]]])
}

# maternal grandfather, paternal grandmother
MGF.pr.M <- outer(A.dad/sum(A.dad), A.mum/sum(A.mum), "*")
MGF.pr <- numeric(MaxT)
for (i in 0:MaxT) {
  MGF.pr[i+1] <- sum(MGF.pr.M[AgeL[[i+1]]])
}



################
### age of aunts/uncles ###
# for simplicity, assume equal prob. mat-mat, mat-pat, pat-mat & pat-pat

AgeD.M <- outer(Age.M,x,"-")

AgeDL <- list()
for (i in 0:MaxT) {
  AgeDL[[i+1]] <- which(abs(AgeD.M)==i,arr.ind=T)
}


# pat sib of father
UA.pr.4D <- array(dim=c(rep(MaxT+1,3),4))
UA.pr.4D[,,,1] <- outer(PGF.pr.M, A.dad/sum(A.dad), "*")
UA.pr.4D[,,,2] <- outer(MGF.pr.M, A.dad/sum(A.dad), "*")
UA.pr.4D[,,,3] <- outer(MGF.pr.M, A.mum/sum(A.mum), "*")
UA.pr.4D[,,,4] <- outer(MGM.pr.M, A.mum/sum(A.mum), "*")
UA.pr.3D <- apply(UA.pr.4D, 1:3, mean)

UA.pr <- numeric(MaxT)
for (i in 0:MaxT) {
  UA.pr[i+1] <- sum(UA.pr.3D[AgeDL[[i+1]]])  # sum over age combo probs
}

################
# combine into 1 matrix & write to file

#===============
AgePrior <- cbind(MS = LR.MS,
                  PS = LR.PS,
                  MGM = c(MGM.pr/P.all),
                  PGF = c(PGF.pr/P.all),
                  MGF = c(MGF.pr/P.all),
                  UA = c(UA.pr/P.all),
                  M = LR.mum,
                  P = LR.dad)
AgePrior <- round(AgePrior, 3)
AgePrior[is.na(AgePrior)] <- 0    # age combination not obs between genotyped indiv.


# elongate the tails where AgePrior==0
MinP <- 0.001   # depends on rounding
for (x in 1:ncol(AgePrior)) {
  for (i in 1:(MaxT/2)) {
    if (x==1) break
    AgePrior[i,x] <- ifelse(AgePrior[i,x]>0, AgePrior[i,x],
                        ifelse(AgePrior[i+1,x]>MinP, -99, AgePrior[i,x]))
  }
  if (MaxT>2) {
    for (i in 3:MaxT) {
      AgePrior[i,x] <- ifelse(AgePrior[i,x]>0, AgePrior[i,x],
                         ifelse(AgePrior[i-1,x]>=MinP | AgePrior[i-2,x]>MinP, -99,
                            AgePrior[i,x]))
    }
  }
}
AgePrior[AgePrior==-99] <- MinP
AgePrior[1, c("M", "P")] <- 0
AgePrior[1:2, c("MGM", "PGF", "MGF")] <- 0

AgePrior
}

#=====================================================================
#=====================================================================



