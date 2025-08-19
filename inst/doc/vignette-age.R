## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
library(sequoia)
knitr::opts_chunk$set(echo = TRUE, 
                      fig.height=4, fig.width=6, 
                      fig.pos="ht!")

## ----DeerAP, echo=FALSE, fig.cap="Example of 'ageprior' distribution in a red deer population", out.width="90%"----
knitr::include_graphics("Deer-AgePriorExtra.png")   

## ----pipeline, echo=FALSE, fig.cap="Pipeline overview", out.width="90%"-------
knitr::include_graphics("agepriors-pipeline.png")

## ----Glossary, echo=FALSE, results="asis"-------------------------------------
Glossary <- matrix(c("$BY_i$", "Birth year", "Time unit (year, decade, month) of $i$'s birth/ hatching",
         "$A_{i,j}$", "Age difference",  "$BY_i - BY_j$, i.e. if $j$ is older than $i$, than $A_{i,j}$ is positive",
         "$R_{i,j}$", "Relationship", "Relationship between $i$ and $j$, e.g. parent-offspring",
         "$\\alpha_{A,R}$", "Ageprior", "Probability ratio of relationship $R$ given $A$, versus $R$ for a random pair"),
       ncol=3, byrow=TRUE)

knitr::kable(Glossary, col.names = c("Symbol", "Term", "Definition"),
            caption="Terms and abbreviations", escape=FALSE, booktabs=TRUE)

## ----AP-grif1, out.width="70%"------------------------------------------------
MakeAgePrior(Pedigree = Ped_griffin, Smooth=FALSE, Flatten=FALSE)  # Note: Ped_griffin includes a birth year column

## ----RelAbbr, echo=FALSE------------------------------------------------------
as.DF <- function(V, cn=c("name", "x")) {  # named vector --> data.frame
  setNames(data.frame(c1 = names(V), c2=V), cn)
}

Rels <- c(M = "Mother - offspring",
          P = "Father - offspring",
          FS = "Full siblings",   # From version XX
          MS = "Maternal siblings (full + half)",
          PS = "Paternal siblings (full + half)",
          MGM = "Maternal grandmother",
          PGF = "Paternal grandfather",
          MGF = "Maternal grandfather (+ paternal grandmother)",
          UA = "avuncular (niece/nephew -- aunt/uncle, full + half)",
          "(M/P)(M/P/F)A" = "Avuncular; Mother's/Father's Maternal/Paternal/Full sibling")

Rels <- cbind(as.DF(Rels, c("Column (R)", "Meaning")),
              "Version" = c(rep("all", 2), "from 1.0", "all", "all",
                            rep("up to 1.3 *", 3), "up to 1.3",
                            "from 2.0 *"))

knitr::kable(Rels, row.names=FALSE, booktabs=TRUE,
             caption = "AgePriors column names per sequoia version. *: From sequoia v2.0 only in AgePriorExtra")

## ----AP-grif2-----------------------------------------------------------------
# sequoia() output for the griffin data is included in the package:
# see ?SeqOUT_griffin on how it was generated
round(SeqOUT_griffin$AgePriorExtra, 2)
PlotAgePrior(SeqOUT_griffin$AgePriorExtra)

## ----AP-default, fig.cap="Default minimally informative ageprior", out.width = "70%"----
MakeAgePrior(LifeHistData = LH_HSg5)

## ----AP-default-2-------------------------------------------------------------
table(LH_HSg5$BirthYear)

## ----APdiscrete, fig.cap="Discrete generations"-------------------------------
# running sequoia to also get the extra ageprior columns: 
Seq_HSg5 <- sequoia(SimGeno_example, LH_HSg5, Module="ped",
                    args.AP=list(Discrete = TRUE),
                    CalcLLR=FALSE, Plot=FALSE, quiet=TRUE)
PlotAgePrior(Seq_HSg5$AgePriorExtra) 

## ----APDiscrete2--------------------------------------------------------------
MakeAgePrior(Discrete = TRUE, MaxAgeParent = 2, Plot=FALSE)

## ----pedHead------------------------------------------------------------------
tail(Ped_griffin)      # the pedigree

## ----APall, fig.cap="Ageprior for overlapping generations (griffin example)", fig.pos="h!", out.width = "70%"----
AP.griffin <- MakeAgePrior(Ped_griffin, Smooth=FALSE, Flatten=FALSE, 
                           Return="all")   
names(AP.griffin)

## ----BYrange------------------------------------------------------------------
AP.griffin$BirthYearRange

## ----MAP----------------------------------------------------------------------
AP.griffin$MaxAgeParent

## ----IncrMAP, out.width = "70%"-----------------------------------------------
# Specify maximum age for fathers, but estimate max mother age from pedigree 
APX <- MakeAgePrior(Ped_griffin, Smooth=FALSE, MaxAgeParent = c(NA, 5))

## ----tblAR--------------------------------------------------------------------
AP.griffin[["tblA.R"]]

## ----tblAR-2------------------------------------------------------------------
AgeDifM <- outer(Ped_griffin$birthyear, Ped_griffin$birthyear, "-")
# closest relationship for each pair:
RelM <- GetRelM(Ped_griffin, patmat=TRUE, GenBack=1)
table(c(AgeDifM), c(RelM))  # vectorising first speeds up table() a lot

## ----tblAR-2a-----------------------------------------------------------------
# - 'MaxT' ensures sufficient, constant number of rows across intermediary tables
MaxAgeParent <- c(NA, NA)
MaxAgePO <- ifelse(!is.na(MaxAgeParent), MaxAgeParent, diff(AP.griffin$BirthYearRange)+1)  # in actual code with thorough input checks
MaxT <- max(MaxAgePO+1, diff(AP.griffin$BirthYearRange))

# - factor() ensures levels with 0 count are included
# - drop negative age differences: all pairs are included 2x
tblA.R <- table(factor(AgeDifM, levels = 0:MaxT),
                factor(RelM, levels = c("M", "P", "FS", "MHS", "PHS")))

## ----tblAR-2b-----------------------------------------------------------------
RelA <- GetRelM(Ped_griffin, patmat=TRUE, GenBack=1, Return="Array")
tblA.R <- sapply(c("M", "P", "FS", "MHS", "PHS"),
                 function(r) table(factor(AgeDifM[RelA[,,r]==1], 
                                          levels = 0:MaxT)))
head(tblA.R)

## ----tblAR-2c-----------------------------------------------------------------
tblA.R <- cbind(tblA.R,
                  "MS" = tblA.R[,"MHS"] + tblA.R[,"FS"],
                  "PS" = tblA.R[,"PHS"] + tblA.R[,"FS"])
tblA.R <- tblA.R[, c("M", "P", "FS", "MS", "PS")]  # drop MHS & PHS columns


# add column for age diff. across all relationships (except 'Self' pairs)
tblA.R <- cbind(tblA.R,
                "X" = table(factor(AgeDifM[c(RelM)!="S"], levels = 0:MaxT))) 

# All pairs are included 2x in AgeDifM/RelM/RelA, to fix this:
# - drop the negative age differences  (done with factor() above)
# - divide counts for row A=0 by 2  (including 'X' column)
tblA.R["0", ] <- tblA.R["0", ] / 2  

# check that this is really the output from MakeAgePrior():
all(tblA.R == AP.griffin[["tblA.R"]])

## ----calcPAR------------------------------------------------------------------
NAK.R <- colSums(tblA.R, na.rm=TRUE)  # total count per relationship
NAK.R
PA.R <- sweep(tblA.R, 2, STATS=NAK.R, FUN="/")
head(round(PA.R, 2))

## ----FSuseHS------------------------------------------------------------------
MinPairs.AgeKnown <- 20   # not user-setable.
# Use the product of the distributions if there are many MS & PS,
# and the average across the two distributions if there are a medium number:
if (NAK.R["FS"] / min(NAK.R[c("MS", "PS")]) < 0.5 & 
    all(NAK.R[c("MS", "PS")] > MinPairs.AgeKnown)) {
  if (all(NAK.R[c("MS", "PS")] > 5*MinPairs.AgeKnown)) {
    FS.tmp <- PA.R[,"MS"] * PA.R[,"PS"]
    FS.tmp <- FS.tmp/sum(FS.tmp)
  } else {
    FS.tmp <- apply(PA.R[, c("MS", "PS")], 1, mean)
  }
  PA.R[,"FS"] <- (PA.R[,"FS"] + FS.tmp)/2
}

round(PA.R[1:5, ], 2)  # note change in 'FS' column

## ----calcLRAR-----------------------------------------------------------------
LR.RU.A <- sweep(PA.R[, 1:5], 1, STATS=PA.R[,"X"], FUN="/")
head(round(LR.RU.A, 2))

# check that this is identical to output from MakeAgePrior():
LR.RU.A[!is.finite(LR.RU.A)] <- 0
all(abs(LR.RU.A  - AP.griffin[["LR.RU.A.raw"]]) < 0.001)  # allow for rounding errors

## ----griffinSmoothFlatten, echo=FALSE, fig.cap="Smoothed & Flattened agepriors for overlapping generations (griffin example). Note the different y-axes; non-displayed age differences have an ageprior of zero for all R.", fig.height=8, fig.width=8, message=FALSE----
AP_list <- list()
AP_list[["Original"]] <- MakeAgePrior(Ped_griffin, Smooth=FALSE, Flatten=FALSE, Plot=FALSE)
AP_list[["Smooth"]] <- MakeAgePrior(Ped_griffin, Smooth=TRUE, Flatten=FALSE, Plot=FALSE)
AP_list[["Flatten"]] <- MakeAgePrior(Ped_griffin, Smooth=FALSE, Flatten=TRUE, Plot=FALSE)
AP_list[["Flatten + Smooth"]] <- MakeAgePrior(Ped_griffin, Smooth=TRUE, Flatten=TRUE, Plot=FALSE)

par(mfcol=c(2,2), mai=c(.7,.8,.4,.1))
for (x in names(AP_list)) {
  PlotAgePrior(AP_list[[x]], legend=FALSE)
  mtext(x, side=3, line=0.5, cex=1.5)
}

## ----tailstretch--------------------------------------------------------------
SmoothAP <- function(V, MinP) {
  Front <- max(1, min(which(V > MinP)), na.rm=TRUE)
  End <- min(max(which(V > MinP)), length(V), na.rm=TRUE)
  if (Front > 1 & V[Front] > 3*MinP)    V[Front -1] <- MinP
  if (End < length(V))   							  V[End +1] <- V[End]/2
  if ((End+1) < length(V) & V[End+1] > 3*MinP)  V[End+2] <- MinP
  # ... dip-fixing ...
  V
}

LR.RU.A <- apply(LR.RU.A, 2, SmoothAP, MinP = 0.001)

## ----PedComp-ex1, echo=TRUE, eval=TRUE----------------------------------------
AP.pedigree <- MakeAgePrior(Ped_griffin, Smooth=FALSE, Flatten=FALSE, 
                            quiet=TRUE, Plot=FALSE)
AP.default <- MakeAgePrior(MaxAgeParent = 3, quiet=TRUE, Plot=FALSE)

knitr::kable(list(AP.pedigree, AP.default), 
             caption = "Pedigree-based ageprior (left) and default, flat ageprior (right).",
             booktabs=TRUE)

## ----PedComp-ex1b, echo=TRUE, eval=TRUE---------------------------------------
AP.flattened <- MakeAgePrior(Ped_griffin, Smooth=FALSE, Flatten=TRUE, quiet=TRUE, Plot=FALSE)

## ----PedComp-ex1c, echo=FALSE, eval=TRUE--------------------------------------
if (is_html_output()) {
knitr::kable(AP.flattened, caption = "Flattened ageprior", format="html", 
             table.attr = "style='width:50%;'")
} else {
  knitr::kable(AP.flattened, caption = "Flattened ageprior", booktabs=TRUE)
}

## ----plotFlattenNW, echo=FALSE, fig.cap="Weights versus number of pairs with relationship R, for weighed average between pedigree-derived ageprior and flat 0/1 ageprior", fig.pos="hbt", fig.height=5, , out.width="60%"----
N <- 0:310
lambdaNW <- c(1e-4, -log(0.5) / c(50,100,200), 0.2)
W <- matrix(NA, 5, length(N))
for (i in 1:5) {
  W[i,] <- 1 - exp(-lambdaNW[i] * N)
}

par(mai=c(.8,.8,.1,.1))
plot(N, W[1,], type="n", ylim=c(0,1), las=1,
     xlab="Number of pairs (N_R)", ylab = "Weight (W_R)")
abline(h=seq(0,1,.2), col="lightgrey")
abline(h=0.5, col="darkgrey")
segments(x0=c(100,50,200), y0=-1, y1=0.5, col="darkgrey", 
         lty=1:3, lwd=c(1,1,1.5))
axis(side=2, at=0.5, labels=0.5, las=1)
for (i in 1:3) {
  lines(N, W[i+1,], lwd=3, col=c(4,1,4)[i], lty=c(2,1,3)[i])
}
for (i in c(1,5)) {
  lines(N, W[i,], lwd=2, col="darkgreen", lty=ifelse(i==1, 3, 2))
}

legend(210,0.45, legend=c(0.2, paste("-log(0.5)/", c(50,100,200)), 0.0001),
       col=c("darkgreen", 4,1,4, "darkgreen"), lty=c(2,2,1,3,3), lwd=c(2,3,3,3,2),
       title="lambdaNW", bg="white", inset=.03)

## ----griffinweights-----------------------------------------------------------
AP.griffin$Weights

## ----griffinweights2----------------------------------------------------------
NAK.R   # counts per relationship (column sums of tblA.R)

lambdaNW <- -log(0.5)/100     # default input value
W.R <-  1 - exp(-lambdaNW * NAK.R[1:5])
round(W.R, 4)

## ----timelineAU, echo=FALSE, fig.cap="Age difference between individual i and aunt/uncle j, and possible 'birth years' for i's parent k.", out.width="60%"----
knitr::include_graphics("agedif-AU.png")

## ----skip, fig.cap="", out.width = "70%"--------------------------------------
AP <- MakeAgePrior(MaxAgeParent = c(12, 8), Plot=FALSE)
AP[as.character(seq(1,12,by=2)), c("MS", "FS")] <- 0
PlotAgePrior(AP)

# check if valid ageprior:
chk <- sequoia:::CheckAP(AP) 

