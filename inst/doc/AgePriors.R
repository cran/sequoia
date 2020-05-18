## ----setup, include=FALSE-----------------------------------------------------
# output: 
#   bookdown::pdf_book: 
#   pdf_document:
library(knitr)
library(sequoia)
knitr::opts_chunk$set(echo = TRUE, 
                      fig.path='figs/', fig.height=4, fig.width=6,
                      fig.pos="h!")

# superseded: \newcommand{\grey}[1]{\textcolor{gray}{#1}}  (in main text, pdf only)

# adapted from: https://bookdown.org/yihui/rmarkdown-cookbook/font-color.html
MkGrey <- function(x) {    
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{gray}{%s}", x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: gray;'>%s</span>", x)
  } else x
}

## ----pipeline, echo=FALSE, fig.cap="Pipeline overview", out.width="90%"-------
if (knitr::is_latex_output()) {
  knitr::include_graphics("figs/agepriors_pipeline.pdf")
} else if (knitr::is_html_output()) {
  knitr::include_graphics("figs/agepriors_pipeline.png")
}

## ----Glossary, echo=FALSE, results="asis"-------------------------------------
Glossary <- matrix(c("$BY_i$", "Birth year", "Time unit (year, decade, month) of $i$'s birth/ hatching",
         "$A_{i,j}$", "Age difference",  "$BY_i - BY_j$, i.e. if $j$ is older than $i$, than $A_{i,j}$ is positive",
         "$R_{i,j}$", "Relationship", "Relationship between $i$ and $j$, e.g. parent-offspring",
         "$\\alpha_{A,R}$", "Ageprior", "Probability ratio of relationship $R$ given $A$, versus $R$ across all age differences"),
       ncol=3, byrow=TRUE)

knitr::kable(Glossary, col.names = c("Symbol", "Term", "Definition"),
            caption="Terms and abbreviations", escape=FALSE, booktabs=TRUE)

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
          "(M/P)(M/P/F)A" = "Mother's/Father's Maternal/Paternal/Full sibling")

Rels <- cbind(as.DF(Rels, c("Column", "Meaning")),
              "Version" = c(rep("all", 2), "from 1.0", "all", "all",
                            rep("up to 1.3 *", 3), "up to 1.3",
                            "from 2.0 *"))

knitr::kable(Rels, row.names=FALSE, booktabs=TRUE,
             caption = "AgePriors column names. *: From sequoia v2.0 only in AgePriorExtra")

## ----APDiscrete, fig.cap="Ageprior for discrete generations", out.width="70%"----
MakeAgePrior(Discrete = TRUE)

## ----APDiscrete2--------------------------------------------------------------
MakeAgePrior(Discrete = TRUE, MaxAgeParent = 2, Plot=FALSE)

## ----APMaxPO, fig.cap="Ageprior for overlapping generations, MaxAgeParent is 3 for dams and 2 for sires.", out.width="70%"----
MakeAgePrior(MaxAgeParent = c(3, 2))

## ----pedHead------------------------------------------------------------------
data(Ped_griffin)
tail(Ped_griffin)

## ----APall, fig.cap="Ageprior for overlapping generations (griffin example)", fig.pos="h!"----
AP.griffin <- MakeAgePrior(Ped_griffin, Smooth=FALSE, Flatten=FALSE, Return="all")

## ----APnames------------------------------------------------------------------
names(AP.griffin)

## ----tblAR--------------------------------------------------------------------
AP.griffin[["tblA.R"]]

## ----matAge-------------------------------------------------------------------
for (x in c("id", "dam", "sire")) {
  Ped_griffin[,x] <- as.character(Ped_griffin[,x])
}
Ped.g2 <- merge(Ped_griffin, 
                setNames(Ped_griffin[, c("id", "birthyear")], c("dam", "BY.dam")),
                all.x=TRUE)
Ped.g2$Age.dam <- with(Ped.g2, birthyear - BY.dam)
table(Ped.g2$Age.dam)

## ----sibAge-------------------------------------------------------------------
# this is part of the code in 'MakeAgePrior':
Ped.R <- setNames(Ped_griffin, c("id", "dam", "sire", "BirthYear"))  
RCM <- sapply(seq_along(Ped.R$id), GetRelCat, Ped.R, GenBack=1)  
# NOTE: GetRelCat is slow for very large pedigrees
AgeDifM <- outer(Ped.R$BirthYear, Ped.R$BirthYear, "-")
diag(AgeDifM) <- NA

MaxT <- 9   # determined internally based on LifeHistData
RR <- c("M", "P", "FS", "MS", "PS")
tblA.R <- matrix(NA, MaxT+1, length(RR)+1, dimnames=list(0:MaxT, c(RR, "X")))

tblA.R[, "FS"] <- table(factor(abs(AgeDifM[RCM == "FS"]), levels=0:MaxT)) /2
tblA.R[, "MS"] <- table(factor(abs(AgeDifM[RCM %in% c("FS", "MHS")]), levels=0:MaxT)) /2
tblA.R[, "PS"] <- table(factor(abs(AgeDifM[RCM %in% c("FS", "PHS")]), levels=0:MaxT)) /2

# Reference: age difference distribution across all pairs of individuals
tblA.R[, "X"] <- table(factor(abs(AgeDifM), levels=0:MaxT))/2
  
tblA.R

## ----calcPAR------------------------------------------------------------------
tblA.R <- AP.griffin[["tblA.R"]]    # with dam & sire columns
PA.R <- sweep(tblA.R, 2, STATS=colSums(tblA.R), FUN="/")
round(PA.R, 2)

## ----calcLRAR-----------------------------------------------------------------
LR.A.R <- sweep(PA.R[, 1:5], 1, STATS=PA.R[,"X"], FUN="/")
round(LR.A.R, 2)

## ----griffin_sim, eval=FALSE--------------------------------------------------
#  GenoS <- SimGeno(Ped_griffin, nSnp=400, ParMis=0.4)
#  
#  Ped_griffin$id <- as.character(Ped_griffin$id)
#  griffin.sex <- sapply(Ped_griffin$ID, function(x) substr(x, start=nchar(x), stop=nchar(x)))
#  LH.griffin <- data.frame(ID = Ped_griffin$ID,
#                           Sex = ifelse(griffin.sex=="F", 1, 2),
#                           BirthYear = Ped_griffin$birthyear)
#  
#  SeqOUT.griffin <- sequoia(GenoS, LH.griffin,
#                    MaxSibIter = 10,
#                    args.AP = list(Smooth = FALSE))

## ----griffinAPxtra, fig.cap="Ageprior with extra columns for Griffin example", fig.width=8----
data(SeqOUT_griffin)   # example output included with package
PlotAgePrior(SeqOUT_griffin$AgePriorExtra)

## ----OwnAP, eval=FALSE, error=FALSE-------------------------------------------
#  MyAgePrior <- MakeAgePrior(LifeHistData)   # with or without scaffold/old pedigree
#  # (...)  # edit by hand
#  SeqOUT <- sequoia(Genotypes, LifeHistData, SeqList(AgePriors = MyAgePrior))

## ----plotFlattenNW, echo=FALSE, fig.cap="Weights versus number of pairs with relationship R, for weighed average between pedigree-derived ageprior and flat 0/1 ageprior", fig.pos="hbt", fig.height=5, , out.width="60%"----
N <- 0:310
lambdaNW <- c(1e-4, -log(0.5) / c(50,100,200), 0.2)
W <- matrix(NA, 5, length(N))
for (i in 1:5) {
  W[i,] <- 1 - exp(-lambdaNW[i] * N)
}

par(mai=c(.8,.8,.1,.1))
plot(N, W[1,], type="n", ylim=c(0,1), las=1,
     xlab="Number of pairs N_R", ylab = "Weight")
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

## ----griffinSmoothFlatten, echo=FALSE, fig.cap="Smoothed & Flattened agepriors for overlapping generations (griffin example). Note the different y-axes; non-displayed age differences have an ageprior of zero for all R.", fig.height=8, fig.width=8, message=FALSE----
AP_list <- list()
AP_list[["Original"]] <- MakeAgePrior(Ped_griffin, Smooth=FALSE, Flatten=FALSE, Plot=FALSE)
AP_list[["Smooth"]] <- MakeAgePrior(Ped_griffin, Smooth=TRUE, Flatten=FALSE, Plot=FALSE)
AP_list[["Flat"]] <- MakeAgePrior(Ped_griffin, Smooth=FALSE, Flatten=TRUE, Plot=FALSE)
AP_list[["Flat + Smooth"]] <- MakeAgePrior(Ped_griffin, Smooth=TRUE, Flatten=TRUE, Plot=FALSE)

par(mfrow=c(2,2), mai=c(.7,.8,.4,.1))
for (x in names(AP_list)) {
  PlotAgePrior(AP_list[[x]], legend=FALSE)
  mtext(x, side=3, line=0.5, cex=1.5)
}

## ----griffinweights-----------------------------------------------------------
AP.griffin$Weights

## ----griffinNR----------------------------------------------------------------
N.R <- colSums(AP.griffin$tblA.R)
N.R

## ----griffinweights2----------------------------------------------------------
lambdaNW <- -log(0.5)/100   # the default
W.R <-  1 - exp(-lambdaNW * N.R)
round(W.R, 4)

## ----SimNoSire, message=FALSE, results="hide"---------------------------------
data(Ped_HSg5, LH_HSg5)
# Simulate data: 20% of mothers non-genotyped & 100% of fathers
GM <- SimGeno(Ped_HSg5, nSnp=400, ParMis=c(0.2, 1.0))  

ParOUT <- sequoia(GM, LH_HSg5, MaxSibIter = 0, Plot=FALSE)

## ----simNoSirePlot, fig.cap="Example: no fathers genotyped", out.width="70%"----
PlotAgePrior(ParOUT$AgePriors)

## ----NoSire2, eval=FALSE------------------------------------------------------
#  ParOUT <- sequoia(GM, LH_HSg5, MaxSibIter = 0,
#                    args.AP = list(MaxAgeParent=c(NA, 3),
#                                   Smooth = FALSE))

