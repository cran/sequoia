## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
library(sequoia)
knitr::opts_chunk$set(echo = TRUE, eval=TRUE,
                      auto_pdf = TRUE,
                      fig.height=4, fig.width=6, fig.pos="htb!")

## ----rel7, echo=FALSE, eval=TRUE, results="asis"------------------------------
knitr::kable(cbind(" "= c("PO", "FS", "HS", "GP", "FA", "HA", "U"),
                   Relationship = c("Parent - offspring",
                     "Full siblings",
                     "Half siblings",
                     "Grandparental",
                     "Full avuncular (aunt/uncle)",
                     "Half avuncular (+ other 3rd degree)",
                     "Unrelated")),
             caption = "Main relationships considered", escape=FALSE,
             booktabs=TRUE, position="!th")

## ----ErrSimSeq, echo=FALSE, eval=TRUE, fig.cap="(ref:caption-err)", out.width="80%"----
knitr::include_graphics("Errs-effect.png")

## ----ARER-cats, echo=FALSE, eval=TRUE, fig.cap="(ref:caption-ARERcats)", out.width="100%"----
knitr::include_graphics("ARER-cats.png")

## ----Errs-time, echo=FALSE, eval=TRUE, fig.cap="Runtimes are shorter with lower genotyping error rates and more SNPs", out.width="35%", fig.show='hold', fig.pos="!th"----
knitr::include_graphics(c("Errs-time-HSg5.png", "Errs-time-deer.png"))

## ----ErrM, echo=FALSE, eval=TRUE, results="asis"------------------------------
ErrM <- matrix(c("$(1-E/2)^2$", "$E\\times(1-E/2)$", "$(E/2)^2$",
                 "$E/2$", "$1-E$", "$E/2$",
                 "$(E/2)^2$", "$E\\times(1-E/2)$", "$(1-E/2)^2$"),
               3,3, byrow=TRUE,
               dimnames=rep(list(c("aa","Aa","AA")), 2))
knitr::kable(ErrM,
             caption = "Probability of observed genotype (columns) conditional on actual genotype (rows) and per-locus error rate E", escape=FALSE, booktabs=TRUE,
             position="!hb")

## ----AP-pipeline, echo=FALSE, eval=TRUE, fig.cap="Ageprior pipeline overview", out.width="100%"----
knitr::include_graphics("agepriors-pipeline.png")

## ----ClusterMatPat, echo=FALSE, eval=TRUE, fig.cap="Sibship clustering of individuals A-D", out.width="100%"----
knitr::include_graphics("clustering_mat_helps_pat.png")

## ----CoreFun, echo=FALSE, eval=TRUE, results="asis"---------------------------
tbl_core_fun <- rbind(
c("sequoia","Y","","","+","+","Pedigree, SeqList","SeqListSummary"),
c("GetMaybeRel","Y","","","+","+","Pairs","PlotRelPairs"),
c("CalcOHLLR","Y","Y","","+","+","","SeqListSummary"),
c("CalcPairLL","Y","+","Y","+","+","","PlotPairLL "),
c("MakeAgePrior","","+","","","","AgePrior","PlotAgePrior"),
c("SimGeno","","Y","","","","GenoM","SnpStats"),
c("PedCompare","","YY","","","","","PlotPedComp"),
c("ComparePairs","","Y+","+","","","Pairs","PlotRelPairs"))
colnames(tbl_core_fun) <- 
 c("Function","GenoM","Pedigree","Pairs","AgePrior","SeqList","Reusable output","Plot function")

knitr::kable(tbl_core_fun, align='lcccccll',
             caption = "Core functions with selected subset of their required (Y) and optional (+) input, and their output that can be re-used as input in other functions.", escape=FALSE, booktabs=TRUE)

## ----BY-example, echo=FALSE, eval=TRUE----------------------------------------
LH <- data.frame(ID = c("Alpha", "Beta", "Gamma"), 
                 Sex = c(2,2,1),
                 BirthYear = c(NA, 2012, NA),
                 BY.min = c(NA, NA, 2008),
                 BY.max = c(2010, NA, 2009))
knitr::kable(LH,
             caption = "Example LifeHistData with the optional columns for minimum + maximum birth year", booktabs=TRUE, position="b")

## ----argsAP, eval=TRUE--------------------------------------------------------
SeqOUT <- sequoia(GenoM = SimGeno_example, LifeHistData = LH_HSg5, Module='par',
                  args.AP = list(Discrete = TRUE), quiet=TRUE)
PlotAgePrior(SeqOUT$AgePriors)

## ----DoubleRels, echo=FALSE, eval=TRUE, results="asis"------------------------
DR <- rbind(PO = c("--", "--", "Y", "Y","","Y","","","", "Y"),
            FS = c("--", "--", "--", "--", "--", "Y","", "--", "Y", "Y"),
            HS = c("Y", "--", "(FS)", "Y", "Y ^2^", "Y", "Y ^2^","","", "Y"),
            GP = c("Y", "--", "Y", "^1^","","","","","", "Y"),
            FA = c("", "--","","","", "Y","","","", "Y"),
            HA = c("", "Y", "Y ^2^","","","","","","", "Y"),
            GGG = c("", "--","","","","", "^3^","","", "Y"),
            F1C = c("","","","","","","","","", "Y"))
colnames(DR) <- c(rownames(DR), "H1C", "U")
knitr::kable(DR,
             caption = "Double relationships between pairs of individuals. Abbreviations as before, and GGG=great-grandparent, F1C=full first cousins, H1C=half first cousins (parents are HS).",
             format = "pipe", booktabs=TRUE)  # forces ^1^ etc. through pandoc first 

## ----UseAge, echo=FALSE, eval=TRUE, results="asis"----------------------------
UA <- cbind("no" = c("Y", "Y", ""),
            "yes" = c("Y", "Y", "Y"),
            "extra" = c("Y", "", "Y"))
rownames(UA) <- c("$LR_{age} > 0$", "$LLR_{SNP} > T_{assign}$", 
                  "$LLR_{SNP} + LLR_{age} > T_{assign}$")

knitr::kable(UA, caption = "Parameter 'UseAge' (see text for definitions)", align='ccc', escape=FALSE, booktabs=TRUE, position="!ht")

## ----griffin-useAge-1, echo=TRUE, eval=TRUE-----------------------------------
data(Ped_griffin, SeqOUT_griffin, package="sequoia")
GenoX <- SimGeno(Ped_griffin, nSnp = 400, ParMis=0)
LLR_SNP <- CalcPairLL(Pairs = data.frame(ID1="i122_2007_M", ID2="i042_2003_F"), 
                      GenoM = GenoX, Plot=FALSE)
LLR_Age <- GetLLRAge(SeqOUT_griffin$AgePriorExtra, agedif=4, patmat=2)
knitr::kable(rbind(SNP = LLR_SNP[,colnames(LLR_Age)],
                   Age = LLR_Age,
                   "SNP + Age" = LLR_SNP[,colnames(LLR_Age)] +
                     LLR_Age),
             digits=1, booktabs=TRUE, 
             caption = "LLR example for a grandparent - grand-offspring pair")

## ----AP-2, eval=TRUE, out.width="70%", fig.align="center"---------------------
SeqOUT.B <- sequoia(GenoM = SimGeno_example, Err = 0.005,
                    LifeHistData = LH_HSg5, Module="par", Plot=FALSE, quiet=TRUE,
                    args.AP = list(MaxAgeParent = c(3,2), Smooth=FALSE))
PlotAgePrior(SeqOUT.B$AgePriors)

## ----SumSeq-parents, echo=FALSE, eval=TRUE, fig.cap="Number of parents assigned to genotyped individuals, split by category", out.width="80%"----
knitr::include_graphics("SummarySeq_deer.png")

## ----SumSeq-sibsize, echo=FALSE, eval=TRUE, fig.cap="Family sizes, split by genotyped (dark) or dummy (light) parent", out.width="80%"----
knitr::include_graphics("SummarySeq_deer_sibships.png")

## ----Relplot, echo=FALSE, eval=TRUE, fig.cap="Example plot from PlotRelPairs()", out.width="100%"----
# Rel.griffin <- GetRelM(PedPolish(Ped_griffin[Ped_griffin$birthyear <2004, ]), 
#                        patmat=TRUE, GenBack=1)
# PlotRelPairs(Rel.griffin)   # poor figure quality
knitr::include_graphics("RelPlot-griffin.png")

## ----PedComp-ex1, echo=FALSE, eval=TRUE---------------------------------------
Ped1 <- data.frame(id = c("Alpha", "Beta", "Abigail", "Aster", "Blossom", "Bob", "-"),
                    dam = c("","", "Alpha", "Alpha", "Beta", "Beta", ""),
                    sire = c("","", "Zorro", "Yann", "", "Zorro", ""))
Ped2 <- data.frame(id = c("Alpha", "-", "Abigail", "Aster", "Blossom", "Bob", "F0007"),
                   dam = c("", "", "Alpha", "Alpha", "F0007", "F0007", "Alpha"),
                   sire = c("", "", "Zorro", "M0004", "", "Zorro", ""))
knitr::kable(list(Ped1, Ped2), 
             caption = "Example comparison between Pedigree1 (left) and Pedigree2 (right)", booktabs=TRUE)

## ----RunPedComp, eval=TRUE----------------------------------------------------
PCG  <- PedCompare(Ped1 = cbind(FieldMums_griffin,
                                sire = NA),
                   Ped2 = SeqOUT_griffin$Pedigree)

# pedigrees side-by-side (subset of columns because no field-observed sires here)
PCG$MergedPed[127:133, c("id", "dam.1", "dam.2", "dam.r", "id.dam.cat", "dam.class")]
#  Non-genotyped field mums have a two-colour code (e.g. 'BlueRed')

PCG$MergedPed[c(137,138, 6, 128,129,7), c("id", "id.r", "dam.1", "dam.2", "dam.r")]
# column 'id': ids common to both pedigrees, plus those only in Pedigree2 
# column 'id.r': 'consensus' ids, plus those only occurring in Pedigree1  

# dummy individuals from Ped2 with their best-matching non-genotyped individual in Ped1
head(PCG$DummyMatch[, -c(3:5)])
# 'nomatch' in the `id.1` column: none of the siblings in Ped2 had a field-observed
# mother, or there is a mismatch.

# Total number of matches & mismatches:
PCG$Counts
# dim1: category: genotyped/dummy/total
# dim2: classification: total (could-have-had-parent)/ match/ mismatch/ 
#       P1only (no parent in Ped2)/ P2only
# dim3: dam/sire

## ----checks, eval=TRUE--------------------------------------------------------
# ~~ Mismatches ~~
PCG$MergedPed[which(PCG$MergedPed$dam.class == "Mismatch"), 
               c("id", "dam.1", "dam.2", "id.dam.cat")]

PedM <- PCG$MergedPed[, c("id", "dam.1", "dam.2")]   # short-hand to minimise typing

# who are the mismatching individual's siblings according to Ped1 & Ped2? 
PedM[which(PedM$dam.1 == "GreenBlue"), ]
PedM[which(PedM$dam.2 == "i081_2005_F"), ]

# who are their maternal grandparents?
PedM[which(PedM$id.1 == "GreenBlue"), ]
PedM[which(PedM$id.2 == "i081_2005_F"), ]

## ----calcRped, fig.cap="", fig.width=4, out.width="60%"-----------------------
Rped.old <- CalcRped(Ped_griffin, OUT="DF")
Rped.new <- CalcRped(SeqOUT_griffin$Pedigree, OUT="DF")
library(data.table)
Rped.both <- merge(data.table(Rped.old, key=c("IID1", "IID2")),
                   data.table(Rped.new, key=c("IID1", "IID2")), 
                   all=TRUE, suffixes=c(".old", ".new"))

cor(Rped.both$R.ped.new, Rped.both$R.ped.old, use="pairwise.complete")
plot(Rped.both$R.ped.old, Rped.both$R.ped.new, pch=16, cex=0.3,
     xlab = 'R (old ped)', ylab = 'R (new ped)')

## ----eval=TRUE, fig.cap="Relatedness hexbinplot example", out.width="80%", position="!ht"----
knitr::include_graphics("hexbin.png")

## ----parent-selfed, echo=FALSE, eval=TRUE, fig.cap="Two configurations with identical likelihoods", out.width="20%", fig.pos="!ht", fig.align="center"----
knitr::include_graphics("selfed_parent.png")

