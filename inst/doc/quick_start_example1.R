## ----setup-q, include=FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval=TRUE, #cache=TRUE,
                      fig.height=4, fig.width=6, 
#                      fig.path='/figs',   DO NOT USE
                      fig.pos="ht!")
library(sequoia)

## -----------------------------------------------------------------------------
Geno <- Geno_HSg5

## ----example1-e---------------------------------------------------------------
ParOUT <- sequoia(GenoM = Geno,
                  LifeHistData = LH_HSg5,
                  Module = 'par')

## ----example1-f---------------------------------------------------------------
names(ParOUT)

## ----example1-g---------------------------------------------------------------
tail(ParOUT$PedigreePar)

## ----example1-h---------------------------------------------------------------
chk_par <- PedCompare(Ped1 = Ped_HSg5, Ped2 = ParOUT$PedigreePar)
chk_par$Counts["GG",,]   

## ----example1-i2, eval=TRUE---------------------------------------------------
ParOUT$AgePriors

## ----example1-j2, eval=TRUE---------------------------------------------------
data(SeqOUT_HSg5)

## ----example1-k---------------------------------------------------------------
summary_seq1 <- SummarySeq(SeqOUT_HSg5, Panels=c("sibships", "D.parents", "LLR"))
names(summary_seq1)

## ----example1-l---------------------------------------------------------------
chk <- PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqOUT_HSg5$Pedigree)
chk$Counts

