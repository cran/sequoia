## ----setup-q, include=FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval=TRUE, cache=TRUE,
                      fig.height=4, fig.width=6, 
#                      fig.path='/figs',   DO NOT USE
                      fig.pos="ht!")
library(sequoia)

## ----example1-a, eval=FALSE---------------------------------------------------
#  # Install the package. This is only the required the first time, or if you wish to update
#  install.packages("sequoia")
#  
#  # Load the package. This is required at the start of every new R session.
#  library(sequoia)
#  
#  # Load the example pedigree and life history data
#  data(Ped_HSg5, LH_HSg5)
#  
#  # Take a look at the data structure
#  tail(Ped_HSg5)
#  head(LH_HSg5)
#  
#  # or, in Rstudio, view the full dataframe:
#  View(Ped_HSg5)

## ----example1-d---------------------------------------------------------------
Geno <- SimGeno(Ped = Ped_HSg5, nSnp = 200)

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

## ----example1-i, eval=FALSE---------------------------------------------------
#  SeqOUT <- sequoia(GenoM = Geno, LifeHistData = LH_HSg5, Module = "ped")

## ----example1-j2, eval=TRUE---------------------------------------------------
ParOUT$AgePriors

## ----example1-j---------------------------------------------------------------
SeqOUT <- sequoia(GenoM = Geno,
                  SeqList = ParOUT,
                  Module = "ped")

## ----example1-k---------------------------------------------------------------
summary_seq1 <- SummarySeq(SeqOUT, Panels=c("sibships", "D.parents", "LLR"))
names(summary_seq1)

## ----example1-l---------------------------------------------------------------
chk <- PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqOUT$Pedigree)
chk$Counts

## ----example1-m, eval=FALSE---------------------------------------------------
#  save(SeqOUT, Geno, file="Sequoia_output_date.RData")
#  writeSeq(SeqList = SeqOUT, GenoM = Geno, folder = "Sequoia-OUT")

