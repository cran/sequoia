## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval=TRUE,
                      fig.height=4, fig.width=6, fig.pos="ht!") 
library(sequoia)

## ----RunPedComp---------------------------------------------------------------
library(sequoia)
data(SeqOUT_griffin, FieldMums_griffin, package="sequoia")

PCG  <- PedCompare(Ped1 = cbind(FieldMums_griffin,
                                sire = NA),
                   Ped2 = SeqOUT_griffin$Pedigree,
                   SNPd = SeqOUT_griffin$PedigreePar$id,
                   Symmetrical = TRUE, Plot=FALSE)

## ----MergedPed----------------------------------------------------------------
PCG$MergedPed[c(127:133), c("id", "dam.1", "dam.2", "dam.class")]

## ----id-idr-sep---------------------------------------------------------------
# subset some individuals:
these <- c("i177_2009_M", "i179_2009_M", "i165_2009_F", "i166_2009_F", "F0002",
           "F0007", "YellowPink", "PinkBlue")
knitr::kable(list(Ped1 = FieldMums_griffin[FieldMums_griffin$id %in% these, ], 
                  Ped2 = SeqOUT_griffin$Pedigree[SeqOUT_griffin$Pedigree$id %in% these, 1:3]),
             caption = "Subsets of Pedigree1 (left) and Pedigree2 (right)")

## ----id-idr-------------------------------------------------------------------
PCG$MergedPed[PCG$MergedPed$id %in% these, 
              c("id", "id.r", "dam.1", "dam.2", "dam.r")]

## ----DumMatch-----------------------------------------------------------------
head(PCG$DummyMatch[, -c(3:5)], n=6)

## ----Counts-plot--------------------------------------------------------------
PlotPedComp(PCG$Counts)

## ----Counts-------------------------------------------------------------------
PCG$Counts[,,"dam"]

## ----Mismatch1----------------------------------------------------------------
PCG$MergedPed[which(PCG$MergedPed$dam.class == "Mismatch"), c("id", "dam.1", "dam.2", "id.dam.cat")]

## ----MisMatch1----------------------------------------------------------------
PedM <- PCG$MergedPed[, c("id", "dam.1", "dam.2")]   # short-hand to minimise typing

# does the mismatch affect all of GreenBlue's offspring?
PedM[which(PedM$dam.1 == "GreenBlue"), ]
# > yes, these 4 are all of her known offspring

# does genetic mother i081_2005_F have any field-observed offspring?
PedM[which(PedM$dam.1 == "i081_2005_F"), ]
# no.

# does i081_2005_F have any other genetic offspring?
PedM[which(PedM$dam.2 == "i081_2005_F"), ]
# no.

## ----MisMatch2----------------------------------------------------------------
# why is this flagged as a mismatch?
PedM[which(PedM$dam.1 == "OrangeGreen"), ]
# all of OrangeGreen's offspring are in sibship F0001

PedM[which(PedM$dam.2 == "F0001"), c("id", "dam.1", "dam.2")]
# but sibship F0001 is split across two field mothers

## ----MisMatch3----------------------------------------------------------------
# as before
PedM[which(PedM$dam.1 == "YellowBlue"), ]
# something odd going on involving sibships F0003 & F0005

PedM[which(PedM$dam.2 %in% c("F0003", "F0005")), ]

## ----Ped1only-----------------------------------------------------------------
PCG$MergedPed[which(PCG$MergedPed$dam.class == "P1only"), 
              c("id", "id.r", "dam.1", "dam.2", "id.dam.cat")]

## -----------------------------------------------------------------------------
SeqOUT_griffin$DummyIDs[c(6,7), c("id", "dam", "BY.est", "NumOff", "O1", "O2", "O3", "O4")]

## -----------------------------------------------------------------------------
PedX <- SeqOUT_griffin$Pedigree
# INSTEAD, SOMETHING LIKE THIS:
Name_match <- matrix(c(PedX$dam[PedX$id=='i123_2007_F'], 'OrangeGreen',  
                       PedX$dam[PedX$id=='i165_2009_F'], 'PinkBlue',  
                       PedX$dam[PedX$id=='i121_2007_M'], 'GreenYellow'),
                     ncol = 2, byrow=TRUE)
Name_match

for (i in 1:nrow(Name_match)) {
  PedX$id[ PedX$id == Name_match[i,1] ] <- Name_match[i,2] 
  PedX$dam[ PedX$dam == Name_match[i,1] ] <- Name_match[i,2] 
}
  

## -----------------------------------------------------------------------------
head(SeqOUT_griffin$DummyIDs[, c('id','dam','sire','NumOff','O1','O2','O3')])

