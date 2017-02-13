### R code from vignette source 'sequoia.Rnw'

###################################################
### code chunk number 1: sequoia.Rnw:38-67 (eval = FALSE)
###################################################
## install.packages("sequoia")  # only required first time
## setwd("E:/Sequoia/test")     # set the working directory
## library(sequoia)             # load the package
## #
## # get the example pedigree and life history data
## data(Ped_HSg5, LH_HSg5)
## tail(Ped_HSg5)
## #
## # simulate genotype data for 200 SNPs
## Geno <- SimGeno(Ped = Ped_HSg5, nSnp = 200)
## #
## # run sequoia - duplicate check & parentage assignment only
## ParOUT <- sequoia(GenoM = Geno,
##                   LifeHistData = LH_HSg5,
##                   MaxSibIter = 0)
## #
## # run sequoia - sibship clustering & grandparent assignment
## SeqOUT <- sequoia(GenoM = Geno,
##                   SeqList = ParOUT,
##                   MaxSibIter = 5)
## #
## # compare the assigned real and dummy parents to the true pedigree
## PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqOUT$Pedigree)
## #
## # save results
## save(SeqOUT, file="Sequoia_output_date.RData")
## write.table(SeqOUT$Pedigree, file="Pedigree_date.txt",
##             sep="\t", row.names=FALSE, quote=FALSE)
## 


###################################################
### code chunk number 2: sequoia.Rnw:104-106 (eval = FALSE)
###################################################
## system("cmd", input = "plink --file mydata --extract plink.prune.in
##   --recodeA --out inputfile_for_sequoia")


###################################################
### code chunk number 3: sequoia.Rnw:109-111 (eval = FALSE)
###################################################
## GenoM <- GenoConvert(InFile = "inputfile_for_sequoia.raw",
##                      OutFile = NA)


###################################################
### code chunk number 4: sequoia.Rnw:165-166 (eval = FALSE)
###################################################
## GenoM <- as.matrix(read.table("SimGeno.txt", row.names=1))


###################################################
### code chunk number 5: sequoia.Rnw:172-176 (eval = FALSE)
###################################################
## ParOUT <- sequoia(GenoM = GenoM,
##                   LifeHistData = read.table("LifeHistFile.txt",
##                                             header=TRUE),
##                   MaxSibIter=0)


###################################################
### code chunk number 6: sequoia.Rnw:204-206 (eval = FALSE)
###################################################
## load("Sequoia_output_date.RData")  # if it was saved to disk
## SeqOUT2 <- sequoia(SeqList = list(Specs = SeqOUT1$Specs))


###################################################
### code chunk number 7: sequoia.Rnw:212-213 (eval = FALSE)
###################################################
## SeqOUT2 <- sequoia(SeqList = SeqOUT1)


###################################################
### code chunk number 8: sequoia.Rnw:286-293 (eval = FALSE)
###################################################
## AP <- as.matrix(SeqOUT$AgePriors)
## AP[AP>0] <- 0
## AP[1,c("MS", "PS")] <- 1
## AP[2,c("M", "P", "UA")] <- 1
## AP[3,c("MGM", "PGF", "MGF")] <- 1
## SeqOUT2 <- sequoia(SeqList=list(Specs=SeqOUT$Specs, AgePriors=AP),
##                    MaxSibIter = 0)


###################################################
### code chunk number 9: sequoia.Rnw:346-351 (eval = FALSE)
###################################################
## TLL <- c(SeqOUT$TotLikParents, NA, SeqOUT$TotLikSib)
## xv <- c(paste0("p", 1:length(SeqOUT$TotLikParents)-1), " ",
##         paste0("s", 1:length(SeqOUT$TotLikSib)-1))
## plot(TLL, type="b", xaxt="n", xlab="Round")
## axis(1, at=1:length(TLL), labels=xv)


###################################################
### code chunk number 10: sequoia.Rnw:359-360 (eval = FALSE)
###################################################
## compareOUT <- PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqOUT$Pedigree)


###################################################
### code chunk number 11: sequoia.Rnw:393-410 (eval = FALSE)
###################################################
## data(SimGeno_example, LH_HSg5, Ped_HSg5)
## SeqOUTX <- sequoia(GenoM = SimGeno_example, LifeHistData = LH_HSg5)
## #
## # compare
## compare <- PedCompare(Ped1 = Ped_HSg5,
##                       Ped2 = SeqOUTX$Pedigree)
## compare$Counts
## # , , dam
## #    Total Match Mismatch P1only P2only
## # GG   130   127        2      1      0
## # GD    54    52        2      0      0
## # GT   182   179        2      1      0
## # DG     0     0        0      0      0
## # DD     0     0        0      0      0
## # DT     0     0        0      0      0
## # TT   960   179        2    779      0
## # ...


###################################################
### code chunk number 12: sequoia.Rnw:416-436 (eval = FALSE)
###################################################
## with(compare$MergedPed,
##      compare$MergedPed[which(dam.1!=dam.2 & dam.1!=dam.r), ])
## #     id  dam.1 sire.1 dam.2 sire.2   id.r   dam.r sire.r
## # b01137 a00001 b00009 F0007 b00009 b01137 nomatch   <NA>
## # a01139 a00001 b00009 F0007 b00009 a01139 nomatch   <NA>
## #
## compare$MergedPed[which(compare$MergedPed$dam.2=="F0007"), ]
## # --> only a01139 and b01137 are in this half-sibship
## compare$MergedPed[which(compare$MergedPed$dam.1=="a00001"), ]
## # --> a01139 and b01137 are the only ones not assigned to a00001
## #
## SeqOUTX$MaybeParent
## #       ID1    ID2 Sex1 Sex2 AgeDif Relx LLR_Rx_U TopRel LLR_R1_R2 OH
## # (...)
## # 10 a01063 b00013    1    2      1   PO    10.63     FS      1.37  1
## # 3  a01139 a00001    1    1      1   PO     7.32    2nd      4.79  2 <--
## # 8  b01168 b00002    2    2      1   PO     4.05    2nd      1.23  2
## # 2  b01137 a00001    2    1      1   PO     8.90    2nd      0.91  1 <--
## # 6  b01165 b00002    2    2      1   PO     4.14    2nd      0.45  2
## # 1  a01142 a00018    1    1      1   PO     9.89    2nd      0.17  2


###################################################
### code chunk number 13: sequoia.Rnw:442-445 (eval = FALSE)
###################################################
## SeqOUTX$Pedigree[SeqOUTX$Pedigree$id=="a01063", ]
## #        id   dam   sire LLRdam LLRsire LLRpair
## # 90 a01063 F0001 b00013   3.88    2.47    6.05


###################################################
### code chunk number 14: sequoia.Rnw:453-457 (eval = FALSE)
###################################################
## BestConfig <- read.table("Colony/file/file.BestConfig",
##                          header=T, sep="", comment.char="")
## PedCompare(PedFile1 = "ExistingPedigree.txt",
##            Ped2 = BestConfig)


###################################################
### code chunk number 15: sequoia.Rnw:473-488 (eval = FALSE)
###################################################
## data(LH_HSg5, Ped_HSg5)
## ARER <- array(dim=c(10,7,2))
## for (x in 1:10) {
##   cat(x, "\t", format(Sys.time(), "%H:%M:%S"), "\n")
##   GM <- SimGeno(Ped = Ped_HSg5, nSnp = 100)
##   SimOUT <- sequoia(GenoM = GM,  LifeHistData = LH_HSg5, quiet=TRUE)
##   CountX <- PedCompare(Ped1 = Ped_HSg5, Ped2 = SimOUT$Pedigree)$Counts
##   ARER[x,,1] <- rowMeans(CountX[, "Match", ] / CountX[, "Total", ])
##   ARER[x,,2] <- rowMeans((CountX[, "Mismatch", ] +
##                           CountX[, "P2only", ]) / CountX[, "Total", ])
## }
## dimnames(ARER) <- list(c(1:10), dimnames(CountX)[[1]], c("AR", "ER"))
## #
## # average assignment rate (AR) & error rate (ER) per category:
## round(apply(ARER, c(2:3), mean), 4)


###################################################
### code chunk number 16: sequoia.Rnw:492-512 (eval = FALSE)
###################################################
## ConfProb <- 1 - round(colMeans(ARER[, , "ER"]),3)
## #
## PedC <- PedCompare(Ped1 = Ped_HSg5,
##                    Ped2 = SeqOUT$Pedigree)$ConsensusPed
## PedC$dam.prob <- ConfProb[as.character(PedC$dam.cat)]
## PedC$sire.prob <- ConfProb[as.character(PedC$sire.cat)]
## #
## # or, when assuming that replacement of dummies by IDs of
## # non-genotyped individuals is free from error,
## PedC$dam.prob <- with(PedC,
##                       ifelse(dam.cat=="GG",
##                                ConfProb["GG"],
##                                ifelse(dam.cat %in% c("GD", "GR"),
##                                   ConfProb["GD"],
##                                   ifelse(dam.cat %in% c("DG", "RG"),
##                                      ConfProb["DG"],
##                                      ifelse(dam.cat %in% c("DD", "DR",
##                                                            "RD", "RR"),
##                                         ConfProb["DD"],
##                                         NA)))))


###################################################
### code chunk number 17: sequoia.Rnw:521-550 (eval = FALSE)
###################################################
## Rel.snp <- read.table("GT.grm.gz")
## Rel.id <- read.table("GT.grm.id", stringsAsFactors=FALSE)
## Rel.snp[,1] <- as.character(factor(Rel.snp[,1], labels=Rel.id[,2]))
## Rel.snp[,2] <- as.character(factor(Rel.snp[,2], labels=Rel.id[,2]))
## names(Rel.snp) <- c("IID2", "IID1", "SNPS", "R.SNP")
## Rel.snp <- Rel.snp[Rel.snp$IID1 != Rel.snp$IID2,]
## #
## library(pedantics)
## PedStats <- pedigreeStats(SeqOUT$Pedigree[,1:3], graphicalReport=FALSE,
##                           includeA=TRUE)
## Rel.ped <- as.data.frame.table(PedStats$Amatrix)
## names(Rel.ped) <- c("IID1", "IID2", "R.seq")
## #
## library(data.table)
## Rel.snp <- data.table(Rel.snp, key=c("IID1", "IID2"))
## Rel.ped <- data.table(Rel.ped, key=c("IID1", "IID2"))
## Rel.gt <- merge(Rel.snp[,c(1,2,4)], Rel.ped, all.x=TRUE)
## Rel.gt <- as.data.frame(Rel.gt)
## rm(PedStats, Rel.snp, Rel.ped)
## #
## round(cor(Rel.gt[, 3:4], use="pairwise.complete"),4)
## #
## library(hexbin)
## ColF <- function(n) rev(rainbow(n, start=0, end=4/6,
##                                 s=seq(.9,.6,length.out=n),v=.8))
## hexbinplot(Rel.gt$R.SNP~Rel.gt$R.ped, xbins=100, aspect=1, maxcnt=10^6.5,
##            trans=log10,inv=function(x) 10^x, colorcut=seq(0,1,length=14),
##            xlab="Pedigree relatedness", ylab="Genomic relatedness",
##            xlim=c(-.1,.9), ylim=c(-.1, .9), colramp=ColF, colorkey = TRUE)


