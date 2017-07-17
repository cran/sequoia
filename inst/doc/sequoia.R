### R code from vignette source 'sequoia.Rnw'

###################################################
### code chunk number 1: sequoia.Rnw:31-65 (eval = FALSE)
###################################################
## install.packages("sequoia")  # only required first time
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
## # (maximum number of sibship-clustering iterations = 0)
## ParOUT <- sequoia(GenoM = Geno,
##                   LifeHistData = LH_HSg5,
##                   MaxSibIter = 0)
## names(ParOUT)
## # [1] "Specs"  "AgePriors"  "LifeHist"  "PedigreePar"  "MaybeParent"
## # "TotLikParents"
## #
## # run sequoia - sibship clustering & grandparent assignment
## # use parents assigned above (in 'ParOUT$PedigreePar')
## SeqOUT <- sequoia(GenoM = Geno,
##                   SeqList = ParOUT,
##                   MaxSibIter = 5)
## #
## # compare the assigned real and dummy parents to the true pedigree
## chk <- PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqOUT$Pedigree)
## chk$Counts
## #
## # save results
## save(SeqOUT, file="Sequoia_output_date.RData")
## write.table(SeqOUT$Pedigree, file="Pedigree.txt",
##             sep="\t", row.names=FALSE, quote=FALSE)
## 


###################################################
### code chunk number 2: sequoia.Rnw:86-88 (eval = FALSE)
###################################################
## GenoM <- as.matrix(read.table("MyGenoData.txt",
##                               row.names=1, header=FALSE))


###################################################
### code chunk number 3: sequoia.Rnw:103-104 (eval = FALSE)
###################################################
## system("cmd", input = "plink --file mydata --maf 0.4 --indep 50 5 2")


###################################################
### code chunk number 4: sequoia.Rnw:116-117 (eval = FALSE)
###################################################
## GenoM <- GenoConvert(InFile = "inputfile_for_sequoia.raw")


###################################################
### code chunk number 5: sequoia.Rnw:185-198 (eval = FALSE)
###################################################
## load("Sequoia_output_date.RData")  # if it was saved to disk
## ParOUT$Specs
## #   NumberIndivGenotyped NumberSnps GenotypingErrorRate MaxMismatch
## # 1                  920        200               1e-04           3
## #   Tfilter Tassign nAgeClasses MaxSibshipSize MaxSibIter
## # 1      -2     0.5           6            100          0
## #   DummyPrefixFemale DummyPrefixMale Complexity FindMaybeRel CalcLLR
## # 1                 F               M       full         TRUE    TRUE
## ParOUT$Specs$DummyPrefixFemale <- "D-FEM"
## ParOUT$Specs$DummyPrefixMale <- "D-MALE"
## SeqOUTX <- sequoia(GenoM = Geno,
##                   SeqList = list(Specs = ParOUT$Specs),
##                   MaxSibIter = 10)


###################################################
### code chunk number 6: sequoia.Rnw:202-203 (eval = FALSE)
###################################################
## SeqOUT <- sequoia(SeqList = ParOUT)


###################################################
### code chunk number 7: sequoia.Rnw:271-278 (eval = FALSE)
###################################################
## AP <- as.matrix(SeqOUT1$AgePriors)
## AP[AP>0] <- 0
## AP[1,c("MS", "PS")] <- 1
## AP[2,c("M", "P", "UA")] <- 1
## AP[3,c("MGM", "PGF", "MGF")] <- 1
## SeqOUT2 <- sequoia(SeqList=list(Specs=SeqOUT1$Specs, AgePriors=AP),
##                    MaxSibIter = 0)


###################################################
### code chunk number 8: sequoia.Rnw:417-422 (eval = FALSE)
###################################################
## TLL <- c(SeqOUT$TotLikParents, SeqOUT$TotLikSib)
## xv <- c(paste("p", 1:length(SeqOUT$TotLikParents)-1),
##         paste("s", 1:length(SeqOUT$TotLikSib)-1))
## plot(TLL, type="b", xaxt="n", xlab="Round")
## axis(1, at=1:length(TLL), labels=xv)


###################################################
### code chunk number 9: sequoia.Rnw:437-438 (eval = FALSE)
###################################################
## compareOUT <- PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqOUT$Pedigree)


###################################################
### code chunk number 10: sequoia.Rnw:445-453 (eval = FALSE)
###################################################
## compareOUT2 <- PedCompare(Ped1 = Ped_HSg5, Ped2 = ParOUT$Pedigree)
## compareOUT2$Counts["GG",,]
## #          dam sire
## # Total    130  170
## # Match    128  166
## # Mismatch   0    0
## # P1only     2    4
## # P2only     0    0


###################################################
### code chunk number 11: sequoia.Rnw:467-491 (eval = FALSE)
###################################################
## data(LH_HSg5, Ped_HSg5)
## GM <- SimGeno(Ped = Ped_HSg5, nSnp = 100, ErHQ = 1e-3)
## SeqX <- sequoia(GenoM = GM,  LifeHistData = LH_HSg5, MaxSibIter = 5)
## comp <- PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqX$Pedigree)
## comp$Counts
## # , , dam
## #    Total Match Mismatch P1only P2only
## # GG   561   550        6      5      0
## # GD   334   323        5      3      0
## # GT   890   873        9      8      0
## # DG    42    39        1      2      0
## # DD    28    28        0      0      0
## # DT    70    67        1      2      0
## # TT   964   940       10     10      4
## #
## # , , sire
## #    Total Match Mismatch P1only P2only
## # GG   502   498        3      1      0
## # GD   391   382        4      2      0
## # GT   890   880        7      3      0
## # DG    41    39        0      2      0
## # DD    29    28        0      1      0
## # DT    70    67        0      3      0
## # TT   965   947        7      6      5


###################################################
### code chunk number 12: sequoia.Rnw:495-503 (eval = FALSE)
###################################################
## with(comp$MergedPed,
##      comp$MergedPed[which(dam.1!=dam.2 & dam.1!=dam.r), ])
## #         id  dam.1 sire.1 dam.2 sire.2 id.r   dam.r sire.r
## # 205 b02032 a01166 b01133 F0045 b01165 <NA> nomatch   <NA>
## # 206 a02031 a01166 b01133 F0045 b01165 <NA> nomatch   <NA>
## # 720 b03011 a02106 b02181 F0047  M0020 <NA> nomatch b02181
## # 835 b04097 a03098 b03079 F0048  M0032 <NA> nomatch b03079
## # 897 a02030 a01166 b01133 F0045   <NA> <NA> nomatch   <NA>


###################################################
### code chunk number 13: sequoia.Rnw:507-513 (eval = FALSE)
###################################################
## comp$MergedPed[which(comp$MergedPed$dam.2=="F0045"), ]
## # -> only a02030, a02031 and b02032 are in this half-sibship
## # -> not wrongly merged
## comp$MergedPed[which(comp$MergedPed$dam.1=="a01166"), ]
## # -> true offspring of a01166 split over dummy mothers
## # F0027 and F0045


###################################################
### code chunk number 14: sequoia.Rnw:517-537 (eval = FALSE)
###################################################
## comp$MergedPed[which(comp$MergedPed$dam.2=="F0048"), ]
## # -> b04097 is sole individual in this sibship
## SeqX$MaybeParent[SeqX$MaybeParent$ID1=="b04097", ]
## #        ID1    ID2 Sex1 Sex2 AgeDif TopRel  LLR OH
## # 229 b04097 a03098    2    1      1     GP 0.29  1
## # 310 b04097 b04098    2    2      0     FA 1.17  1
## # -> dam in pedigree1 (a03098) more likely to be PO than U
## # (otherwise not in MaybeParent), but even more likely GP
## # Genotyping errors resulted in pair being OH at 1 SNP.
## 
## comp$MergedPed[which(comp$MergedPed$id=="F0048"), ]
## # -> parents of F0048 = grandparents of b04097: a02106 and b02158
## comp$MergedPed[which(comp$MergedPed$id=="a03098"), ]
## #         id  dam.1 sire.1  dam.2 sire.2 id.r dam.r sire.r
## # 288 a03098 a02106 b02158 a02106 b02158 <NA>  <NA>   <NA>
## # -> correct grandparents assigned, even if mother not assigned.
## 
## SeqX$MaybeRel[SeqX$MaybeRel$ID1=="b04097", ]
## #        ID1    ID2 Sex1 Sex2 AgeDif TopRel  LLR
## # 124 b04097 a03098    2    1      1     PO 1.31


###################################################
### code chunk number 15: sequoia.Rnw:544-548 (eval = FALSE)
###################################################
## BestConfig <- read.table("Colony/file/file.BestConfig",
##                          header=T, sep="", comment.char="")
## PedCompare(PedFile1 = "ExistingPedigree.txt",
##            Ped2 = BestConfig)


###################################################
### code chunk number 16: sequoia.Rnw:564-579 (eval = FALSE)
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
### code chunk number 17: sequoia.Rnw:583-598 (eval = FALSE)
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
##                       ifelse(dam.cat=="GG", ConfProb["GG"],
##                         ifelse(dam.cat %in% c("GD", "GR"), ConfProb["GD"],
##                           ifelse(dam.cat %in% c("DG", "RG"), ConfProb["DG"],
##                             ifelse(dam.cat %in% c("DD", "DR", "RD", "RR"),
##                                   ConfProb["DD"], NA)))))


###################################################
### code chunk number 18: sequoia.Rnw:607-636 (eval = FALSE)
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


