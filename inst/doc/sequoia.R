### R code from vignette source 'sequoia.Rnw'

###################################################
### code chunk number 1: sequoia.Rnw:29-62 (eval = FALSE)
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
## writeSeq(SeqList = SeqOUT, GenoM = Geno, PedComp = chk,
##           folder = "Sequoia-OUT")


###################################################
### code chunk number 2: sequoia.Rnw:69-124 (eval = FALSE)
###################################################
## install.packages("sequoia")  # only required first time
## library(sequoia)             # load the package
## #
## # save the genotype data with 1 row per individual, 1 column per SNP
## # (0/1/2/NA), e.g. in PLINK:
## # plink --file mydata --geno 0.1 --maf 0.3 --indep 50 5 2
## # plink --file mydata --extract plink.prune.in --out
## #
## # read in genotype data
## # if already coded as 0/1/2, with missing=-9:
## Geno <- as.matrix(read.csv("mydata.csv", header=FALSE, row.names=1))
## # for many other  input formats:
## Geno <- GenoConvert(InFile = "mydata.ped", InFormat="ped")
## #
## # read in lifehistory data: order ID-Sex-birthyear, column names ignored
## LH <- read.table("LifeHistoryData.txt", header=T)
## #
## # duplicate check & parentage assignment (takes few minutes)
## # (maximum number of sibship-clustering iterations = 0)
## # if genotyping error rate is unknown, start of high
## ParOUT <- sequoia(GenoM = Geno,  LifeHistData = LH_HSg5,
##                   MaxSibIter = 0, Err=0.01, MaxMismatch=10)
## #
## # inspect duplicates (intentional or accidental)
## ParOUT$DupGenotype
## # (...)
## #
## # check if distr. of age-differences for each relative type is sensible
## PlotAgePrior(ParOUT$AgePriors)
## #
## # compare assigned parents to field pedigree (check column order!)
## FieldPed <- read.table("FieldPed.txt", header=T)
## PC.par <- PedCompare(Ped1 = FieldPed[, c("id", "dam", "sire")],
##            Ped2 = ParOUT$PedigreePar)
## PC.par$Counts["TT",,]
## # (...)
## #
## # calculate Mendelian errors per SNP (works also w field pedigree)
## stats <- SnpStats(Geno, ParOUT$PedigreePar)
## #
## # polish dataset: remove one indiv. from each duplicate pair
## # (1st one, or one w lowest call rate) & drop high error rate SNPs
## Geno2 <- Geno[!rownames(Geno) %in% ParOUT$DupGenotype$ID2, ]
## Geno2 <- Geno2[, -which(stats[,"ER"]>50)]
## # (check histogram for sensible threshold)
## #
## # iterate the above as necessary
## #
## # run full pedigree reconstruction (may take few hours)
## SeqOUT <- sequoia(GenoM = Geno2,
##                   MaxSibIter = 20,
##                   Err=0.001)
## # inspect no. assigned parents, proportion dummy parents, etc.
## SummarySeq(SeqOUT)
## # (see Example 1 for saving results)


###################################################
### code chunk number 3: sequoia.Rnw:157-162 (eval = FALSE)
###################################################
## GenoM <- as.matrix(read.table("MyGenoData.txt",
##                               row.names=1, header=FALSE))
## # or
## GenoM <- as.matrix(read.csv("MyGenoData.csv",
##                               row.names=1, header=FALSE))


###################################################
### code chunk number 4: sequoia.Rnw:178-179 (eval = FALSE)
###################################################
## system("cmd", input = "plink --file mydata --maf 0.3 --indep 50 5 2")


###################################################
### code chunk number 5: sequoia.Rnw:193-194 (eval = FALSE)
###################################################
## GenoM <- GenoConvert(InFile = "inputfile_for_sequoia.raw")


###################################################
### code chunk number 6: sequoia.Rnw:287-300 (eval = FALSE)
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
### code chunk number 7: sequoia.Rnw:304-306 (eval = FALSE)
###################################################
## SeqOUT <- sequoia(GenoM = Geno,
##                   SeqList = ParOUT)


###################################################
### code chunk number 8: sequoia.Rnw:399-409 (eval = FALSE)
###################################################
## AP <- matrix(0, nrow=5, ncol=9,
##              dimnames=list(0:4,
##                            c("M","P","MGM","PGF","MGF","FS","MS","PS","UA")))
## AP[1,c("MS", "PS")] <- 1
## AP[2,c("M", "P", "UA")] <- 1
## AP[3,c("MGM", "PGF", "MGF")] <- 1
## SeqOUT <- sequoia(GenoM = Geno,
##                   LifeHistData = LH,
##                   SeqList=list(AgePriors=AP),
##                    MaxSibIter = 0)


###################################################
### code chunk number 9: sequoia.Rnw:413-419 (eval = FALSE)
###################################################
## AP <- MakeAgePrior(Ped=OldPed, LifeHistData = LH,
##                    Flatten = FALSE, Smooth = FALSE)
## SeqOUT <- sequoia(GenoM = Geno,
##                   LifeHistData = LH,
##                   SeqList=list(AgePriors=AP),
##                    MaxSibIter = 0)


###################################################
### code chunk number 10: sequoia.Rnw:423-427 (eval = FALSE)
###################################################
## SeqOUT <- sequoia(GenoM = Geno,
##                   LifeHistData = LH,
##                   args.AP = list(Flatten = FALSE, Smooth = FALSE),
##                   MaxSibIter = 10)


###################################################
### code chunk number 11: sequoia.Rnw:431-442 (eval = FALSE)
###################################################
## AP <- MakeAgePrior(Ped=OldPed, LifeHistData = LH,
##                    Flatten = FALSE, Smooth = FALSE)
## SeqOUT1 <- sequoia(GenoM = Geno,
##                   LifeHistData = LH,
##                   SeqList=list(AgePriors=AP),
##                    MaxSibIter = 0)
## SeqOUT2 <- sequoia(GenoM = Geno,
##                    LifeHistData = LH,
##                    SeqList = list(AgePriors = AP,
##                                   PedigreePar = SeqOUT1$PedigreePar),
##                    MaxSibIter = 10)


###################################################
### code chunk number 12: sequoia.Rnw:558-563 (eval = FALSE)
###################################################
## TLL <- c(SeqOUT$TotLikParents, SeqOUT$TotLikSib)
## xv <- c(paste("p", 1:length(SeqOUT$TotLikParents)-1),
##         paste("s", 1:length(SeqOUT$TotLikSib)-1))
## plot(TLL, type="b", xaxt="n", xlab="Round")
## axis(1, at=1:length(TLL), labels=xv)


###################################################
### code chunk number 13: sequoia.Rnw:576-577 (eval = FALSE)
###################################################
## save(SeqList, LHdata, Geno, file="Sequoia_output_date.RData")


###################################################
### code chunk number 14: sequoia.Rnw:581-583 (eval = FALSE)
###################################################
## load("Sequoia_output_date.RData")
## # 'SeqList' and 'LHdata' will appear in R environment


###################################################
### code chunk number 15: sequoia.Rnw:588-589 (eval = FALSE)
###################################################
## writeSeq(SeqList, GenoM = Geno, folder=paste("Sequoia_OUT", Sys.Date()))


###################################################
### code chunk number 16: sequoia.Rnw:594-595 (eval = FALSE)
###################################################
## writeSeq(SeqList, OutFormat="xls", file="Sequoia_OUT.xlsx")


###################################################
### code chunk number 17: sequoia.Rnw:598-601 (eval = FALSE)
###################################################
## library(xlsx)
## write.xlsx(Geno, file = "Sequoia_OUT.xlsx", sheetName="Genotypes",
##       col.names=FALSE, row.names=TRUE, append=TRUE, showNA=FALSE)


###################################################
### code chunk number 18: sequoia.Rnw:610-611 (eval = FALSE)
###################################################
## compareOUT <- PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqOUT$Pedigree)


###################################################
### code chunk number 19: sequoia.Rnw:618-626 (eval = FALSE)
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
### code chunk number 20: sequoia.Rnw:640-675 (eval = FALSE)
###################################################
## data(LH_HSg5, Ped_HSg5)
## 
## GM <- SimGeno(Ped = Ped_HSg5, nSnp = 200, ErHQ = 1e-3)
## #
## LH <- LH_HSg5
## LH$BY[sample.int(nrow(LH), round(nrow(LH)*0.2))] <- NA
## LH$Sex[sample.int(nrow(LH), round(nrow(LH)*0.2))] <- NA
## #
## # run sequoia, with max 5 iterations of full pedigree reconstruction
## SeqX <- sequoia(GenoM = GM,  LifeHistData = LH, MaxSibIter = 5)
## #
## #check the number of mismatches in the full pedigree
## comp <- PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqX$Pedigree)
## comp$Counts
## # , , dam
## #
## #    Total Match Mismatch P1only P2only
## # GG   529   522        4      3      0
## # GD   367   363        4      0      0
## # GT   892   885        4      3      0
## # DG    39    39        0      0      0
## # DD    29    29        0      0      0
## # DT    68    68        0      0      0
## # TT   961   953        4      3      1
## #
## # , , sire
## #
## #    Total Match Mismatch P1only P2only
## # GG   550   549        1      0      0
## # GD   343   337        5      1      0
## # GT   892   886        5      1      0
## # DG    38    38        0      0      0
## # DD    30    30        0      0      0
## # DT    68    68        0      0      0
## # TT   960   954        5      1      0


###################################################
### code chunk number 21: sequoia.Rnw:680-687 (eval = FALSE)
###################################################
## # error rate:
## (4+1+5+0)/(2*960)
## #[1] 0.005208333
## 
## # correct assignment rate
## (953+954)/(2*960)
## #[1] 0.9932292


###################################################
### code chunk number 22: sequoia.Rnw:691-702 (eval = FALSE)
###################################################
## comp$Mismatch
##  #     id  dam.1 sire.1 dam.2 sire.2   id.r  dam.r  sire.r Cat Parent
##  # b05019 a04004 b04002 F0003  M0004 b05019 a04001  b04002  GG    dam
##  # b05018 a04004 b04002 F0003  M0004 b05018 a04001  b04002  GG    dam
##  # a05017 a04004 b04002 F0003  M0004 a05017 a04001  b04002  GG    dam
##  # b05020 a04004 b04002 F0003  M0004 b05020 a04001  b04002  GG    dam
##  # b05164 a04053 b04048 F0047  M0031 b05164 a04053 nomatch  GG   sire
##  # a05090 a04053 b04164 F0047  M0031 a05090 a04053 nomatch  GD   sire
##  # b05092 a04053 b04164 F0047  M0031 b05092 a04053 nomatch  GD   sire
##  # a05091 a04053 b04164 F0047  M0031 a05091 a04053 nomatch  GD   sire
##  # a04004 a03173 b03044 F0031  M0009 a04004 a03173  b03093  GD   sire


###################################################
### code chunk number 23: sequoia.Rnw:709-714 (eval = FALSE)
###################################################
## SeqX$DummyIDs[SeqX$DummyIDs$id=="F0003", ]
## #      id   dam  sire LLRdam LLRsire LLRpair sex BY.est BY.min BY.max NumOff     O1
## # 3 F0003 F0031 M0007   9.34   10.79    4.11   1      5      5      5     12 b05019
## #       O2     O3     O4     O5     O6     O7     O8     O9    O10    O11    O12
## # 3 a05017 b05173 a05174 b05175 a05176 b05037 b05038 b05040 a05039 b05020 b05018


###################################################
### code chunk number 24: sequoia.Rnw:717-722 (eval = FALSE)
###################################################
## Ped_HSg5[Ped_HSg5$id %in% c("a04001", "a04004", "b04002"), ]
## #         id    dam   sire
## # 617 a04001 a03173 b03044
## # 618 b04002 a03173 b03044
## # 620 a04004 a03173 b03044


###################################################
### code chunk number 25: sequoia.Rnw:728-736 (eval = FALSE)
###################################################
## PedM <- comp$MergedPed  # just to save typing
## #
## PedM[which(PedM$sire.2=="M0031"), ]
## #         id  dam.1 sire.1 dam.2 sire.2 id.r  dam.r  sire.r
## # 877 a05090 a04053 b04164 F0047  M0031 <NA> a04053 nomatch
## # 878 b05164 a04053 b04048 F0047  M0031 <NA> a04053 nomatch
## # 879 b05092 a04053 b04164 F0047  M0031 <NA> a04053 nomatch
## # 880 a05091 a04053 b04164 F0047  M0031 <NA> a04053 nomatch


###################################################
### code chunk number 26: sequoia.Rnw:742-756 (eval = FALSE)
###################################################
## PedM[which(PedM$sire.1=="b04164"), ]
## #         id  dam.1 sire.1  dam.2 sire.2 id.r  dam.r  sire.r
## # 846 b05175 a04001 b04164  F0003  M0028 <NA> a04001  b04164
## # 847 a05166 a04122 b04164 a04122  M0028 <NA>   <NA>  b04164
## # 848 a05176 a04001 b04164  F0003  M0028 <NA> a04001  b04164
## # 849 a05089 a04053 b04164  F0047  M0028 <NA> a04053  b04164
## # 850 a05167 a04122 b04164 a04122  M0028 <NA>   <NA>  b04164
## # 851 a05174 a04001 b04164  F0003  M0028 <NA> a04001  b04164
## # 852 b05173 a04001 b04164  F0003  M0028 <NA> a04001  b04164
## # 853 b05165 a04122 b04164 a04122  M0028 <NA>   <NA>  b04164
## # 877 a05090 a04053 b04164  F0047  M0031 <NA> a04053 nomatch
## # 879 b05092 a04053 b04164  F0047  M0031 <NA> a04053 nomatch
## # 880 a05091 a04053 b04164  F0047  M0031 <NA> a04053 nomatch
## # 897 a05168 a04122 b04164   <NA>   <NA> <NA>   <NA>    <NA>


###################################################
### code chunk number 27: sequoia.Rnw:765-775 (eval = FALSE)
###################################################
## PedM[which(PedM$dam.1=="a03173"), ]
## #         id  dam.1 sire.1 dam.2 sire.2   id.r  dam.r sire.r
## # 619 a04003 a03173 b03044 F0031  M0007   <NA> a03173 b03044
## # 636 b04080 a03173 b03093 F0031  M0009   <NA> a03173 b03093
## # 639 b04079 a03173 b03093 F0031  M0009   <NA> a03173 b03093
## # 640 a04004 a03173 b03044 F0031  M0009   <NA> a03173 b03093 <--
## # 643 a04078 a03173 b03093 F0031  M0009   <NA> a03173 b03093
## # 645 a04077 a03173 b03093 F0031  M0009   <NA> a03173 b03093
## # 968  M0004 a03173 b03044 F0031  M0007 b04002 a03173 b03044
## # 970  F0003 a03173 b03044 F0031  M0007 a04001 a03173 b03044


###################################################
### code chunk number 28: sequoia.Rnw:782-788 (eval = FALSE)
###################################################
## DyadCompare(Ped_HSg5, SeqX$PedigreePar)
## #     RC.2
## # RC.1     FS     HS      U
## #   FS    561    325    500
## #   HS      0   2131   2579
## #   U       0      0 416644


###################################################
### code chunk number 29: sequoia.Rnw:796-800 (eval = FALSE)
###################################################
## BestConfig <- read.table("Colony/file/file.BestConfig",
##                          header=T, sep="", comment.char="")
## PedCompare(Ped1 = ExistingPedigree,
##            Ped2 = BestConfig)


###################################################
### code chunk number 30: sequoia.Rnw:816-833 (eval = FALSE)
###################################################
## data(SimGeno_example, LH_HSg5, package="sequoia")
## SeqOUT <- sequoia(GenoM = SimGeno_example[, 1:100],
##                    LifeHistData = LH_HSg5, MaxSibIter = 5)
## ConfPr <- EstConf(Ped = SeqOUT$PedigreePar,
##                   LifeHistData = LH_HSg5,
##                   Specs = SeqOUT$Specs, Full = TRUE,
##                   nSim = 3, ParMis = 0.4)
## # , , mean
## #      GG    GD    GT  DG  DD  DT    TT
## # dam   1 0.950 0.980 NaN NaN NaN 0.980
## # sire  1 0.986 0.995 NaN NaN NaN 0.995
## #
## # , , min
## #      GG    GD    GT  DG  DD  DT    TT
## # dam   1 0.896 0.957 NaN NaN NaN 0.957
## # sire  1 0.957 0.986 NaN NaN NaN 0.986
## 


###################################################
### code chunk number 31: sequoia.Rnw:839-852 (eval = FALSE)
###################################################
## PedC <- PedCompare(Ped1 = Ped_HSg5,
##                    Ped2 = SeqOUT$Pedigree)$ConsensusPed
## ConfProb <- cbind(ConfPr[,,"mean"],
##                   "U" = NA,  # Ungenotyped, parent taken from Ped1
##                   "X" = NA)  # no parent in either pedigree
## 
## PedC$dam.cat2 <- PedC$dam.cat
## PedC$dam.cat2[PedC$dam.cat == "GR"] <- "GD"
## PedC$dam.cat2[PedC$dam.cat == "RG"] <- "DG"
## PedC$dam.cat2[PedC$dam.cat %in% c("DD", "DR", "RD", "RR")] <- "DD"
## PedC$dam.prob <- ConfProb["dam", as.character(PedC$dam.cat2)]
## 
## # and analogously for sires.


###################################################
### code chunk number 32: sequoia.Rnw:861-890 (eval = FALSE)
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


###################################################
### code chunk number 33: sequoia.Rnw:930-946 (eval = FALSE)
###################################################
## cand.dams <- read.table("Candidate_dams.txt", header=TRUE,
##                         stringsAsFactors=FALSE)
## # cdam has columns 'id' and 'dam', and does not need entries for all ids
## cand.par <- cbind(cand.dams, sire=NA)
## 
## par.herm <- sequoia(GenoM = Geno,
##                     LifeHistData = LH,
##                     SeqList = list(PedigreePar = cand.par),
##                     MaxSibIter=0)
## # In combination with maxSibIter=0, PedigreePar is a pedigree prior
## 
## # re-use all settings (including the newly assigned parents):
## seq.herm <- sequoia(GenoM = Geno,
##                     LifeHistData = LH,
##                     SeqList = par.herm,
##                     MaxSibIter = 10)


