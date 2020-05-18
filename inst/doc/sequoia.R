### R code from vignette source 'sequoia.Rnw'

###################################################
### code chunk number 1: sequoia.Rnw:28-30
###################################################
options(prompt=" ",continue=" ")
library(sequoia)


###################################################
### code chunk number 2: sequoia.Rnw:40-72 (eval = FALSE)
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
##                   MaxSibIter = 10)
## #
## # compare the assigned real and dummy parents to the true pedigree
## chk <- PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqOUT$Pedigree)
## chk$Counts["TT",,]
## #
## # save results
## save(SeqOUT, file="Sequoia_output_date.RData")
## writeSeq(SeqList = SeqOUT, GenoM = Geno, folder = "Sequoia-OUT")


###################################################
### code chunk number 3: sequoia.Rnw:79-129 (eval = FALSE)
###################################################
## # read in genotype data if already coded as 0/1/2, with missing=-9:
## Geno <- as.matrix(read.csv("mydata.csv", header=FALSE, row.names=1))
## CheckGeno(Geno)
## # read in many other input formats (not .vcf (yet)):
## Geno <- GenoConvert(InFile = "mydata.ped", InFormat="ped")
## #
## # read in lifehistory data: ID-Sex-birthyear, column names ignored
## # optional: minimum & maximum birth year, when not exactly known
## LH <- read.table("LifeHistoryData.txt", header=T)
## #
## # duplicate check & parentage assignment (takes few minutes)
## # (maximum number of sibship-clustering iterations = 0)
## ParOUT <- sequoia(GenoM = Geno,  LifeHistData = LH_HSg5,
##                   MaxSibIter = 0, Err=0.005,
##                   quiet = FALSE, Plot = TRUE)
## #
## # inspect duplicates (intentional or accidental)
## ParOUT$DupGenotype
## #
## # compare assigned parents to field pedigree (check column order!)
## FieldPed <- read.table("FieldPed.txt", header=T)
## PC.par <- PedCompare(Ped1 = FieldPed[, c("id", "dam", "sire")],
##                      Ped2 = ParOUT$PedigreePar)
## PC.par$Counts["TT",,]
## #
## # calculate Mendelian errors per SNP (works also w field pedigree)
## stats <- SnpStats(Geno, ParOUT$PedigreePar)
## MAF <- ifelse(stats[,"AF"] <= 0.5, stats[,"AF"], 1-stats[,"AF"])
## #
## # ..........................................................
## # polish dataset: remove one indiv. from each duplicate pair
## # & drop low call rate samples
## # & drop SNPs with high error rate and/or low MAF
## Geno2 <- Geno[!rownames(Geno) %in% ParOUT$DupGenotype$ID2, ]
## Geno2 <- Geno2[, -which(stats[,"Err.hat"]>0.05 | MAF < 0.1)]
## #
## Indiv.Mis <- apply(Geno2, 1, function(x) sum(x == -9)) / ncol(Geno2)
## Geno2 <- Geno2[Indiv.Mis < 0.2, ]
## # check histograms for sensible thresholds, iterate if necessary
## #
## # run full pedigree reconstruction (may take up to a few hours)
## # including re-run of parentage assignment
## SeqOUT <- sequoia(GenoM = Geno2,
##                   LifeHistData = LH_HSg5,
##                   MaxSibIter = 20,
##                   Err = 0.001)
## #
## # inspect assigned parents, proportion dummy parents, etc.
## SummarySeq(SeqOUT)
## # (see Example 1 for saving results)


###################################################
### code chunk number 4: sequoia.Rnw:387-389 (eval = FALSE)
###################################################
## SeqOUT <- sequoia(GenoM = Geno, LifeHistData = LH,
##                   args.AP = list(Discrete = TRUE))


###################################################
### code chunk number 5: sequoia.Rnw:394-400 (eval = FALSE)
###################################################
## APfromOld <- MakeAgePrior(Pedigree = MyOldPedigree,
##                           LifeHistData = LH,
##                           Smooth = TRUE)
## SeqOUT <- sequoia(GenoM = Geno,
##                   SeqList = list(AgePriors = APfromOld),
##                   MaxSibIter = 10)


###################################################
### code chunk number 6: sequoia.Rnw:406-407 (eval = FALSE)
###################################################
## APdraft <- MakeAgePrior(MaxAgeParent = c(4, 5))  # dams, sires


###################################################
### code chunk number 7: sequoia.Rnw:458-460 (eval = FALSE)
###################################################
## GenoM <- as.matrix(read.table("MyGenoData.txt",
##                               row.names=1, header=FALSE))


###################################################
### code chunk number 8: sequoia.Rnw:476-477 (eval = FALSE)
###################################################
## system("cmd", input = "plink --file mydata --maf 0.3 --indep 50 5 2")


###################################################
### code chunk number 9: sequoia.Rnw:491-492 (eval = FALSE)
###################################################
## GenoM <- GenoConvert(InFile = "inputfile_for_sequoia.raw", InFormat="raw")


###################################################
### code chunk number 10: sequoia.Rnw:624-631 (eval = FALSE)
###################################################
## load("Sequoia_output_date.RData")  # if it was saved to disk
## ParOUT$Specs$DummyPrefixFemale <- "D-FEM"
## ParOUT$Specs$DummyPrefixMale <- "D-MALE"
## SeqOUTX <- sequoia(GenoM = Geno,
##                   SeqList = list(Specs = ParOUT$Specs,
##                                  PedigreePar = ParOUT$PedigreePar),
##                   MaxSibIter = 10)


###################################################
### code chunk number 11: sequoia.Rnw:702-703 (eval = FALSE)
###################################################
## save(SeqList, LHdata, Geno, file="Sequoia_output_date.RData")


###################################################
### code chunk number 12: sequoia.Rnw:707-709 (eval = FALSE)
###################################################
## load("Sequoia_output_date.RData")
## # 'SeqList' and 'LHdata' will appear in R environment


###################################################
### code chunk number 13: sequoia.Rnw:714-715 (eval = FALSE)
###################################################
## writeSeq(SeqList, GenoM = Geno, folder=paste("Sequoia_OUT", Sys.Date()))


###################################################
### code chunk number 14: sequoia.Rnw:720-721 (eval = FALSE)
###################################################
## writeSeq(SeqList, OutFormat="xls", file="Sequoia_OUT.xlsx")


###################################################
### code chunk number 15: sequoia.Rnw:724-727 (eval = FALSE)
###################################################
## library(xlsx)
## write.xlsx(Geno, file = "Sequoia_OUT.xlsx", sheetName="Genotypes",
##       col.names=FALSE, row.names=TRUE, append=TRUE, showNA=FALSE)


###################################################
### code chunk number 16: sequoia.Rnw:780-801 (eval = FALSE)
###################################################
## data(LH_HSg5, Ped_HSg5)
## 
## GM <- SimGeno(Ped = Ped_HSg5, nSnp = 100, SnpError = 0.01)
## #
## LH <- LH_HSg5
## LH$BirthYear[sample.int(nrow(LH), round(nrow(LH)*0.2))] <- NA
## LH$Sex[sample.int(nrow(LH), round(nrow(LH)*0.2))] <- NA
## #
## # run sequoia, parentage assignment only
## SeqX <- sequoia(GenoM = GM,  LifeHistData = LH, MaxSibIter = 0,
##                 Err = 0.005)  # assumed genotyping error rate
## #
## # check if there were any assignment errors
## compPar <-  PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqX$PedigreePar)
## compPar$Counts["GG",, ]  # GG = Genotyped offspring, Genotyped parent
## #         dam sire
## # Total    550  553
## # Match    494  478
## # Mismatch  11   12
## # P1only    45   63
## # P2only     0    0


###################################################
### code chunk number 17: sequoia.Rnw:806-814 (eval = FALSE)
###################################################
## print(head(compPar$Mismatch[, c(1:5, 9:12)]), row.names=F)  # subset columns for easier viewing
##  #    id  dam.1 sire.1  dam.2 sire.2 id.dam.cat id.sire.cat dam.class sire.class
##  # a02015 a01002 b01151 a01001 b01151         GG          GG  Mismatch      Match
##  # a02082 a01070 b01069 a02084 b01072         GG          GG  Mismatch   Mismatch
##  # a02084 a01070 b01069 a01071 b01072         GG          GG  Mismatch   Mismatch
##  # a03052 a02112 b02158 a02110   <NA>         GG          GD  Mismatch     P1only
##  # a03112 a02107 b02187 a03111 b02187         GG          GG  Mismatch      Match
##  # a03157 a02051 b02061 a04122   <NA>         GG          GD  Mismatch     P1only


###################################################
### code chunk number 18: sequoia.Rnw:821-831 (eval = FALSE)
###################################################
## SeqXF <- sequoia(GenoM = GM,  LifeHistData = LH, MaxSibIter = 10,
##                 Err = 0.005)
## comp <-  PedCompare(Ped1 = Ped_HSg5, Ped2 = SeqXF$Pedigree)
## comp$Counts["GG",, ]   # let's look at `GG' again first
## #          dam sire
## # Total    549  550
## # Match    530  526
## # Mismatch  10    8
## # P1only     9   16
## # P2only     0    0


###################################################
### code chunk number 19: sequoia.Rnw:836-843 (eval = FALSE)
###################################################
## comp$Counts["TT",, ]   # Totals
## #          dam sire
## # Total    962  962
## # Match    885  880
## # Mismatch  20   21
## # P1only    55   59
## # P2only     2    2


###################################################
### code chunk number 20: sequoia.Rnw:847-857 (eval = FALSE)
###################################################
## tmp <- apply(comp$Counts, 1:2, mean)
## round(tmp / tmp[,"Total"], 3)
## #    Total Match Mismatch P1only P2only
## # GG     1 0.961    0.016  0.023  0.000
## # GD     1 0.850    0.027  0.123  0.000
## # GT     1 0.918    0.021  0.062  0.000
## # DG     1 0.939    0.045  0.015  0.000
## # DD     1 0.946    0.018  0.036  0.000
## # DT     1 0.943    0.033  0.025  0.000
## # TT     1 0.917    0.021  0.059  0.002


###################################################
### code chunk number 21: sequoia.Rnw:861-865 (eval = FALSE)
###################################################
## print(comp$P2only[, 1:6], row.names=F)
##  #    id dam.1 sire.1  dam.2 sire.2    id.r
##  # M0033  <NA>   <NA> a02116 b02009 nomatch
##  # M0034  <NA>   <NA>  F0022 b03005 nomatch


###################################################
### code chunk number 22: sequoia.Rnw:870-876 (eval = FALSE)
###################################################
## print(comp$MergedPed[which(comp$MergedPed$sire.2 == "M0033"), 1:8], row.names=F)
##  #     id  dam.1 sire.1  dam.2 sire.2   id.r  dam.r  sire.r
##  # a04033 a03177 b03092 a03177  M0033 a04033   <NA> nomatch
##  # a04045 a03083 b03092  F0045  M0033 a04045 a03083 nomatch
##  # b04081 a03165 b03092 a03165  M0033 b04081   <NA> nomatch
##  # b04084 a03165 b03092 a03165  M0033 b04084   <NA> nomatch


###################################################
### code chunk number 23: sequoia.Rnw:879-894 (eval = FALSE)
###################################################
## print(comp$MergedPed[which(comp$MergedPed$sire.1 == "b03092"),
##                      c("id", "sire.1", "sire.2", "sire.r")], row.names=F)
##  #    id sire.1 sire.2  sire.r
##  # a04033 b03092  M0033 nomatch
##  # a04034 b03092  M0027  b03092
##  # a04035 b03092  M0027  b03092
##  # a04036 b03092  M0027  b03092
##  # a04045 b03092  M0033 nomatch
##  # a04046 b03092  M0027  b03092
##  # a04082 b03092  M0027  b03092
##  # b04048 b03092  M0027  b03092
##  # b04081 b03092  M0033 nomatch
##  # b04083 b03092  M0027  b03092
##  # b04084 b03092  M0033 nomatch
##  #  M0010 b03092  M0027  b03092


###################################################
### code chunk number 24: sequoia.Rnw:902-915 (eval = FALSE)
###################################################
## comp.p  <- ComparePairs(Ped1 = Ped_HSg5, Ped2 = SeqXF$Pedigree,
##                         GenBack = 1,
##                         patmat = TRUE)  # distinguish between maternal & paternal relatives
## comp.p
## #      Ped2
## # Ped1       M      P      O     FS    MHS    PHS      U      X
## #   M      530      0      2      0      0      0     11    417
## #   P        0    526      1      0      0      0     16    417
## #   FS       3      0      3   1197     35     22    145    198
## #   MHS      0      0      0      0   1342      0    200    218
## #   PHS      0      0      0      0      0   2758    476    446
## #   U        6      6     12      1     31     13 415419  75064
## #   X        0      0      0      0      0      0      0      0


###################################################
### code chunk number 25: sequoia.Rnw:923-927 (eval = FALSE)
###################################################
## BestConfig <- read.table("Colony/file/file.BestConfig",
##                          header=T, sep="", comment.char="")
## PedCompare(Ped1 = ExistingPedigree,
##            Ped2 = BestConfig)


###################################################
### code chunk number 26: sequoia.Rnw:936-955 (eval = FALSE)
###################################################
## Ped <- PedPolish(SeqOUT$Pedigree[, 1:3],
##                  FillParents = TRUE)  # add fake single parents
## Ped <- merge(Ped, MyLifeHistData, by.x="id", by.y="ID", all.x=TRUE)
## 
## library(kinship2)
## # change sex coding: for kinship2 1=male, 2=female
## Ped.fix <- with(Ped_, kinship2::fixParents(id=id, dadid=sire, momid=dam,
##                                sex=c("female", "male", "unknown")[Sex]))
## # create a pedigree object:
## ped.k <- with(Ped.fix, kinship2::pedigree(id, dadid, momid, sex, missid=0))
## # calculate kinship matrix:
## kin <- kinship(ped.k)
## # turn matrix into dataframe:
## kin.df <- data.frame("IID1" = rep(rownames(kin.M), each=ncol(kin.M)),
##                      "IID2" = rep(colnames(kin.M), times=nrow(kin.M)),
##                      "k.ped" = c(t(kin.M)),
##                      stringsAsFactors=FALSE)
## # calculate pedigree relatedness:
## kin.df$R.ped <- kin.df$k.ped *2


###################################################
### code chunk number 27: sequoia.Rnw:962-992 (eval = FALSE)
###################################################
## Rel.snp <- read.table("GT.grm.gz")
## Rel.id <- read.table("GT.grm.id", stringsAsFactors=FALSE)
## Rel.snp[,1] <- as.character(factor(Rel.snp[,1], labels=Rel.id[,2]))
## Rel.snp[,2] <- as.character(factor(Rel.snp[,2], labels=Rel.id[,2]))
## names(Rel.snp) <- c("IID1", "IID2", "SNPS", "R.SNP")
## Rel.snp <- Rel.snp[Rel.snp$IID1 != Rel.snp$IID2,]
## #
## library(data.table)
## Rel.gt <- merge(data.table(Rel.snp[,c(1,2,4)], key=c("IID1", "IID2")),
##                 data.table(kin.df, key=c("IID1", "IID2")), all.x=TRUE)
## Rel.gt <- as.data.frame(Rel.gt)  # turn back into regular dataframe
## rm(PedStats, Rel.snp, kin.df)   # clean up large dataframes
## #
## round(cor(Rel.gt[, c("R.SNP", "R.ped")], use="pairwise.complete"),4)
## #
## library(hexbin)
## ColF <- function(n) rev(rainbow(n, start=0, end=4/6,
##                                 s=seq(.9,.6,length.out=n),v=.8))
## hexbinplot(Rel.gt$R.SNP~Rel.gt$R.ped, xbins=100, aspect=1, maxcnt=10^6.5,
##            trans=log10,inv=function(x) 10^x, colorcut=seq(0,1,length=14),
##            xlab="Pedigree relatedness", ylab="Genomic relatedness",
##            xlim=c(-.1,.9), ylim=c(-.1, .9), colramp=ColF, colorkey = TRUE)
## #
## # or, if you want to add a diagonal line to the plot:
## hb <- hexbin(Rel.gt$R.SNP~Rel.gt$R.ped, xbins=100,
##              xbnds=c(-.1,.9), ybnds=c(-.1, .9),
##              xlab="Pedigree relatedness", ylab="Genomic relatedness")
## hbp <- plot(hb, maxcnt=10^6.5, colramp=ColF, trans=log10,
##            inv=function(x) 10^x, colorcut=seq(0,1,length=14))
## hexVP.abline(hbp$plot.vp, a=0, b=1)


###################################################
### code chunk number 28: sequoia.Rnw:1046-1062 (eval = FALSE)
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


