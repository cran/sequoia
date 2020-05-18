#' @title Summarise sequoia output or pedigree
#'
#' @description Number of assigned parents and grandparents and sibship sizes,
#'   split by genotyped, dummy, and 'observed'.
#'
#' @param SeqList the list returned by \code{\link{sequoia}}. Only elements
#'   'Pedigree' or 'PedigreePar' and 'AgePriors' are used.
#' @param Pedigree  Dataframe, pedigree with the first three columns being id - dam -
#'   sire. Column names are ignored, as are additional columns.
#' @param DumPrefix character vector of length 2 with prefixes for dummy dams
#'   (mothers) and sires (fathers). Will be read from \code{SeqList}'s 'Specs'
#'   if provided. Used to distinguish between dummies and non-dummies.
#' @param SNPd character vector with ids of SNP genotyped individuals. Only when
#'   \code{Pedigree} is provided instead of \code{SeqList}, then used to distinguish
#'   between genetically assigned parents and 'observed' parents (e.g. observed
#'   in the field, or assigned previously using microsatellites). Will be read
#'   from \code{SeqList}'s 'PedigreePar' if provided.
#' @param Plot  Show barplots and histograms of the results, as well as of the
#'   parental LLRs, Mendelian errors, and agepriors, if present.
#' @param Panels  character vector with panel(s) to plot. Choose from 'all',
#'   'G.parents' (parents of genotyped individuals), 'D.parents' (parents of
#'   dummy individuals), 'sibships' (distribution of sibship sizes), 'LLR'
#'   (log10-likelihood ratio parent/otherwise related), 'OH' (count of opposite
#'   homozygote SNPs).
#'
#' @return A list with the following elements:
#'   \item{PedSummary}{a 1-column dataframe with basic summary statistics, as
#'   used to be returned by \code{Pedantics}' \code{pedStatSummary} (now
#'   archived on CRAN)}
#'   \item{ParentCount}{a 2x3x2x4 array with the number of assigned parents,
#'   split by D1: genotyped vs dummy individuals; D2: female, male and
#'   unknown-sex individuals; D3: dams vs sires; D4: genotyped, dummy, observed
#'   vs no parent}
#'   \item{GPCount}{a 4x4 matrix with for all genotyped individuals the number
#'   of assigned grandparents, split by D1: Maternal grandmother, maternal
#'   grandfather, paternal grandmother, paternal grandfather; D2: genotyped,
#'   dummy, observed vs no grandparent}
#'   \item{SibSize}{a list with as first element a table of maternal sibship
#'   sizes, and as second element a table of paternal sibship sizes. Each table
#'   is a matrix with a number of rows equal to the  maximum sibship size, and 3
#'   columns, splitting by the type of parent: genotyped, dummy, or observed.}
#'
#' @seealso \code{\link{sequoia}} for pedigree reconstruction;
#'   \code{\link{CalcOHLLR}} to (re-)calculate opposite homozygosity & parental
#'   LLR; \code{\link{PlotAgePrior}} to visualise just the ageprior.
#'
#' @examples
#' \dontrun{
#' data(SimGeno_example, LH_HSg5, package="sequoia")
#' SeqOUT <- sequoia(GenoM = SimGeno_example,
#'                    LifeHistData = LH_HSg5, MaxSibIter = 10)
#' Ped_example <- SeqOUT$Pedigree
#' Ped_example$dam[1:20] <- paste0("Mum", 1:20)  # some field mums
#' SummarySeq(SeqOUT, Pedigree=Ped_example)
#' }
#'
#' @export

SummarySeq <- function(SeqList = NULL,
                       Pedigree = NULL,
                       DumPrefix = c("F0", "M0"),
                       SNPd = NULL,
                       Plot = TRUE,
                       Panels = "all")
{
  if(!is.null(Pedigree)) {
    PedIN <- Pedigree
  } else if ("Pedigree" %in% names(SeqList)) {
    PedIN <- SeqList$Pedigree
  } else if ("PedigreePar"  %in% names(SeqList)) {
    PedIN <- SeqList$PedigreePar
  } else if (is.data.frame(SeqList) & any(c("id", "ID") %in% names(SeqList))) {
    PedIN <- SeqList
    SeqList <- NULL
  } else {
    stop("Please provide dataframe 'Pedigree' or SeqList with element Pedigree or PedigreePar")
  }

  # check input
  Ped <- PedPolish(PedIN, ZeroToNA=TRUE)
  Ped$Sex <- ifelse(Ped$id %in% Ped$dam, "Female",
                ifelse(Ped$id %in% Ped$sire, "Male", "Unknown"))
  Ped$Sex <- factor(Ped$Sex, levels=c("Female", "Male", "Unknown"))

  if ("Specs" %in% names(SeqList)) {
    DumPrefix <- c(SeqList$Specs$DummyPrefixFemale, SeqList$Specs$DummyPrefixMale)
    if (nchar(DumPrefix[1])==1)  DumPrefix <- paste0(DumPrefix, "0")
  } else {
    if (length(DumPrefix)!=2)  stop("DumPrefix should be a character vector of length 2")
  }
  if (is.null(SNPd) & "PedigreePar" %in% names(SeqList)) {
    SNPd <- SeqList$PedigreePar$id
  }

  for (x in c("id", "dam", "sire")) {
    Ped[, paste0("GDO.", x)] <- getGDO(Ped[, x], SNPd, DumPrefix)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~
  ParentCount <- array(dim=c(2,3,2,4),
                       dimnames=list(c("G", "D"),   # type of indiv
                                     c("Female", "Male", "Unknown"),   # sex of indiv
                                     c("Dam", "Sire"),   # sex of parent
                                     c("Genotyped", "Dummy", "Observed", "None")))  # type of parent
  ParentCount[,,"Dam",] <- with(Ped, table(GDO.id, Sex, GDO.dam))[c("Genotyped", "Dummy"),,]
  ParentCount[,,"Sire",] <- with(Ped, table(GDO.id, Sex, GDO.sire))[c("Genotyped", "Dummy"),,]

  #~~~~~~~~~~~~~~~~~~~~~~~~~
  PedGP <- merge(Ped[Ped$GDO.id=="Genotyped", 1:3],
               setNames(Ped[, c("id","dam","sire","GDO.dam","GDO.sire")],
                        c("dam","MGM","MGF","GDO.MGM","GDO.MGF")), all.x=TRUE)
  PedGP <- merge(PedGP,
                 setNames(Ped[, c("id","dam","sire","GDO.dam","GDO.sire")],
                          c("sire","PGM","PGF","GDO.PGM","GDO.PGF")), all.x=TRUE)
  for (x in paste0("GDO.",c("MGM","MGF","PGM","PGF"))) {
    PedGP[is.na(PedGP[,x]), x] <- "None"
  }
  GPCount <- matrix(NA,4,4,
                    dimnames = list(c("MGM", "MGF", "PGM", "PGF"),
                                    c("Genotyped", "Dummy", "Observed", "None")))
  for (x in c("MGM","MGF","PGM","PGF")) {
    GPCount[x,] <- table(PedGP[, paste0("GDO.",x)])
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~
  Ped$ParentPair <- with(Ped, ifelse(is.na(dam) | is.na(sire), NA,
                                     paste0(dam, "_-_", sire)))
  MaxSibSize <- c("mat" = ifelse(any(!is.na(Ped$dam)), max(table(Ped$dam)), 0),
                  "pat" = ifelse(any(!is.na(Ped$sire)), max(table(Ped$sire)), 0),
                  "full" = ifelse(any(!is.na(Ped$ParentPair)), max(table(Ped$ParentPair)), 0))
  SibSize <- list()
  for (p in c("mat", "pat")) {
    if (MaxSibSize[p] == 0)  next
    SibSize[[p]] <- matrix(NA, MaxSibSize[p], 3,
                           dimnames=list(1:MaxSibSize[p],
                                         c("Genotyped", "Dummy", "Observed")))
    for (x in c("Genotyped", "Dummy", "Observed")) {
      if (p=="mat") {
        tbl <- with(Ped, table(dam[GDO.dam==x]))
      } else if (p=="pat") {
        tbl <- with(Ped, table(sire[GDO.sire==x]))
      }
      SibSize[[p]][, x] <- table(factor(tbl, levels=1:MaxSibSize[p]))
    }
  }
  SibSize[["full"]] <- table(factor(table(Ped$ParentPair), levels=1:MaxSibSize["full"]))

  #~~~~~~~~~~~~~~~~~~~~~~~~~
  SibCount <- setNames(rep(NA,3), c("full", "mat", "pat"))
  for (x in c("full", "mat", "pat")) {
    p <- c(full = "ParentPair", mat = "dam", pat = "sire")[x]
    tbl <- table(table(Ped[, p]))
    SibS <- as.numeric(rownames(tbl))
    SibCount[x] <- sum(tbl * SibS * (SibS -1) / 2)
  }
  PedSummary <- c("records" = nrow(Ped),
                  "maternities" = sum(!is.na(Ped$dam)),
                  "paternities" = sum(!is.na(Ped$sire)),
                  "full sibs" = SibCount["full"],
                  "maternal half sibs" = SibCount["mat"] - SibCount["full"],
                  "paternal half sibs" = SibCount["pat"] - SibCount["full"],
                  rowSums(GPCount[, 1:3]) )
  names(PedSummary)[7:10] <- c("maternal grandmothers",
                               "maternal grandfathers",
                               "paternal grandmothers",
                               "paternal grandfathers")
  PedSummary["maximum pedigree depth"] <- max(getGenerations(Ped)[,1])
  PedSummary["founders"] <- sum(is.na(Ped$dam) & is.na(Ped$sire))
  #~~~~~~~~~~~~~~~~~~~~~~~~~

  SummaryOUT <- list(PedSummary = as.data.frame(PedSummary),
                     ParentCount = ParentCount,
                     GPCount = GPCount,
                     SibSize = SibSize)

  if (Plot) {
    plot.seq(SummaryOUT, Ped, Panels)
    if (!is.null(SeqList) & any(Panels=="all")) {
      inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
			if ("AgePriorExtra" %in% names(SeqList)) {
				PlotAgePrior(SeqList$AgePriorExtra)
			} else {
				PlotAgePrior(SeqList$AgePriors)
			}
    }
  }

  invisible( SummaryOUT )
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Plot summary overview of sequoia output
#'
#' @description visualise the numbers of assigned parents, sibship sizes, and parental LLRs
#'
#' @param SeqSum  list output from \code{\link{SummarySeq}}
#' @param Pedigree dataframe with at least id, dam and sire in columns 1-3,
#'   respectively. If columns with parental LLRs and/or Mendelian errors are
#'   present, these will be plotted as well.
#' @param Panels  character vector with panel(s) to plot. Choose from 'all',
#'   'G.parents', 'D.parents', 'sibships', 'LLR', 'OH'.
#'
#' @importFrom graphics par barplot hist text axis mtext
#'
#' @keywords internal

plot.seq <- function(SeqSum, Pedigree, Panels="all")
{
  col.damsire <- matrix(c("darkred", "firebrick2", "pink", "lightgrey",
                         "darkblue", "dodgerblue", "lightblue","lightgrey"), 4,2)
  AllPanels <- c('G.parents', 'D.parents', 'sibships', 'LLR', 'OH')
  if (Panels[1]=="all") {
    Panels <- AllPanels
  } else if (!all(Panels %in% AllPanels)) {
    stop("Invalid value for 'Panels'")
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # barplots proportion of parents genotyped / dummy / none
  oldpar <- par(no.readonly = TRUE)
  par(mai=c(.9, 1.8, 1.4,.3), mfrow=c(1,1))

  ParCount.G <- apply(SeqSum$ParentCount["G",,,], 2:3, sum)
  N.G <- sum(ParCount.G[1,])
  GPCount <- SeqSum$GPCount

  # genotyped individuals
	if ('G.parents' %in% Panels) {
		bp <- barplot(t(rbind(ParCount.G,
													GPCount)[6:1,]), horiz=TRUE, las=1,
									space=c(rep(.2,4), 1,.2),
									main="No. parents and grandparents \nassigned to genotyped individuals",
									xlab = "No. genotyped individuals",
									names.arg=c("Pat. grandfather", "Pat. grandmother",
															"Mat. grandfather","Mat. grandmother",
															"Sire (father)", "Dam (mother)"))
		axis(side=3, at=c(0:10)*N.G/10,
				 labels=paste0(seq(0,100,10),"%"), col="darkgrey", col.axis="darkgrey")
		abline(v=c(0:10)*N.G/10, col="grey", lty=3, xpd=FALSE)

		barplot(t(ParCount.G["Dam",,drop=FALSE]), horiz=TRUE, space=7, las=1,
						col=col.damsire[,1], add=TRUE, axes=F, names.arg="")
		barplot(t(ParCount.G["Sire",,drop=FALSE]), horiz=TRUE, space=5.8, las=1,
						col=col.damsire[,2], add=TRUE, axes=F,names.arg="")

		for (p in 1:2) {
			for (x in 1:4) {
				if (ParCount.G[p, x]>0) {
					rot <- ifelse(ParCount.G[p, x]/N.G < 0.05, TRUE, FALSE)
					xx <- ifelse(x==1, ParCount.G[p, x]/2,
											 ParCount.G[p, x]/2 + sum(ParCount.G[p, 1:(x-1)]))
					text(xx, bp[7-p], colnames(ParCount.G)[x], col=ifelse(x==1, 0, 1),
							 srt=ifelse(rot, 90, 0), cex=ifelse(rot, 0.8, 1))
				}
			}
		}

		barplot(t(GPCount[4:1,]), horiz=TRUE, las=1, space=.2,
						add=TRUE, axes=F, names.arg=rep("",4))
		for (g in 1:4) {
			for (x in 1:4) {
				if (GPCount[g, x]>0) {
					rot <- ifelse(GPCount[g, x]/N.G < 0.05, TRUE, FALSE)
					xx <- ifelse(x==1, GPCount[g, x]/2,
											 GPCount[g, x]/2 + sum(GPCount[g, 1:(x-1)]))
					text(xx, rev(bp[1:4])[g], colnames(GPCount)[x], col=ifelse(x==1, 0, 1),
							 srt=ifelse(rot, 90, 0), cex=ifelse(rot, 0.8, 1))
				}
			}
		}
	}

  #~~~~~~~~~~~~~~~~~
  # dummy individuals
  ParCount.D <- SeqSum$ParentCount["D",1:2,,]
  if (sum(ParCount.D)>0 & 'D.parents' %in% Panels) {
		if ('G.parents' %in% Panels) {
			inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
		}
    bp <- barplot(t(rbind(ParCount.D["Female",,],
                          ParCount.D["Male",,])[4:1,]), horiz=TRUE, las=1, space=c(.2,.2,1,.2),
                  main="No. parents assigned to \ndummy individuals",
                  xlab = "No. dummy individuals",
                  names.arg=rep(c("Sire", "Dam"), 2))
    mtext(c("Dummy females","Dummy males"), side=2, at=c(4.5,1.3), line=3, las=1, cex=1.0)
    barplot(t(ParCount.D[2:1,"Dam",]), horiz=TRUE, space=c(1.4,2.2), las=1,
            col=col.damsire[,1], add=TRUE, axes=F,
            names.arg=rep("",2))
    barplot(t(ParCount.D[2:1,"Sire",]), horiz=TRUE, space=c(.2,2.2), las=1,
            col=col.damsire[,2], add=TRUE, axes=F,
            names.arg=rep("",2))
    bpM <- matrix(rev(bp),2,2)
    maxN <- max(apply(ParCount.D, 1:2, sum))
    for (s in 1:2) {
      for (p in 1:2) {
        for (x in 1:4) {
          if (ParCount.D[s,p, x]>0) {
            rot <- ifelse(ParCount.D[s,p, x]/maxN < 0.1, TRUE, FALSE)
            xx <- ifelse(x==1, ParCount.D[s,p, x]/2,
                         ParCount.D[s,p, x]/2 + sum(ParCount.D[s,p, 1:(x-1)]))
            text(xx, bpM[p,s], colnames(ParCount.G)[x], col=ifelse(x==1, 0, 1),
                 srt=ifelse(rot, 90, 0), cex=ifelse(rot, 0.8, 1))
          }
        }
      }
    }
  } else if ('D.parents' %in% Panels & length(Panels) < length(AllPanels)) {
    warning("No 'D.parents' panel, because no dummy individuals", immediate.=TRUE)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # hist sibship sizes
	if ('sibships' %in% Panels) {
		if (any(c("G.parents", "D.parents") %in% Panels)) {
			inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
		}
		np <- par(mai=c(.85, .8,1,.2), mfrow=c(1,3))
		for (p in 1:2) {
			barplot(t(SeqSum$SibSize[[p]]), col=col.damsire[,p], las=1, space=0,
						xlab="Sibship size", ylab="Count", main=paste0(c("M","P")[p], "aternal sibships"))
		}
		barplot(SeqSum$SibSize[["full"]], col="forestgreen", las=1, space=0,
						xlab="Sibship size", ylab="Count", main="Full sibships")  #  (size >1)
	}

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # LLR
  LLRc <- c("LLRdam", "LLRsire", "LLRpair")
  if (any(LLRc %in% names(Pedigree)) & 'LLR' %in% Panels) {
    for (i in 1:3) {
      Pedigree[which(Pedigree[,LLRc[i]]==-999), LLRc[i]] <- NA
      Pedigree[which(Pedigree[,LLRc[i]]==999), LLRc[i]] <- NA
    }
    if (any(!is.na(unlist(Pedigree[,LLRc])))) {
			if (any(c("G.parents", "D.parents", "sibships") %in% Panels)) {
				inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
			}
      np <- par(mfrow=c(1,3), mai=c(.8,.4,.4,.1))
			brks.LLR <- pretty(x=unlist(Pedigree[, LLRc]), n=50)
			Mainz.LLR <- c("LLR dam / not dam", "LLR sire / not sire", "LLR parent pair",
             "Offspring - dam\nopposing homozygotes", "Offspring - sire\nopposing homozygotes ",
             "Offspring - dam - sire\nMendelian errors")
      for (i in 1:3) {
        if (!LLRc[i] %in% names(Pedigree))  next
        if (all(is.na(Pedigree[,LLRc[i]]))) next
        hist(Pedigree[,LLRc[i]], breaks=brks.LLR, col="grey",
             main=Mainz.LLR[i], xlab="Log10 likelihood ratio", ylab="")
      }
    } else if (length(Panels) < length(AllPanels)) {
      warning("No 'LLR' panel, because LLR columns are all NA", immediate.=TRUE)
    }
  } else if ('LLR' %in% Panels & length(Panels) < length(AllPanels)) {
    warning("No 'LLR' panel, because no LLR columns", immediate.=TRUE)
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # OH
	OHc <- c("OHdam", "OHsire", "MEpair")
	if (any(c("OH_Mother", "OH_Father", "ME_Pair") %in% names(Pedigree))) {
		these <- match(c("OH_Mother", "OH_Father", "ME_Pair"), names(Pedigree))
		names(Pedigree)[na.exclude(these)] <- OHc[!is.na(these)]
	}
  if (any(OHc %in% names(Pedigree)) & 'OH' %in% Panels) {
		for (i in 1:3) {
			Pedigree[which(Pedigree[,OHc[i]]==-9), OHc[i]] <- NA
		}
    if (any(!is.na(unlist(Pedigree[,OHc])))) {
			if (length(Panels)>1) {
				inp <- readline(prompt = "Press <Enter> to continue to next plot ...")
			}
			np <- par(mfrow=c(1,3), mai=c(.8,.4,.4,.1))
      Mainz.OH <- c("Offspring - dam\nopposing homozygotes", "Offspring - sire\nopposing homozygotes ",
               "Offspring - dam - sire\nMendelian errors")
      for (i in 1:3) {
        if (!OHc[i] %in% names(Pedigree))  next
        if (all(is.na(Pedigree[,OHc[i]]))) next
        hist(Pedigree[,OHc[i]], breaks=c(0:(max(Pedigree[, OHc[i]],na.rm=TRUE)+2))-.5, col="grey",
             main=Mainz.OH[i], xlab="Count", ylab="")
      }
    } else if (length(Panels) < length(AllPanels)) {
      warning("No 'OH' panel, because OH columns are all NA", immediate.=TRUE)
    }
  } else if ('OH' %in% Panels & length(Panels) < length(AllPanels)) {
      warning("No 'OH' panel, because no OH columns", immediate.=TRUE)
	}

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # total log-likelihood per iteration
  # nIt <- with(SeqOUT, c(length(TotLikPar), length(TotLikSib)))
  # with(SeqOUT, plot(1:nIt[1], TotLikPar, type="b", lwd=2, col="darkgrey",
  #                   xlim=c(1, sum(nIt)), ylim=c(min(TotLikSib), max(TotLikPar)),
  #                   xlab="Iteration", ylab="Total log-likelihood"))
  # with(SeqOUT, lines((nIt[1]-1) + 1:nIt[2], TotLikSib, type="b", lwd=2))

#   TLL <- c(SeqOUT$TotLikParents, SeqOUT$TotLikSib)
#   xv <- c(paste("p", 1:length(SeqOUT$TotLikParents)-1),
#           paste("s", 1:length(SeqOUT$TotLikSib)-1))
#   plot(TLL, type="b", xaxt="n", xlab="Round")
#   axis(1, at=1:length(TLL), labels=xv)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~
  par(oldpar)
}


