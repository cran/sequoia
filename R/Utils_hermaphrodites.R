herm_clone_LH <- function(LH, herm.suf=c("f", "m")) {
  hermLH.1 <- LH[which(LH$Sex==4), ]
  hermLH.1$ID <- paste(hermLH.1$ID, herm.suf[1], sep="_")
  hermLH.1$Sex = 1
  hermLH.2 <- LH[which(LH$Sex==4), ]
  hermLH.2$ID <- paste(hermLH.2$ID, herm.suf[2], sep="_")
  hermLH.2$Sex = 2
  unique(rbind(LH, hermLH.1, hermLH.2))
}


#=============================================================
herm_clone_Geno <- function(Geno, LH, herm.suf=c("f", "m")) {
  hermID <- LH$ID[which(LH$Sex==4 & LH$ID %in% rownames(Geno))]
  if (length(hermID)>0) {  # 0: no genotyped hermaprhodites / already duplicated
    herm.Geno.2 <- Geno[hermID, ]
    rownames(herm.Geno.2) <- paste(hermID, herm.suf[2], sep="_")
    these <- which(rownames(Geno) %in% hermID)
    rownames(Geno)[these] <- paste(rownames(Geno)[these], herm.suf[1], sep="_")
    Geno <- rbind(Geno, herm.Geno.2)
  }
  Geno
}


#=============================================================
herm_clone_Ped <- function(Ped, LH, herm.suf=c("f", "m")) {
  Ped <- AddParPed(Ped)
  hermID <- as.character(LH$ID[which(LH$Sex==4 & LH$ID %in% Ped$id)])
  these <- which(Ped$id %in% hermID)
  if (length(these)>0) {
    Ped.herm2 <- Ped[these, ]
    Ped$id[these] <- paste(Ped$id[these], herm.suf[1], sep="_")
    Ped.herm2$id <- paste(Ped.herm2$id, herm.suf[2], sep="_")
    Ped <- rbind(Ped, Ped.herm2)
    these.dams <- which(Ped$dam %in% hermID)
    if (length(these.dams)>0) {
      Ped$dam[these.dams] <- paste(Ped$dam[these.dams], herm.suf[1], sep="_")
    }
    these.sires <- which(Ped$sire %in% hermID)
    if (length(these.sires)>0) {
      Ped$sire[these.sires] <- paste(Ped$sire[these.sires], herm.suf[2], sep="_")
    }
  }
  Ped
}


#=============================================================
herm_unclone_Ped <- function(Ped, LH, herm.suf=c("f", "m")) {
  hermID <- as.character(LH$ID[which(LH$Sex==4)])
  names(Ped)[1:3] <- c("id", "dam", "sire")
  hermPed.1 <- Ped[which(substr(Ped$id, nchar(Ped$id)-1, nchar(Ped$id)) == paste0("_",herm.suf[1]) &
                           chop(Ped$id, herm.suf[1]) %in% hermID), ]
  hermPed.2 <- Ped[which(substr(Ped$id, nchar(Ped$id)-1, nchar(Ped$id)) == paste0("_",herm.suf[2]) &
                           chop(Ped$id, herm.suf[2]) %in% hermID), ]
  hermPed.1$id.orig <- chop(hermPed.1$id, herm.suf[1])
  hermPed.2$id.orig <- chop(hermPed.2$id, herm.suf[2])
  hermPed.b <- merge(hermPed.1, hermPed.2, by="id.orig", all=TRUE)
  if(any(!eqv(hermPed.b$dam.x, hermPed.b$dam.y, xNA=TRUE)) |
      any(!eqv(hermPed.b$sire.x, hermPed.b$sire.y, xNA=TRUE))) {
    warning("Different parents assigned to hermaphrodite in-silico clones")
    PedOUT <- Ped
  } else {
    Hped <- hermPed.1
    Hped <- Hped[, names(Hped)!="id"]
    names(Hped)[names(Hped)=="id.orig"] <- "id"
    Hped <- Hped[, names(Ped)]
    PedOUT <- rbind(Ped[!Ped$id %in% hermPed.1$id & !Ped$id %in% hermPed.2$id, ],
                    Hped)
    PedOUT$dam <- chop(PedOUT$dam, herm.suf[1])
    PedOUT$sire <- chop(PedOUT$sire, herm.suf[2])
  }
  PedOUT
}


#=============================================================
# remove duplicates (e.g. after simulating)
herm_unclone_Geno <- function(Geno, LH, herm.suf=c("f", "m")) {
  hermID <- as.character(LH$ID[which(LH$Sex==4)])
  Geno <- Geno[!rownames(Geno) %in% paste(hermID, herm.suf[2], sep="_"), ]
  these <- which(rownames(Geno) %in% paste(hermID, herm.suf[1], sep="_"))
  rownames(Geno)[these] <- chop(rownames(Geno)[these], herm.suf[1])
  Geno
}

#=============================================================
herm_unclone_MaybeRel <- function(MR, Ped, LH, herm.suf=c("f", "m")) {
  hermID <- as.character(LH$ID[which(LH$Sex==4)])
  H1 <- substr(MR$ID1, nchar(MR$ID1)-1, nchar(MR$ID1)) %in% paste0("_", herm.suf) &
                           chop(MR$ID1, herm.suf[MR$Sex1]) %in% hermID
  H2 <- substr(MR$ID2, nchar(MR$ID2)-1, nchar(MR$ID2)) %in% paste0("_", herm.suf) &
                           chop(MR$ID2, herm.suf[MR$Sex2]) %in% hermID
  MR$ID1 <- ifelse(H1, chop(MR$ID1, suf=herm.suf[MR$Sex1]), MR$ID1)
  MR$ID2 <- ifelse(H2, chop(MR$ID2, suf=herm.suf[MR$Sex2]), MR$ID2)
  MR$Sex1[H1] <- 4
  MR$Sex2[H2] <- 4
  tmp <- merge(MR, Ped[, 1:3], by.x="ID1", by.y="id")
  if(any(tmp$ID2==tmp$dam | tmp$ID2==tmp$sire, na.rm=TRUE)) {
    MR <- with(tmp, tmp[-which(ID2==dam | ID2==sire), names(MR)])
  } else {
    MR <- tmp[, names(MR)]
  }
  MR <- unique(MR[MR$ID1 != MR$ID2, ])
  MR
}

#=============================================================
herm_unclone_Trios <- function(trios, LH, herm.suf=c("f", "m")) {
  hermID <- as.character(LH$ID[which(LH$Sex==4)])
  for (x in 1:3) {
    trios[,x] <- chop(trios[,x], suf=herm.suf[1])
    trios[,x] <- chop(trios[,x], suf=herm.suf[2])
  }
  trios[!duplicated(trios[,1]), ]
}

#=============================================================
# chop suffix from end of character string
chop <- function(x, suf, sep="_") {
   suf <- paste0(sep, suf)
   ifelse(substr(x, start=nchar(x)-nchar(suf)+1, stop=nchar(x)) == suf,
          substr(x, start=1, stop=nchar(x)-nchar(suf)), x)
}


#=============================================================

