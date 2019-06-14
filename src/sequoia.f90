! Author: Jisca Huisman, post doc in Evolutionary biology at 
! University of Edinburgh, United Kingdom
! jisca.huisman@gmail.com
!
! This code is available under GNU General Public License v3
! 
! The program is described in the paper 
! "Pedigree reconstruction from SNP data: 
! Parentage assignment, sibship clustering, and beyond", 
! in Molecular Ecology Resources, 2017 

! ####################################################################
! @@@@   MODULES   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! ####################################################################
module Global 
implicit none

integer :: nInd, nSnp, nIndLH, maxSibSize, MaxMismatch, maxOppHom, &
 nC(2), nAgeClasses, nPairs, BY1, maxAgePO, XP, Complx, Error=0, quiet, AgePhase                             
  integer, parameter :: mxA=32 ! max no. ancestors considered when testing for pedigree loop
 logical :: Hermaphrodites                          
integer, allocatable, dimension(:) :: Sex, BY, PairType, nFS
integer,allocatable,dimension(:,:) :: Genos, AgeDiff, Parent, OppHomM,&
  nS, PairID, FSID
integer, allocatable, dimension(:,:,:) :: SibID, GpID
double precision :: TF, TA, Er, OcA(3,3), AKA2P(3,3,3), OKA2P(3,3,3)
double precision, allocatable, dimension(:) ::  Lind, PairDLLR, AF
double precision, allocatable, dimension(:,:) :: AHWE, OHWE, LLR_O, &
  LR_parent, AgePriorM, CLL
double precision, allocatable, dimension(:,:,:) :: AKAP, OKAP, OKOP, &
  LR_GP, LindG, PHS, PFS, DumBY, LindX
double precision, allocatable, dimension(:,:,:,:) :: DumP
double precision, allocatable, dimension(:,:,:,:,:) :: XPr                                                        

  contains
pure function MaxLL(V)
double precision, intent(IN) :: V(:)
double precision :: MaxLL
MaxLL = 999D0
if (ANY(V < 0 .and. V>-HUGE(0.0D0))) then
  MaxLL = MAXVAL(V, mask = (V<0 .and. V>-HUGE(0.0D0)), DIM=1)
else
  MaxLL = MINVAL(V, mask = (V>-HUGE(0.0D0)), DIM=1)  
  ! 777: can't do; 888: already is; 999: not calc'd
  ! V should not ever be -INF
endif
end function MaxLL

pure function which(V, x)
integer, intent(IN) :: V(:), x
integer :: which
integer :: i

which = 0        
do i = 1, size(V)
    if (V(i) .eq. x) then
        which = i
        exit
    endif
end do

end function which

end module Global

! ####################################################################
! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing
! Made F conformant by Walt Brainerd

! Adapted by J Huisman (j.huisman@ed.ac.uk) to output rank, in order to
! enable sorting of parallel vectors, and changed to decreasing rather 
! than increasing order

module qsort_c_module
implicit none
public :: QsortC
private :: Partition

 contains
recursive subroutine QsortC(A, Rank)
  double precision, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: Rank
  integer :: iq

  if(size(A) > 1) then
   call Partition(A, iq, Rank)
   call QsortC(A(:iq-1), Rank(:iq-1))
   call QsortC(A(iq:), Rank(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker, Rank)
  double precision, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: Rank
  integer, intent(out) :: marker
  integer :: i, j, TmpI
  double precision :: temp
  double precision :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1
  do
   j = j-1
   do
    if (j < 1) exit              
    if (A(j) <= x) exit
    j = j-1
   end do
   i = i+1
   do
    if (i >= size(A)) exit                      
    if (A(i) >= x) exit
    i = i+1
   end do
   if (i < j) then
    ! exchange A(i) and A(j)
    temp = A(i)
    A(i) = A(j)
    A(j) = temp 
    
    TmpI = Rank(i) 
    Rank(i) = Rank(j)
    Rank(j) = TmpI
   elseif (i == j) then
    marker = i+1
    return
   else
    marker = i
    return
   endif
  end do

end subroutine Partition

end module qsort_c_module


! ####################################################################

! @@@@   PROGRAMS   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ####################################################################

subroutine makeped(Ng, SpecsInt, SpecsDbl, GenoFR, &
  SexRF,  BYRF, APRF, parentsRF, LrRF, OhRF, & 
  Nd, DumParRF, DumLrRF, DumBYRF,  totLL, &  
  nDupGenos,DupGenosFR,nMisMFR, DupGenoLR)
use Global
implicit none

integer, intent(IN) :: Ng, SpecsInt(12)
double precision, intent(IN) :: SpecsDbl(3), APRF(9*SpecsInt(6))
integer, intent(IN) :: SexRF(Ng), BYRF(Ng), GenoFR(Ng*SpecsInt(3))
integer, intent(INOUT) :: parentsRF(2*Ng), OhRF(3*Ng), & 
  Nd(2), DumParRF(2*Ng), DumBYRF(3*Ng), &
  nDupGenos, DupGenosFR(2*Ng), nMisMFR(Ng)
double precision, intent(INOUT) :: LrRF(3*Ng), DumLrRF(3*Ng), TotLL(42), DupGenoLR(Ng)
integer :: ParSib, MaxSibIter, i, j,k,l, maybe, topX, x, &
  CalcLLR, UseAge, nAmbMax
integer, allocatable, dimension(:,:,:) :: DumBYmm


call Initiate(Ng, SpecsInt, SpecsDbl, GenoFR, SexRF,  BYRF, APRF, &  ! IN
  parentsRF, Nd, DumParRF, &   ! IN
  ParSib, MaxSibIter, CalcLLR, UseAge, nAmbMax)  ! OUT

 !=========================
if (ParSib == 0) then
  call duplicates(nDupGenos,DupGenosFR, nMisMFR, DupGenoLR)
!=========================
else if (ParSib == 1) then
  call parents(TotLL)
!=========================  
else if (ParSib == 2) then
  call sibships(MaxSibIter, UseAge, TotLL) 
  
  ! dummies
  Nd = nC
  allocate(DumBYmm(3, nInd/2, 2))
!  if(CalcLLR==1) then
    if(quiet<1) call intpr ("Estimating dummy birth years ... ",-1,0,0)
    call DumBYrange(DumBYmm)  ! slow
!  endif
 
 call AtoVi(GpID, 2,nInd/2, nC, DumParRF)
 call AtoVi(DumBYmm, 3,nInd/2, nC, DumBYRF)
! call AtoVi(SibID, maxSibSize, nInd/2, nC, DumOff)
 deallocate(DumBYmm)
endif

!========================= 
OhRF = -9
do i=1,nInd
  if (ParSib == 1) then
    if (Parent(i,1)>0) OhRF(i) = OppHomM(i, Parent(i,1))
    if (Parent(i,2)>0) OhRF(Ng+i) = OppHomM(i, Parent(i,2))
  else
    if (Parent(i,1)>0)  call CalcOH(i, Parent(i,1), OhRF(i))
    if (Parent(i,2)>0)  call CalcOH(i, Parent(i,2), OhRF(Ng+i))
  endif
  if (Parent(i,1)>0 .and. Parent(i,2)>0) then
    call CalcTrioErr(i, OhRF(2*Ng +i))
  endif
enddo

!=========================
allocate(LR_parent(nInd,3))
LR_parent = 999
allocate(LR_GP(3, nInd/2,2))
LR_GP = 999
if(CalcLLR==1 .and. ParSib>0) then 
  if (quiet<1)  call intpr ("Calculating parental LLR ... ", -1, 0, 0)
  call rchkusr()
  call UpdateAllProbs(2) 
  call CalcParentLLR
  call UpdateAllProbs(2) 
endif
parentsRF = 0
LrRF = 999D0
do i=1,nInd 
  parentsRF(i) = Parent(i,1)
  parentsRF(nInd+i) = Parent(i,2)
enddo                                             
do i=1,nInd 
  LrRF(i) = LR_parent(i, 1)
  LrRF(nInd+i) = LR_parent(i, 2)
  LrRF(2*nInd+i) = LR_parent(i, 3)
enddo
if (ParSib == 2)  call AtoVd(LR_GP, 3,nInd/2, nC, DumLrRF)   

 call DeAllocAll

end subroutine makeped

! ####################################################################

subroutine Initiate(Ng, SpecsInt, SpecsDbl, GenoFR, SexRF,  BYRF, APRF, &
  parentsRF, Nd, DumParRF, &
  ParSib, MaxSibIter, CalcLLR, UseAge, nAmbMax)
use Global
implicit none

integer, intent(IN) :: Ng, SpecsInt(12)
double precision, intent(IN) :: SpecsDbl(3), APRF(9*SpecsInt(6))
integer, intent(IN) :: GenoFR(Ng*SpecsInt(3)), SexRF(Ng), BYRF(Ng), &
  parentsRF(2*Ng), Nd(2), DumParRF(2*Ng)
integer, intent(OUT) :: ParSib, MaxSibIter, CalcLLR, UseAge, nAmbMax
integer :: i, j, k, l, ParTmp(2), x, s


! set global parameters & return values
nInd = Ng
ParSib = SpecsInt(1)
MaxSibIter = SpecsInt(2)
nSnp = SpecsInt(3)
MaxMisMatch = SpecsInt(4)
maxSibSize = SpecsInt(5)
nAgeClasses = SpecsInt(6)
 Complx = SpecsInt(7)
! FindMaybe = SpecsInt(8)
 CalcLLR = SpecsInt(9)
quiet = SpecsInt(10)
nAmbMax = SpecsInt(11)
UseAge = SpecsInt(12)
Er = SpecsDbl(1)
TF = SpecsDbl(2)
TA = SpecsDbl(3)

if (Complx == 4) then
  Hermaphrodites = .TRUE.
  Complx = 2
else
  Hermaphrodites = .FALSE.
endif

allocate(Sex(nInd))
allocate(BY(nInd))
Sex = SexRF
BY = BYRF       

if (nAgeClasses>1) then
  allocate(AgePriorM(nAgeClasses, 9))
else 
  allocate(AgePriorM(2, 9))
endif
AgePriorM = 1
k = 0
do i=1,9
  do j=1, nAgeClasses
    k = k+1
    AgePriorM(j,i) =  APRF(k)
  enddo
enddo
if (nAgeClasses==1) then  ! shouldn't be necessary anymore - TODO CHECK
  AgePriorM(2,:) = 0
  AgePriorM(2, 6:8) = 1   ! UA, dam, sire
endif

allocate(Genos(nSnp, nInd))
Genos = -9
j = 0
do l=1,nSnp      
  do i=1, nInd
    j = j+1
    if (GenoFR(j)/=-9) then
      Genos(l,i) = GenoFR(j)+1                                                                                 
    endif
  enddo    
enddo      

nC = 0           
 call PrepData   ! calc agediff, allocate arrays
 call PrecalcProbs
 call UpdateAllProbs(1)

maxOppHom = MaxMismatch - FLOOR(-nSNP * Er)   ! round up to nearest int

 call rchkusr()
!=========================
! prep (pedigree prior or prev. assigned parents)
allocate(CLL(nInd/2,2))
 CLL = 999D0
allocate(nS(nInd/2,2))
ns = 0
allocate(SibID(maxSibSize, nInd/2, 2))
SibID = 0
allocate(DumBY(nAgeClasses +maxAgePO, nInd/2, 2)) 
DumBY = 0

do i=1,nInd
  ParTmp(1) = parentsRF(i)
  ParTmp(2) = parentsRF(nInd + i)
  do k=1,2
    if (ParTmp(k) > 0) then
      Parent(i,k) = ParTmp(k)
      call CalcLind(i)
    else if (ParTmp(k) < 0) then
      s = -ParTmp(k)
      if (nC(k) < s)  nC(k) = s
      if (ns(s,k) == 0) then
        Parent(i,k) = ParTmp(k)
        nS(s,k) = ns(s,k) +1
        SibID(ns(s,k), s, k) = i
        call CalcCLL(s,k)
        call CalcLind(i)
      else
        call NewPar(i, k, ParTmp(k))
      endif
    endif
  enddo
enddo

!nC = Nd
! reverse  AtoVi(GpID, 2,nInd/2, nC, DumParRF)
do k=1,2
  do s=1, nC(k)
    do i=1,2
      x = (k-1)*2*nInd/2 + (i-1)*nC(1) + s
      GpID(j, s, k) = DumParRF(x)
    enddo
  enddo
enddo

do i=1, nInd
  if (Sex(i)==3) then
    if (ANY(Parent(:,1) == i)) then
      Sex(i) = 1
    else if (ANY(Parent(:,2) == i)) then
      Sex(i) = 2
    endif                                     
  endif
enddo    

! find current FS 
do i=1,nInd-1
  do j=i+1,nInd
    if (Parent(i,1)==Parent(j,1) .and. Parent(i,1)/=0 .and. &
      Parent(i,2)==Parent(j,2) .and. Parent(i,2)/=0) then
      call MakeFS(i, j)
    endif
  enddo
enddo

call UpdateAllProbs(3)

end subroutine Initiate

! ####################################################################

subroutine duplicates(nDupGenos,DupGenosFR,nMisMFR, DupGenoLR)
use Global
implicit none

integer, intent(INOUT) :: nDupGenos, DupGenosFR(2*nInd), nMisMFR(nInd)
double precision, intent(INOUT) :: DupGenoLR(nInd)
integer :: i, j, l, Match, DupGenos(nInd,2), nMisMatch(nInd), CountMismatch
double precision :: LL(7), LLtmp(2), LLX(7), LRdup(nInd)

!====================
do i=1,nInd-1
  do j=i+1, nInd
    Match=1
    CountMismatch=0
    do l=1, nSnp
      if (Genos(l,i)==-9 .or. Genos(l,j)==-9) cycle
      if (Genos(l,i) /= Genos(l,j)) then
        CountMismatch=CountMismatch+1
        if (CountMismatch > MaxMismatch) then
          Match=0
          exit
        endif
      endif
    enddo
    if (Match==1) then
      nDupGenos = nDupGenos + 1
      DupGenos(nDupGenos,1) = i
      DupGenos(nDupGenos,2) = j
      nMisMatch(nDupGenos) = CountMismatch
    endif
    if (nDupGenos==nInd) then
        if(quiet<1) call rwarn("reached max for duplicates") 
      exit
    endif
    if (nDupGenos==nInd) exit
  enddo
  if (nDupGenos==nInd) exit
enddo
 
!================
! call intpr("n dup: ", -1, nDupGenos, 1)
if (nDupGenos > 0) then
  maxOppHom = MaxMismatch - FLOOR(-nSNP * Er)   ! round up to nearest int

  do i=1, nDupGenos
    LLtmp = 999
    call PairSelf(DupGenos(i,1),DupGenos(i,2), LLtmp(1))
    call CheckPair(DupGenos(i,1),DupGenos(i,2),3,1,LL, LLX)
    LLtmp(2) = MaxLL(LL)
    if (LLtmp(1) < 0 .and. LLtmp(2)<0) then
      LRdup(i) = LLtmp(1) - LLtmp(2)
    else 
      LRdup(i) = 111
    endif
  enddo
endif  

! ##########################
! copy to vectors to send to R

DupGenosFR = 0
nMisMFR = 0
DupGenoLR = 0
do i=1,nDupGenos
  DupGenosFR(i) = DupGenos(i,1)
  DupGenosFR(nInd+i) = DupGenos(i,2)
  nMisMFR(i) = nMisMatch(i)
  DupGenoLR(i) = LRdup(i)
enddo        

 end subroutine duplicates

! ####################################################################

subroutine findambig(Ng, SpecsInt, SpecsDbl, GenoFR, &
    SexRF, BYRF, APRF, parentsRF, Nd, DumParRF, &
    nAmb, AmbigID, AmbigRel, AmbigLR, AmbigOH, &
    ntrio, trioIDs, trioLR, trioOH)
use Global
implicit none

integer, intent(IN) :: Ng, SpecsInt(12)
double precision, intent(IN) :: SpecsDbl(3), APRF(9*SpecsInt(6))
integer, intent(IN) :: SexRF(Ng), BYRF(Ng), GenoFR(Ng*SpecsInt(3))
integer, intent(INOUT) :: parentsRF(2*Ng), Nd(2), DumParRF(2*Ng), &
  nAmb, AmbigID(2*SpecsInt(11)), AmbigRel(2*SpecsInt(11)), &
  AmbigOH(SpecsInt(11)), ntrio, trioIDs(3*Ng), trioOH(3*Ng)
double precision, intent(INOUT) :: AmbigLR(2*SpecsInt(11)), trioLR(3*Ng) 
integer :: ParSib, MaxSibIter, CalcLLR, UseAge, nAmbMax, &
  i, j, k, x, topX, Anc(2,mxA), ADX, maybe, OH, Lboth
double precision :: LL(7), LLtmp(7,3), dLL(2), LRR(3), LLX(7)

!open(unit=101, file="log.txt", status="unknown")
!write(101, '("Welcome to FindAmbig")') 
!close(101)

call Initiate(Ng, SpecsInt, SpecsDbl, GenoFR, SexRF,  BYRF, APRF, &  ! IN
  parentsRF, Nd, DumParRF, &   ! IN
  ParSib, MaxSibIter, CalcLLR, UseAge, nAmbMax)  ! OUT

nAmb = 0
AmbigID = 0
AmbigLR = 999D0
AmbigOH = -9
AmbigRel = 0   
 
ntrio = 0
trioLR = 999D0 
trioOH = -9       

if(quiet<1) then
  if (ParSib ==1) then
    call intpr("Checking for non-assigned Parent-Offspring pairs ... ",&
      -1, 0, 0)
  else if (ParSib ==2) then
    call intpr ( "Checking for non-assigned relatives ... ", -1, 0, 0)
  endif
endif

allocate(OppHomM(nInd, nInd))
OppHomM = 999 
allocate(LLR_O(nInd, nInd))
LLR_O = 999D0
call CalcOppHom   ! also checks no. SNPs typed in both
if (ParSib == 1) then
  do i=1, nInd-1
    do j=i+1,nInd 
      if (OppHomM(i,j) > maxOppHom) cycle
        call CalcPO(i, j, LLR_O(i,j))  ! LLR PO/U
        call CalcPO(j, i, LLR_O(j,i))
    enddo
  enddo
endif

do i=1,nInd-1
  if (Error/=0) return
  if (MODULO(i,500)==0)   call rchkusr()
  if (quiet<1 .and. nInd>1500) then 
    if (MODULO(i,500)==0) call intpr (" ", 1, i, 1)
  endif
  do j=i+1,nInd
    Lboth = COUNT(Genos(:,i)/=-9 .and. Genos(:,j)/=-9)    
    if (Lboth < nSnp/2.0)   cycle   ! >1/2th of markers missing
    if (ANY(Parent(i,:)==j) .or. ANY(Parent(j,:)==i)) cycle  ! PO
    if (ALL(Parent(i,:)/=0)) then
      if (Parent(i,1)==Parent(j,1) .and. Parent(i,2)==Parent(j,2)) cycle  ! FS
      call GetAncest(i,1,Anc)
      if (ANY(Anc(:,3)==j) .and. ANY(Anc(:,4)==j)) cycle  ! double GP
    endif
    if (Parent(j,1)/=0 .and. Parent(j,2)/=0) then
      call GetAncest(j,1,Anc)
      if (ANY(Anc(:,3)==i) .and. ANY(Anc(:,4)==i)) cycle  ! double GP
    endif
    
    LL = 999D0
    topX = 0
    LLtmp = 999D0
    if (ParSib <2 .or. All(Parent(i,:)/=0) .and. ALL(Parent(j,:)/=0)) then   ! check if they're not PO only
      if (OppHomM(i,j) > maxOppHom .or. OppHomM(i,j)<0)  cycle
      if (LLR_O(i,j)==999 .and. LLR_O(j,i)==999)  cycle
      if (ParSib < 2 .and. MaxLL((/LLR_O(i,j), LLR_O(j,i)/)) < TA)  cycle
      if (ParSib ==2 .and. MaxLL((/LLR_O(i,j), LLR_O(j,i)/)) < 2*TA)  cycle 
    endif
    if (ParSib < 2) then
      ADX = AgeDiff(i,j)
      AgeDiff(i,j) = 999
      AgeDiff(j,i) = 999
      call CheckPair(i, j, Sex(j), 1, LLtmp(:,1), LLX)
      call CheckPair(j, i, Sex(i), 1, LLtmp(:,2), LLX)
      do k=1,7  
        LL(k) = MaxLL(LLtmp(k,1:2)) 
      enddo
      call BestRel2(LL, topX, dLL)
      AgeDiff(i,j) = ADX
      if (ADX /= 999)  AgeDiff(j,i) = -ADX
      if (topX==6 .or. topX==7) cycle   ! conditionally unrelated
      if (COUNT(Parent == 0) > 0.95*nInd .and. topX>2)  cycle  ! else will exceed nAmbMax
    
    else if (ParSib == 2) then 
      maybe = 0
      LRR = 999D0
      LRR(1) = MaxLL((/LLR_O(i,j), LLR_O(j,i)/))
      topX = 0
      do k=1,2 
        if (Parent(i,k)/=0 .and. Parent(i,k)==Parent(j,k)) cycle
        if (Parent(i,k)>0) then
            if (ANY(Parent(Parent(i,k), :)==j)) cycle
        else if (Parent(i,k)<0) then
            if (ANY(GpID(:, -Parent(i,k), k)==j)) cycle
        endif
        if (Parent(j,k)>0) then
            if (ANY(Parent(Parent(j,k), :)==i)) cycle
        else if (Parent(j,k)<0) then
            if (ANY(GpID(:, -Parent(j,k), k)==i)) cycle
        endif                                                                                                              
        call PairQFS(i, j, LRR(2)) 
        call PairQHS(i, j, LRR(3)) 
        maybe = 0
        do x=1,3
          if (LRR(x) > 2*TA .and. LRR(x) < 999)  maybe=1  
        enddo
        if (maybe==0)  cycle
        if (AgeDiff(i,j)>=0) then
          call CheckPair(i, j, k, 7, LL, LLX)  
        else
          call CheckPair(j, i, k, 7, LL, LLX)
        endif
        call BestRel2(LL, topX, dLL)
        if (COUNT(Parent == 0) > 0.95*nInd .and. topX>2) then
          cycle  ! else will exceed nAmbMax)
        else if (topX==6 .or. topX==7 .or. topX==8) then  ! .or. dLL(2)<TA
          maybe = 0
          cycle
        else
          exit
        endif
      enddo
      if (maybe==0) cycle
    endif
    
    nAmb = nAmb + 1
    AmbigID(nAmb) = i
    AmbigID(nAmbMax + nAmb) = j
    AmbigOH(nAmb) = OppHomM(i,j)
    if (ParSib==1) then
      AmbigLR(nAmb) = MIN(LLR_O(i,j), LLR_O(j,i))
      AmbigRel(nAmb) = 1
    else if (ParSib==2) then
      AmbigLR(nAmb) = MAXVAL(LRR, MASK=LRR<222)
      AmbigRel(nAmb) = MAXLOC(LRR, MASK=LRR<222, DIM=1)
    endif
    AmbigRel(nAmbMax + nAmb) = TopX
    AmbigLR(nAmbMax + nAmb) = dLL(1)
    if (nAmb==nAmbMax) then
      if(quiet<1) then
!        call rwarn("reached max for maybe-rel")  ! doesn't return any output
        call intpr(" ",-1, 0, 0)
        if (ParSib == 1) then
          call intpr("WARNING - reached max for maybe-par, truncated!",-1, 0, 0)
        else
          call intpr("WARNING - reached max for maybe-rel, truncated!",-1, 0, 0)
        endif
        call intpr(" ",-1, 0, 0)
      endif
      exit
    endif
  enddo
  if (nAmb==nAmbMax)  exit                                                                    
enddo

if (nAmb>1) then  
  if (COUNT(AmbigRel((nAmbMax+1) : (2*nAmbMax)) == 1) > 1) then
    if(quiet<1) then
      call intpr("Checking for Parent-Parent-Offspring trios ... ",-1, 0, 0)
    endif
    call triads(nAmbMax, AmbigID, AmbigRel((nAmbMax+1) : (2*nAmbMax)),&
   ntrio, trioIDs, trioLR, trioOH)                
  endif
endif

call DeAllocAll

end subroutine findambig

! #####################################################################
                                                                    
subroutine Erstop(message)
use Global
implicit none

 character(len=*), intent(IN) :: message
! Error = 1
 call DeAllocAll
 call rexit("  ERROR! ***"//message//"*** Results are not reliable!")

end subroutine Erstop

! ####################################################################

subroutine parents(TotLL)
use qsort_c_module 
use Global
implicit none

double precision, intent(INOUT) :: TotLL(42)
integer :: i, j, k, Round, isP(2), PriorPed(nInd, 2)
integer, allocatable, dimension(:) :: BYRank
double precision, allocatable, dimension(:) :: SortBY

if (Error/=0) return
 call rchkusr()     
 
allocate(OppHomM(nInd, nInd))
OppHomM = 999 
allocate(LLR_O(nInd, nInd))
LLR_O = 999D0
PriorPed = Parent
Parent = 0
!============================

 call UpdateAllProbs(2)
if(quiet<1)  call intpr ("Parentage ... ", -1, 0, 0)
if(quiet<1)  call dblepr("Initial total LL : ", -1, SUM(Lind), 1) 
 call CalcOppHom   ! also checks no. SNPs typed in both
 
do i=1, nInd-1
  do j=i+1,nInd 
    if (OppHomM(i,j) > maxOppHom) cycle
      call CalcPO(i, j, LLR_O(i,j))  ! LLR PO/U
      call CalcPO(j, i, LLR_O(j,i))
  enddo
enddo

! get birthyear ranking (increasing)
allocate(SortBY(nInd))
allocate(BYRank(nInd))
SortBY = REAL(BY, 8)
WHERE (SortBY < 0) SortBY = HUGE(0.0D0) 
BYRank = (/ (i, i=1, nInd, 1) /)
if(ANY(BY>=0)) call QsortC(SortBy, BYRank)

TotLL = 0
 call UpdateAllProbs(2)
TotLL(1) = SUM(Lind)
do Round=1,41
  if (Error/=0) return
  call rchkusr()
  call Parentage(BYrank, PriorPed)   
  call UpdateAllProbs(3)
  do i=1,nInd
    if (Sex(i)==3) then
      isP = 0
      do k=1,2
        do j=1,nInd
          if (Parent(j,k) == i) then
            isP(k) = isP(k) + 1
          endif
        enddo
      enddo
      if (isP(1)>0 .and. isP(2)>0) then
        call rwarn("individual assigned as both dam & sire")
      else
        do k=1,2
          if (isP(k)>1) then
            Sex(i) = k
          endif
        enddo
      endif
    endif
  enddo

  TotLL(Round + 1) = SUM(Lind)
  if (TotLL(Round + 1) - TotLL(Round) < ABS(TF)) exit ! convergence
  if (Round==41) then
    call Erstop("parentage not converging")
  endif
enddo

if(quiet<1) call dblepr("Post-parentage total LL : ", -1, SUM(Lind), 1) 

deallocate(BYRank)
deallocate(SortBY)

end subroutine parents

! ####################################################################

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! ####################################################################

subroutine sibships(Nrounds, UseAge, TotLL)
use Global
implicit none

integer, intent(INOUT) :: Nrounds, UseAge   ! IN
double precision, intent(INOUT) :: TotLL(42)                 
double precision :: CurLL(8), PrevLL(8)  !, LLTMP(4)
integer :: Round, RX, RoundEA, u=3
 character(len=2) :: RoundChars(42)
 if (Error/=0) return
 call rchkusr()

 RX = 1  ! no. of initial rounds, pairs-cluster-merge only                                                                         
RoundEA = 0  ! round number with extra ageprior (for GGpairs)
XP = 5  ! max no. candidate sib pairs

allocate(PairID(XP*nInd, 2)) 
allocate(PairDLLR(XP*nInd))
allocate(PairType(XP*nInd))  ! mat (1), pat (2) or unknown (3)
!=======================

! K.I.S.S.
! intpr() for iteration & dblepr() for LL --> 4 lines/iteration in R console
! w-rite(RoundC, '(i2)') Round --> use of Fortran I/O will when compiled on Windows 
! interfere with C I/O: when the Fortran I/O support code is initialized (typically 
! when the package is loaded) the C stdout and stderr are switched to LF line endings
RoundChars = (/ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', &
  '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', &
  '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', &
  '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42'/)

 call UpdateAllProbs(u)
if(quiet<1) then
  call dblepr("Sibships - Initial Total LL : ", -1, SUM(Lind), 1)
endif

curLL = 999
TotLL = 999
           
if (UseAge==0) then   !  .or. nAgeClasses==1
  AgePhase = 0
else
  AgePhase = 1
endif
if (Nrounds >= 42) then
  call Erstop("too many sib rounds")
endif

if (.not. Hermaphrodites) then
  call MoreParent(.TRUE.)  ! double check parents, using updated ageprior
endif    

do Round=1, Nrounds
  if (Error/=0) return
  call rchkusr()
  call UpdateAllProbs(u)
  TotLL(Round) = SUM(Lind)
  CurLL(1) = SUM(Lind)   
  if(quiet==-1)  call intpr("--- Round "//RoundChars(Round)//" start ---", -1, 0,0)  
  if(quiet==-1)  call intpr("Find pairs ...", -1, 0,0)  
  call FindPairs   ! do not use age priors yet 
!  call intpr("npairs: ", -1, nPairs, 2) 
  if (Error/=0) return
  call rchkusr()
  if(quiet==-1)  call intpr("Clustering ...", -1, 0,0)
  call Clustering       
  call UpdateAllProbs(u)                                                                                                                                                            
  CurLL(2) = SUM(Lind)
  if (Round > RX+1 .and. nAgeClasses > 1) then 
  if(quiet==-1)  call intpr("Grandparents ...", -1, 0,0)  
    call GGpairs(RoundEA)  ! lag one round with extra age prior  
    call UpdateAllProbs(u)
    CurLL(3) = SUM(Lind)                                                                                          
  endif 
  if(quiet==-1)  call intpr("Merge clusters ...", -1, 0,0)  
  call Merging
  call UpdateAllProbs(u)        
  CurLL(4) = SUM(Lind)
  if (Error/=0) return
  call rchkusr() 
  if (Round > RX .or. Round==Nrounds) then
    if (nAgeClasses > 1) then   
      if(quiet==-1)  call intpr("Sibship parent replacement...", -1, 0,0)    
      call SibParent   ! replace dummy parents by indivs
      call UpdateAllProbs(u)
      CurLL(5) = SUM(Lind)                                                
    endif  
    if(quiet==-1)  call intpr("Grow clusters ...", -1, 0,0)    
    call GrowClusters                  
    call UpdateAllProbs(u)
    CurLL(6) = SUM(Lind)                                                                                          
  endif   

  if (nAgeClasses > 1) then                             
    if (.not. Hermaphrodites) then
      if(quiet==-1)  call intpr("Parents ...", -1, 0,0)
      call MoreParent(.FALSE.)  !  assign additional parents to singletons   
    endif  
    call UpdateAllProbs(u)
    CurLL(7) = SUM(Lind)                                                                                                   
    if ((Round > RX .or. Round==Nrounds) .and. nAgeClasses > 1) then
      if(quiet==-1)  call intpr("Sibship grandparents ...", -1, 0,0)
      call SibGrandparents
      call UpdateAllProbs(u)
      CurLL(8) = SUM(Lind)
    endif
  endif
  
  if(quiet<1) then
    call dblepr("Round "//RoundChars(Round)//" end, Total LogLik: ", -1, SUM(Lind), 1)
    call intpr("No. dams, sires (incl. dummies): ", -1, count(Parent/=0, DIM=1) ,2)
    if(quiet==-1)  call intpr("No. dummies: ", -1, nC ,2)
!    LLTMP = (/ MaxLL(curLL), MaxLL(PrevLL(1:7)), curLL(8), MINVAL(curLL(1:7)) /)
!    call dblepr("LL: ", -1, LLTMP , 4)
    if(quiet==-1)  call intpr("---------------", -1, 0,0)
  endif
  
  if (Round == 1 .and. curLL(8) - MaxLL(curLL(1:7)) < ABS(TF)) then
    call Erstop("LL not decreased in round 1, terminating.")
  else if (Round == nRounds) then                                                      
    exit
  else if (Round>1 .and. (MaxLL(curLL) - MaxLL(PrevLL(1:7)) < 2*ABS(TF) .or. &
   (curLL(8) - MINVAL(curLL(1:7))) < 2*ABS(TF))) then  
    if (UseAge==2 .and. AgePhase==1 .and. nAgeClasses>1) then
      AgePhase = 2
    else
      exit
    endif
  endif
  if (AgePhase==2)  RoundEA = RoundEA +1
  PrevLL = CurLL
  curLL = 999
enddo  
 
call UpdateAllProbs(u)
if (Round<nRounds) then
    TotLL(Round+1) = SUM(Lind)      
endif
  
end subroutine Sibships

! #####################################################################

! @@@@   SUBROUTINES   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! #####################################################################

subroutine CalcOppHom  ! nInd x nInd matrix with no. opp. hom. loci
use Global
implicit none

integer :: i, j, Lboth

OppHomM = -999
do i=1, nInd-1
  do j=i+1,nInd
    call CalcOH(i, j, OppHomM(i,j))
    OppHomM(j,i) = OppHomM(i,j)
    if (OppHomM(i,j) <= maxOppHom) then
      Lboth = COUNT(Genos(:,i)/=-9 .and. Genos(:,j)/=-9)    
      if (Lboth < nSnp/2.0) then   ! >1/2th of markers missing
        OppHomM(i,j) = -Lboth 
        OppHomM(j,i) = -Lboth
      endif
    endif
  enddo
enddo

end subroutine CalcOppHom

! #####################################################################

subroutine CalcOH(A,B,OH)
use Global
implicit none

integer, intent(IN) :: A, B
integer, intent(OUT) :: OH
integer :: l

OH = 0
do l=1,nSnp
  if ((Genos(l,A)==1).and.(Genos(l,B)==3)) then
    OH = OH+1
    if (OH > maxOppHom) exit
  endif                       
  if ((Genos(l,A)==3).and.(Genos(l,B)==1)) then
    OH = OH+1
    if (OH > maxOppHom) exit
  endif                       
enddo

end subroutine CalcOH

! #####################################################################

 subroutine CalcTrioErr(A,ME)   ! Mendelian errors in offspring-parent-parent trios
use Global
implicit none

integer, intent(IN) :: A
integer, intent(OUT) :: ME
integer :: l, k, Ecnt(3,3,3)   ! offspring - dam - sire

Ecnt(:,1,1) = (/ 0, 1, 2 /)
Ecnt(:,1,2) = (/ 0, 0, 1 /)
Ecnt(:,1,3) = (/ 1, 0, 1 /)

Ecnt(:,2,1) = (/ 0, 0, 1 /)
Ecnt(:,2,2) = (/ 0, 0, 0 /)
Ecnt(:,2,3) = (/ 1, 0, 0 /)

Ecnt(:,3,1) = (/ 1, 0, 1 /)
Ecnt(:,3,2) = (/ 1, 0, 0 /)
Ecnt(:,3,3) = (/ 2, 1, 0 /)

ME = 0
do l=1,nSnp
  if (Genos(l,A)==-9 .or. ALL(Genos(l, Parent(A,:))==-9)) then
    cycle
  else if (ANY(Genos(l, Parent(A,:))==-9)) then
    do k=1,2
      if (Genos(l, Parent(A,k))==-9) then
        if (((Genos(l,A)==1).and.(Genos(l,Parent(A,3-k))==3)) .or. &
         ((Genos(l,A)==3).and.(Genos(l,Parent(A,3-k))==1))) then
          ME = ME +1
          cycle
        endif
      endif
    enddo
  else
    ME = ME + Ecnt(Genos(l,A), Genos(l, Parent(A,1)), Genos(l, Parent(A,2)))
  endif
enddo

end subroutine CalcTrioErr

! #####################################################################

subroutine AtoVi(A, d1, d2, x, V)
use Global
implicit none

integer, intent(IN) :: d1, d2, x(2)
integer, intent(IN) :: A(d1,d2,2)   ! IN: integer
integer, intent(OUT) :: V(d1*d2*2)
integer :: i, j, k

V = 0
do j=1,d1
  do k=1,2  ! works for d(3)=2
    do i=1, x(k)
      V((j-1)*2*d2 + (k-1)*x(1) + i) = A(j, i, k)
    enddo
  enddo
enddo

end subroutine AtoVi

! #####################################################################

subroutine AtoVd(A, d1, d2, x, V)
use Global
implicit none

integer, intent(IN) ::  d1, d2, x(2)
double precision, intent(IN) :: A(d1,d2,2)  ! IN: double
double precision, intent(OUT) :: V(d1*d2*2)
integer :: i, j, k

V = 0
do j=1,d1
  do k=1,2  ! works for d(3)=2
    do i=1, x(k)
      V((j-1)*2*d2 + (k-1)*x(1) + i) = A(j, i, k)
    enddo
  enddo
enddo

end subroutine AtoVd

! #####################################################################

subroutine CalcPO(A,B, LLR)  
! LLR of A as offspring from B, vs A as random sample from pop
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LLR
integer :: l
double precision :: PrL(nSnp)

LLR = 999D0
PrL = 0D0
do l=1,nSnp
  if (Genos(l,A)/=-9 .and. Genos(l,B)/=-9) then
    PrL(l) = LOG10(OKOP(Genos(l,A), Genos(l,B), l))
  else if (Genos(l,A)/=-9) then
    PrL(l) = LOG10(OHWE(Genos(l,A),l))
  endif
enddo
LLR = SUM(PrL) - Lind(A)

end subroutine CalcPO

! ######################################################################

subroutine CalcPO2(A,B,C, LLR)   ! LL of A, given B & C as parent  
use Global
implicit none

integer, intent(IN) :: A, B, C
double precision, intent(OUT) :: LLR
integer :: l, x, y
double precision :: PrL(nSnp), tmp(3,3)

LLR = 999D0
PrL = 0D0
do l=1,nSnp
  if (Genos(l,A)==-9) cycle
  do x=1,3
    do y=1,3
      if (B>0 .and. C>0) then
        tmp(x,y) = OKA2P(Genos(l,A),x,y) * LindG(x,l,B) * LindG(y,l,C)
      else if (B>0 .and. C==0) then
        tmp(x,y) = OKA2P(Genos(l,A),x,y) * LindG(x,l,B) * AHWE(y,l)
      else if (B==0 .and. C>0) then
        tmp(x,y) = OKA2P(Genos(l,A),x,y) * AHWE(x,l) * LindG(y,l,C)
      else if (B==0 .and. C==0) then
        tmp(x,y) = OKA2P(Genos(l,A),x,y) * AHWE(x,l) * AHWE(y,l)
      else
        call rexit("invalid call to CalcPO2")
            
      endif
    enddo
  enddo
  PrL(l) = LOG10(sum(tmp))
enddo
LLR = SUM(PrL)

end subroutine CalcPO2

! ######################################################################

subroutine Parentage(BYrank, PriorPed)
use Global
implicit none

integer, intent(IN) :: BYrank(nInd), PriorPed(nInd, 2) 
integer :: i, j, x, y, k, CandPar(51, 2), nCP(2), u, v, AncJ(2,mxA),&
 nof, offspr(maxSibSize), sxOff(maxSibSize), tmp(2), curPar(2), tmpPar(2)
logical :: skip
double precision :: LLX(51,51), LLY(51,51), tmpLind

do x=1, nInd
if (Error/=0) return
  if (MOD(x,250)==0) call rchkusr()
!          call intpr ( " ",2, x, 1)
!   endif
  i = BYRank(x)
  nCP = 0
  CandPar = 0
   curPar = Parent(i,:)                   
  do y=1,nInd 
    j = BYRank(y)
    if (i==j) cycle
    if (ANY(Parent(j,:)==i)) cycle
    if (ANY(Parent(i,:)==j)) then
      do k=1,2
        if(Parent(i,k)==j) then
          if (nCP(k)==50) cycle
          nCP(k) = nCP(k) + 1
          CandPar(nCP(k), k) = j
       endif
      enddo  
    else
      if (OppHomM(i,j) > maxOppHom .or. OppHomM(i,j)<0) cycle 
      if (LLR_O(i,j) < TF) cycle  
      call GetAncest(j,1,AncJ)
      if (ANY(AncJ == i)) cycle
      if (AgeDiff(i,j) == 999D0) then
        skip = .FALSE.
        if (BY(i)/=-999 .and. BY(j)==-999) then
          do k=1,2
            do u=2, mxA
              if (AncJ(k,u)>0) then
                if (AgeDiff(i, AncJ(k,u))<= 0) skip = .TRUE.
              endif
            enddo
          enddo
        endif
        if (skip) cycle
        if (BY(i)==-999 .and. BY(j)/=-999 .and. Sex(i)/=3) then
          if (ANY(Parent(:, Sex(i))==i)) then
            call GetOff(i,sex(i),.FALSE., nof, offspr, sxOff)
            if (ANY(AgeDiff(offspr(1:nof), j) <= 0))  cycle
          endif
        endif
      else 
        if (AgeDiff(i,j) <= 0)  cycle
      endif    
      if (.not. (Sex(j)==3 .and. ALL(Parent(i,:)==0))) then
        call CalcPOZ(i,j, .FALSE.)  ! assigns parent as side effect
      endif
      do k=1,2
        if (Sex(j) /= 3 .and. Sex(j)/= k) cycle
        if (nCP(k)==50) cycle
        nCP(k) = nCP(k) + 1
        CandPar(nCP(k), k) = j
      enddo
    endif
  enddo
                
  if (ALL(nCP == 1) .and. ALL(CandPar(1,:) == curPar))  cycle                                                                                                                                    
  if (ALL(nCP>0)) then   ! test combo's
    if (ALL(CandPar >= 0) .and. BY(i)/=-999 .and. &
      ALL(BY(CandPar(1:nCP(1), 1))/=-999) .and. &
      ALL(BY(CandPar(1:nCP(2), 2))/=-999) .and. &
      ALL(Sex(CandPar(1:nCP(1), 1))==1) .and. &
      ALL(Sex(CandPar(1:nCP(2), 2))==2)) then
      LLX = -999
      LLY = -999                    
      do u = 1, nCP(1)+1    
        do v = 1, nCP(2)+1
          call CalcPO2(i,CandPar(u, 1), CandPar(v, 2), LLX(u,v))
        enddo
      enddo      
      if (MAXVAL(LLX - LLX(nCP(1)+1, nCP(2)+1))<2*TA) cycle 
      tmp = MAXLOC(LLX)
      u = tmp(1)
      v = tmp(2)
      LLY = LLX - MAXVAL(LLX)  ! negative values
      if (COUNT(ABS(LLY) < 0.01) < 2) then
        LLY(u,v) = -888
      endif
      LLX = LLX - LLX(nCP(1)+1, nCP(2)+1)
      
      if (MAXVAL(LLY) < -2*TA) then  ! vs next most likely parents  .and. ALL(nCP>1)
        if (CandPar(u, 1)==0 .and. CandPar(v, 2)==0) then
          call NewPar(i,1,0)
          call NewPar(i,2,0)
        else if (CandPar(u, 1)>0 .and. CandPar(v, 2)>0) then
          if (.not. (Parent(i,1)==CandPar(u,1) .and. Parent(i,2)==CandPar(v,2))) then    ! double check
            if (Parent(i,1)==CandPar(u,1)) then
              call CalcPOZ(i,CandPar(v, 2), .FALSE.)
            else if (Parent(i,2)==CandPar(v,2)) then
              call CalcPOZ(i,CandPar(u, 1), .FALSE.)
            else
              call NewPar(i,1, CandPar(u,1))
              call CalcPOZ(i,CandPar(v, 2), .FALSE.)
            endif
          endif
        else if (CandPar(u, 1)>0 .and. Parent(i,1)/=CandPar(u,1)) then
          call CalcPOZ(i,CandPar(u, 1), .FALSE.) 
        else if (CandPar(v, 2)>0 .and. Parent(i,2)/=CandPar(v,2)) then
          call CalcPOZ(i,CandPar(v, 2), .FALSE.)   
        endif
      else if (MAXVAL(LLY) > -TA/2) then ! use pedigree prior
        call NewPar(i,1,0)
        call NewPar(i,2,0)
        u = nCP(1)+1
        v = nCP(2)+1
        if (ANY(CandPar(1:nCP(1), 1) == PriorPed(i,1))) then
          u = MINLOC(abs(CandPar(:, 1) - PriorPed(i,1)), dim=1)
          v = MAXLOC(LLX(u,:), dim=1) 
        else if (ANY(CandPar(1:nCP(2), 2) == PriorPed(i,2))) then
          v = MINLOC(abs(CandPar(:, 2) - PriorPed(i,2)), dim=1) 
          u = MAXLOC(LLX(:,v), dim=1)
        endif
        if (LLX(u,v) > 2*TA .and. (LLY(u,v)==-888 & 
         .or. LLY(u,v)>-TA/2)) then
          call NewPar(i,1, CandPar(u,1))
          call NewPar(i,2, CandPar(v,2))
        endif
      endif
      call CalcLind(i)
    else if (ALL(Parent(i,:)==0) .or. ALL(nCP>1)) then
      do u = 1, nCP(1)    
        do v = 1, nCP(2)
          if (CandPar(u, 1)==CandPar(v, 2))  cycle  !when unk sex.
          if (Sex(CandPar(u,1))==3 .and. Sex(CandPar(v,2))==3) cycle
          tmpPar = Parent(i,:)
          tmpLind = Lind(i)
          if (Sex(CandPar(u,1))/=3 .and. .not. ANY(Parent(i,:)==CandPar(v,2))) then 
            if (Parent(i,2)/=CandPar(u,1))  call NewPar(i,1, CandPar(u,1))
            if (sex(CandPar(v,2))==3)       call NewPar(i, 2, 0) 
            call calcPOZ(i, CandPar(v,2), .FALSE.) 
            if (Parent(i,2)==0) then    ! undo temp (un)assignment 
              call NewPar(i, 1, tmpPar(1))
              call NewPar(i, 2, tmpPar(2))
            endif            
          else if (.not. ANY(Parent(i,:)==CandPar(u,1))) then
            if(Parent(i,1)/=CandPar(v,2))   call NewPar(i,2, CandPar(v,2))
            if (sex(CandPar(u,1))==3)       call NewPar(i, 1, 0) 
            call calcPOZ(i, CandPar(u,1), .FALSE.) 
            if (Parent(i,1)==0) then
              call NewPar(i, 1, tmpPar(1))
              call NewPar(i, 2, tmpPar(2))
            endif
          endif
          if (Lind(i) < tmpLind) then   ! undo
            call NewPar(i, 1, tmpPar(1))
            call NewPar(i, 2, tmpPar(2))
          endif
        enddo
      enddo
  
      if (Parent(i,1)/=0 .and. Parent(i,2)/=0) then  ! double check
        if (Sex(Parent(i,1))==3 .and. Sex(Parent(i,2))==3) then    ! can't know which way around
          call NewPar(i, 1, 0)
          call NewPar(i, 2, 0) 
        else
          tmpPar(1) = Parent(i,1)
          call NewPar(i, 1, 0)
          call calcPOZ(i, tmpPar(1), .FALSE.) 
          if (Parent(i,1)==0 .or. Parent(i,2)==0) then    ! no combination found
            call NewPar(i, 1, 0)
            call NewPar(i, 2, 0) 
          else 
            do k=1,2
              if (Sex(Parent(i,k))==3)  Sex(Parent(i,k)) = k
            enddo
          endif
        endif
      endif
    endif

    do k=1,2
      if (Parent(i,k)==0) cycle
      if (Sex(Parent(i,k))==3) then
        if (Parent(i,3-k)==0) then
          call NewPar(i, 1, 0)
          call NewPar(i, 2, 0) 
        else if (Sex(Parent(i,3-k))==3) then
          call NewPar(i, 1, 0)
          call NewPar(i, 2, 0) 
        else
          Sex(Parent(i,k)) = k
        endif
      endif
    enddo    
  endif  ! ALL(nCP>0)
  call CalcLind(i)
enddo

end subroutine Parentage

! #####################################################################

subroutine CalcPOZ(A, B, UseAge)  ! replace a current parent of A by B? 
use Global
implicit none 

integer, intent(IN) :: A, B
logical, intent(IN) :: UseAge
integer :: m, CurPar(2), TopX, k, CY(3), kY(3), y, cpb(2)
double precision :: LLA(2,7,7), TopLL, LLcp(3), LLBA(7), TopBA, dLL, &
LLtmp(4), ALR(3), ALRtmp(2), LLAp
logical :: OK, AgeOK                    

 CurPar = Parent(A,:)
 cpb = Parent(B,:)
 call CalcLind(A)                 

if (Sex(A)==3 .and. Sex(B)==3 .and. ALL(Parent(A,:)==0) .and. ALL(Parent(B,:)==0)) then
  return  ! can't tell if A or B is dam or sire
endif
if (Sex(B)/=3) then
  k = Sex(B) 
else if (ALL(Parent(A,:)==0) .or. ALL(Parent(A,:)/=0)) then
  return  ! can't tell if B is father or mother
else 
  do m=1,2
    if (Parent(A,m)==0) then
      k = m  ! try B as opposite-sex parent
    endif
  enddo   
endif

LLtmp = 999D0
if (Sex(B)==3 .and. CurPar(3-k)<0) then
  call AddParent(B, -CurPar(3-k),3-k, LLtmp(1))      
  call CalcU(B,3-k, CurPar(3-k),3-k, LLtmp(2))
  if ((LLtmp(1) - LLtmp(2)) > -TA) then
    return  ! B likely to replace CurPar(3-k) instead
  endif
endif
ALRtmp = 999D0
AgeOK = .FALSE.
if (AgeDiff(A,B)==999) then
!   Parent(A,:) = 0  ! do not rely on (temporarily) assigned parent for age info
  call CalcAgeLR(A, 0, B, 0, k, 1, .FALSE., ALRtmp(1))
!  Parent(A,:) = CurPar
  if (ALRtmp(1) < 3*TF .or. ALRtmp(1)==777) return
  if (abs(ALRtmp(1))>0.01) then
    Parent(A,:) = 0
    call CalcAgeLR(B, 0, A, 0, k, 1, .FALSE., ALRtmp(2))
    Parent(A,:) = CurPar                 
    if (ALRtmp(2)/=777 .and. ALRtmp(2) - ALRtmp(1) > ABS(TF)) then
      return
    else if (ALRtmp(2) == 777 .or. ALRtmp(1) - ALRtmp(2) > 2*ABS(TF)) then
      AgeOK = .TRUE.
    endif
  else if (abs(ALRtmp(1)) <= 0.01) then  ! no age info
    do m=1,2
      if ((Sex(A)==m .or. Sex(A)==3) .and. Parent(B,m) < 0) then
        call AddParent(A, -Parent(B,m),m, LLtmp(1)) 
        call CalcU(A,m, Parent(B,m),m, LLtmp(2))
        if (LLtmp(1)<0 .and. (LLtmp(1) - LLtmp(2)) > -TA) then 
          return    ! A likely to replace Parent(B,m) instead
        endif
      endif
    enddo
  endif
endif

if (Parent(A,3-k) < 0) call CalcCLL(-Parent(A,3-k), 3-k)
 call CalcLind(A)
 call CalcLind(B)
TopX = 0
LLA = 999D0
LLBA = 999D0
TopLL = 999D0
TopBA = 999D0
if (ALL(Parent(A,:)==0)) then
  call CheckRel(A, k, B, k, UseAge, 1, LLA(1,:,7))  
  call BestRel(LLA(1,:,7), 1, TopX, dLL)
  TopLL = MaxLL(LLA(1,:,7))                                                       
  if (AgeDiff(A,B)==999 .and. .not. AgeOK) then
    do m=1,2
      if (Sex(A)/=3 .and. Sex(A)/=m) cycle
      if (Parent(B,m)>0) cycle
      if (Parent(B,m)<0)  call NewPar(B, m, 0)
      call CheckRel(B, m, A, m, UseAge, 1, LLBA)
      if (Parent(B,3-m) < 0) then  ! include changes in CLL
        call CalcU(B,m, Parent(B,3-m),3-m, LLtmp(1))
        call NewPar(B, m, A)
        call CalcU(B,m, Parent(B,3-m),3-m, LLtmp(2))
        LLBA(1) = LLBA(1) + (LLtmp(2) - LLtmp(1))
        call NewPar(B, m, cpb(m))            
      endif
      do y=1,2
        call NewPar(B, y, cpb(y))
      enddo
      TopBA = MaxLL(LLBA)  
      if (ABS(TopBA - TopLL) < TA .or. (TopBA > TopLL)) return
    enddo
  endif
  if (TopX==1) then
    call NewPar(A, k, B)
  endif

else 
  LLcp = 0  ! need LL over all (2-4) indiv involved
  TopLL = 999D0
  call CalcU(CurPar(1), 1, CurPar(2), 2, LLcp(3))
  do m=1,2
    call CalcU(CurPar(3-m), 3-m, B, k, LLcp(m))
    if (curPar(3-m) < 0) then
      LLCP(m) = LLCP(m) - Lind(A)
      LLCP(3) = LLCP(3) - Lind(A)  ! never called when both par <0
    endif
  enddo

  ! age based prior prob.
  ALR = 0
  if (UseAge) then
    CY = (/ CurPar(1), CurPar(2), B /)
    kY = (/ 1, 2, k /)
    do y=1,3
      if (y<3)  call NewPar(A, kY(y), 0) 
      call CalcAgeLR(A,Sex(A), CY(y),kY(y), 0,1, .FALSE., ALR(y))
      if (y<3)  call NewPar(A, kY(y), CY(y))
    enddo
  endif
     
  do m=1,2
    if (CurPar(m)==0) cycle  ! LLA(m,:,:) empty if only 1 CurPar
    call CheckRel(A, 0, B, k, UseAge, 1, LLA(m,:,1))
    call NewPar(A, m, 0)                                                       
    call CheckRel(A, 0, B, k, UseAge, 1, LLA(m,:,7))    
    call CheckRel(A, 0, CurPar(m), m, UseAge, 1, LLA(m,7,:))     
    call NewPar(A, k, B)
    LLAp = LLA(m,1,1)
    call CheckRel(A, 0, CurPar(m), m, UseAge, 1, LLA(m,1,:))  
    if (CurPar(m)<0) then 
      LLA(m,7,2) = 333   ! FS does not count here.
      call ReOrderAdd(LLA(m,7,:))
      call ReOrderAdd(LLA(m,1,:))         
    endif
    LLA(m,1,1) = LLAp                 
    if (Parent(B, m)==CurPar(m) .and.  CurPar(m)/= 0) then  
      LLA(m,2:3,7) = 888  ! HS implies curPar = Par
      if (Parent(B, 3-m)==CurPar(3-m) .and.  CurPar(3-m)/= 0) then
        LLA(m,2,7) = 888
      endif
    endif
    
    call NewPar(A, 1, CurPar(1))   ! restore
    call NewPar(A, 2, CurPar(2))
    WHERE (LLA(m,1:6,1)<0) LLA(m,1:6,1) = LLA(m,1:6,1) + LLcp(3) + &
      ALR(1) + ALR(2)
     WHERE (LLA(m,2:6,7)<0) LLA(m,2:6,7)=LLA(m,2:6,7)+LLcp(3)+ALR(3-m)
    WHERE (LLA(m,1,2:7)<0) LLA(m,1,2:7) = LLA(m,1,2:7) + LLcp(m) + ALR(3)
    if (m==k) WHERE (LLA(m,1,2:7)<0) LLA(m,1,2:7) = LLA(m,1,2:7) + ALR(3-m)
    WHERE (LLA(m,7,:)<0) LLA(m,7,:) = LLA(m,7,:) + LLcp(m) +ALR(3-m) 

     if (curPar(m)>0) then
      ! check if FS trio (ignores all parents)
      ALRtmp = 999D0
      call CalcAgeLR(A,Sex(A), B,k, 0,2, .FALSE., ALRtmp(1))
      call CalcAgeLR(A,Sex(A),  CurPar(m),m, 0,2, .FALSE., ALRtmp(2))
      if (all(ALRtmp/=777)) then
        call trioFS(A, B, CurPar(m), LLA(m,2,2)) 
      endif
      !check of both B & curPar(m) are FS of A's true parent
      ALRtmp = 999D0
      call CalcAgeLR(A,Sex(A), B,k, 0,5, .FALSE., ALRtmp(1))
      call CalcAgeLR(A,Sex(A),  CurPar(m),m, 0,5, .FALSE., ALRtmp(2))
      if (all(ALRtmp/=777)) then
        call trioFA(A, B, CurPar(m), LLA(m,5,5))   
      endif
    endif
  enddo  
  
  TopLL = MaxLL(RESHAPE(LLA(:,:,:), (/2*7*7/)))!MAXLOC doesn't do ties

  if (AgeDiff(A,B)==999 .and. .not. AgeOK) then
    LLBA = 999D0
    do m=1,2
      if (Sex(A)/=3 .and. Sex(A)/=m) cycle
      if (Parent(B,m)>0)  cycle
      if (Parent(B,m)<0)  call RemoveSib(B, -Parent(B,m), m)
      call CheckRel(B, 0, A, m, UseAge, 1, LLBA)
      if (Parent(B,3-m) < 0 .and. LLBA(1)<0) then  
        call CalcU(Parent(B,3-m),3-m,0,0, LLtmp(1))
        call NewPar(B, m, A)                 
        call CalcU(Parent(B,3-m),3-m,0,0, LLtmp(2))
        LLBA(1) = LLBA(1) + (LLtmp(2) - LLtmp(1))
        call NewPar(B, m, cpb(m))   ! restore                            
      endif
      do y=1,2
        if (Parent(A,y) < 0 .and. LLBA(1)<0 .and. .not. &
          (y==3-m .and. Parent(A,y)==Parent(B,3-m))) then 
          call CalcU(Parent(A,y),y,B,m, LLtmp(3))
          call NewPar(B, m, A)                      
          call CalcU(Parent(A,y),y,B,m, LLtmp(4))
          LLBA(1) = LLBA(1) + (LLtmp(4) - LLtmp(3))  
          call NewPar(B, m, cpb(m)) 
        endif
      enddo
      WHERE (LLBA<0) LLBA = LLBA + LLcp(3)
      TopBA = MaxLL(LLBA)  
      if (ABS(TopBA - TopLL) < TA .or. (TopBA > TopLL))  return
    enddo
    do y=1,2
      if (cpb(y)<0) call DoAdd(B, -cpb(y), y)
    enddo
    do m=1,2
      call NewPar(B, m, cpb(m))
    enddo
  endif
  
  do m=1,2
    call NewPar(A, m, 0)
  enddo    
  if (LLA(3-k,1,1)==TopLL .or. ANY(LLA(k,1,:) == TopLL)) then   ! TODO: threshold?
    call NewPar(A, k, B)
    call NewPar(A, 3-k, curPar(3-k))
!  else if (ANY(LLA(:,2,2)==TopLL)) then
    ! do nothing         
  else if ((TopLL - MaxLL(RESHAPE(LLA(:,2:7,1), (/2*6/)))) <0.01) then  !  .and. ALL(curPar/=0)
    call NewPar(A, k, curPar(k))
    call NewPar(A, 3-k, curPar(3-k))
  else if (ANY(LLA(:,1,:)==TopLL) .and. ALL(TopLL - LLA(k, 2:7, 2:7) > TF)) then  ! only B
    call BestRel(LLA(k,:,1),1,topX, TopBA)
    if ((MaxLL(LLA(k,1,:)) - TopBA) > TA) then
      call NewPar(A, k, curPar(k))
    endif
  else if (ANY(LLA(k, 2:7, 2:7) == TopLL) .and. ALL(TopLL - LLA(:,1,:) > TF)) then  ! keep CurPar(3-k). TF OR TA? CHECK
    call NewPar(A, 3-k, curPar(3-k))                                  
  else if (ANY(LLA(3-k, 2:7, 2:7) == TopLL)) then  ! keep CurPar(k)
    call NewPar(A, k, curPar(k))
  endif
endif 

do m = 1, 2
  if (CurPar(m)<0 .and. Parent(A,m)==0) then                             
    if (nS(-CurPar(m), m)==0) then ! no sibs left (was GP pair)
      call DoMerge(0, -CurPar(m), m, OK)
    endif                        
  endif
enddo
 call CalcLind(A)

end subroutine CalcPOZ

! #####################################################################

subroutine ReOrderAdd(LL)  
! reorder output from CheckAdd for compatibility withCheckPair (for POZ)
use Global
implicit none

double precision, intent(INOUT) :: LL(7)
double precision :: LLtmp(7)

LLtmp = 999D0
if (LL(4) < 0 .and. (LL(4) - MaxLL(LL(2:3))) > -TA) then  
  LLtmp(1) = 222
else
  LLtmp(1) = MaxLL(LL(2:3))
endif
LLtmp(2:3) = LL(5:6)
LLtmp(7) = LL(7) 

LL = LLtmp

end subroutine ReOrderAdd

! #####################################################################

subroutine FindPairs
use Global
use qsort_c_module
implicit none

logical :: UseAge, cPair
integer :: k, i, j, top, PairTypeTmp(XP*nInd), PairIDtmp(XP*nInd,2), mp(2)
double precision :: dLL, PairLLRtmp(XP*nInd), LL(7), LLg(7), LRS(2)
integer, allocatable, dimension(:) :: Rank
double precision, allocatable, dimension(:) :: SortBy

nPairs = 0
PairID = -9
PairDLLR = 999D0
PairType = 0
UseAge = AgePhase > 0

do i=1,  nInd-1 
  if (Error/=0) return
  if (MODULO(i,500)==0) call rchkusr()  
  if (ALL(Parent(i,:)/=0)) cycle
  do j=i+1,nInd
    if (hermaphrodites .and. ((ANY(parent(i,:)/=0) .and. ALL(parent(j,:)==0)) .or. &
     (ALL(parent(i,:)==0) .and. ANY(parent(j,:)/=0))))  cycle                                                       
    LRS = 0
    mp = 0                                                      
    call PairQHS(i, j, LRS(1))  ! quick check
    if (LRS(1) < TF) cycle  !true FS more likely to be HS than U
    do k=1,2
      if (Parent(i,k)/=0 .or. Parent(j,k)/=0) cycle
      if (Parent(i,k)==j .or. Parent(j,k)==i) cycle
      if (AgeDiff(i,j) /= 999 .and. UseAge) then
        if (AgePriorM(ABS(AgeDiff(i, j))+1, k) == 0.0)  cycle
      endif
      mp(k) = 1
    enddo
    if ((ALL(Parent(i,:)==0) .and. ALL(Parent(j,:)==0) .and. UseAge .and. &
      ALL(mp==1)) .or. (Hermaphrodites .and. (ALL(Parent(i,:)==0) .or. &
       ALL(Parent(j,:)==0)))) then
      call PairQFS(i, j, LRS(2)) 
    endif

    cPair = .FALSE.    
    do k=1,2
      if (mp(k)==0 .or. cPair)  cycle
      if (k==2 .and. mp(1)==1) then
        if (.not. UseAge) cycle
        if (LRS(2) < TF) cycle  
      endif    
      if (hermaphrodites .and. (ALL(Parent(i,:)==0) .or. ALL(Parent(j,:)==0))) then
        if (LRS(2) < TF) cycle 
      endif   
      if (AgeDiff(i,j)>=0) then
        call CheckPair(i, j, k, 3, LLg, LL)
      else
        call CheckPair(j, i, k, 3, LLg, LL)
      endif
      if (UseAge)  call BestRel(LL, 3, top, dLL)
      if (.not. UseAge)  call BestRel(LLg, 3, top, dLL)                          
      if (hermaphrodites .and. (ALL(Parent(i,:)==0) .or. &
       ALL(Parent(j,:)==0)) .and. top/=2)  cycle                                     
      if (top==2 .or. top==3) then  
        if (nPairs >= XP*nInd) cycle  ! do in next round
        nPairs = nPairs+1
        PairID(nPairs, :) = (/ i, j /)
        PairDLLR(nPairs) = dLL
        if (k==1 .and. mp(2)==1) then
          pairType(nPairs) = 3
          cPair = .TRUE.
        else
          PairType(nPairs) = k
        endif
      endif
    enddo                  
  enddo
enddo

! sort by decreasing dLL
PairIDtmp = 0
PairLLRtmp = 0
allocate(Rank(nPairs))
allocate(SortBy(nPairs))
Rank = (/ (i, i=1, nPairs, 1) /)
SortBy = PairDLLR(1:nPairs)
 
 call QsortC(SortBy, Rank(1:nPairs))
do i=1,nPairs
  PairTypeTmp(i) = PairType(Rank(nPairs-i+1))  ! decreasing order
  PairIDtmp(i,1:2) = PairID(Rank(nPairs-i+1), 1:2)  
  PairLLRtmp(i) = PairDLLR(Rank(nPairs-i+1)) 
enddo 

PairType = PairTypeTmp
PairID = PairIDtmp 
PairDLLR = PairLLRtmp
deallocate(Rank)
deallocate(SortBy)

end subroutine FindPairs

! #####################################################################

subroutine triads(nAmbMax, AmbigIDIN, topRel, ntrio, trioIDsOUT, trioLROUT, trioOHOUT)
use Global
implicit none

integer, intent(IN) :: nAmbMax, AmbigIDIN(2*nAmbMax), topRel(nAmbMax)
integer, intent(OUT) :: ntrio, trioIDsOUT(3*nInd), trioOHOUT(3*nInd)
double precision, intent(OUT) :: trioLROUT(3*nInd)
integer :: i, j, u, v, ncp, CandPar(50), m, k, curPar(2), nAmb, &
 AmbigID(nAmbMax,2), trioIDs(nInd, 3), trioOH(nInd, 3) 
double precision :: LLA(7), LLAA(7), LLtmp(2,2,2), trioLR(nInd, 3)

AmbigID(:,1) = AmbigIDIN(1:nAmbMax)
AmbigID(:,2) = AmbigIDIN((nAmbMax+1) : (2*nAmbMax))
nAmb = COUNT(AmbigID(:,1)>0)

ntrio = 0
trioIDs = 0
trioLR = 999D0
do i=1, nInd
  if (Error/=0) return
  if (MODULO(i,500)==0)   call rchkusr()
  if (ntrio == nInd) exit
  if (ANY(Parent(i,:)/=0) .and. .not. Hermaphrodites) cycle
  if (ANY(Parent(i,:)>0) .and. Hermaphrodites) cycle                                                  
  ncp = 0
  CandPar = 0  
  if ((COUNT(AmbigID(:,1) == i .and. topRel == 1) + &
    COUNT(AmbigID(:,2) == i .and. topRel == 1)) < 2) cycle
  do j=1, nAmb
    if (topRel(j) /= 1)  cycle
    if (.not. ANY(AmbigID(j,:) == i))  cycle
    if (ncp == 50) exit
    do m=1,2
      if (AmbigID(j,m) == i) then
        if (AgeDiff(i, AmbigID(j,3-m))>0) then   
          ncp = ncp + 1
          CandPar(ncp) = AmbigID(j,3-m)
        endif
      endif
    enddo
  enddo

  if (ncp > 1) then
    do u=1, ncp-1
      do v=u+1, ncp
        if (Sex(CandPar(u))/=3 .and. Sex(CandPar(u))==Sex(CandPar(v))) cycle      ! hermaphrodites                                                               
        call NewPar(i, 1, CandPar(u))
        call NewPar(i, 2, 0)
        if (CandPar(v) == CandPar(u))  cycle  ! not sure why that should happen?
        call CalcPOZ(i, CandPar(v), .FALSE.)
        if (ANY(Parent(i,:)==0)) then
          call NewPar(i, 1, 0)
          call NewPar(i, 2, 0)
          cycle
        endif
        ntrio = ntrio +1
        trioIDs(ntrio,1) = i
        trioIDs(ntrio, 2:3) = Parent(i,:)
        
        call CalcOH(i, Parent(i,1), trioOH(ntrio,1))
        call CalcOH(i, Parent(i,2), trioOH(ntrio,2))
        call CalcTrioErr(i, trioOH(ntrio,3))
        
        curPar = Parent(i,:)
        call NewPar(i, 1, 0)
        call NewPar(i, 2, 0)                   
        do m=1,2  ! m=1: no opp. sex parent;  m=2: with opp. sex parent
          do k=1,2
            if (m==2) then  
              call NewPar(i, 3-k, CurPar(3-k))          
            endif
            call CheckPair(i, CurPar(k), k, 1, LLA, LLAA)
            LLtmp(1,k,m) = LLA(1)
            LLtmp(2,k,m) = MaxLL(LLA(2:7))
            if (m==2) then  
              call NewPar(i, 3-k, 0)            
            endif
          enddo
        enddo
        do k=1,2
          trioLR(ntrio,k) = LLtmp(1,k,1)-LLtmp(2,k,1)
        enddo
        if (Complx>0) then
          trioLR(ntrio,3) = MIN(LLtmp(1,1,2) -MaxLL((/LLtmp(2,1,2), &
            LLtmp(:,1,1)/)), LLtmp(1,2,2) -MaxLL((/LLtmp(2,2,2), &
            LLtmp(:,2,1)/)))
        else
          trioLR(ntrio,3) = trioLR(ntrio,1)
        endif 
        call NewPar(i, 1, 0)
        call NewPar(i, 2, 0)           

        if (ntrio == nInd) exit        
      enddo
      if (ntrio == nInd) exit
    enddo
  endif
enddo

trioIDsOUT(1:nInd) = trioIDs(:,1)
trioIDsOUT((nInd+1):(2*nInd)) = trioIDs(:,2)
trioIDsOUT((2*nInd+1):(3*nInd)) = trioIDs(:,3)

trioLROUT(1:nInd) = trioLR(:,1)
trioLROUT((nInd+1):(2*nInd)) = trioLR(:,2)
trioLROUT((2*nInd+1):(3*nInd)) = trioLR(:,3)

trioOHOUT(1:nInd) = trioOH(:,1)
trioOHOUT((nInd+1):(2*nInd)) = trioOH(:,2)
trioOHOUT((2*nInd+1):(3*nInd)) = trioOH(:,3)

end subroutine triads

! #####################################################################

subroutine PairQHS(A, B, LR)  !quick check, not conditioning on parents.
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LR
integer :: l
double precision :: PrL(nSnp)

PrL = 0D0
do l=1,nSnp
  if (Genos(l,A)==-9 .or. Genos(l,B)==-9) cycle
  PrL(l) = LOG10(PHS(Genos(l,A), Genos(l,B), l)) ! note: >0 for FS 
enddo
LR = SUM(PrL)

end subroutine PairQHS

! #####################################################################

subroutine PairQFS(A, B, LR)  !quick check, not conditioning on parents.
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LR
integer :: l
double precision :: PrL(nSnp)

PrL = 0D0
do l=1,nSnp
  if (Genos(l,A)==-9 .or. Genos(l,B)==-9) cycle
  PrL(l) = LOG10(PFS(Genos(l,A), Genos(l,B), l))  
enddo
LR = SUM(PrL)

end subroutine PairQFS

! #####################################################################

subroutine CheckPair(A, B, kIN, focal, LLg, LL) 
! joined LL A,B under each hypothesis
use Global
implicit none

integer, intent(IN) :: A,B,kIN, focal                                   
double precision, intent(OUT) :: Llg(7), LL(7)  ! PO,FS,HS,GG,FAU,HAU,U
integer :: x, cgp, k
double precision :: LLtmpA(2,3), LLGGP(5), LLCC, LRS, ALR(7), &
 LLHH(2,3), LLX(5), LLZ(6), LLC(7), LLP(5), LLFC, LLPK(3)

LLg = 999D0
LL = 999D0
LLGGP = 999D0
LRS = 999D0
LLCC = 999D0
LLZ = 999D0           
ALR = 999D0           
if (kIN==3) then
  k = 1
else
  k = kIN
endif

if (focal==1 .and. Sex(B)/=k .and. Sex(B)/=3) then  
  LLg(1) = 777
  LL(1) = 777
  return
endif

 call CalcU(A,k,B,k, LLg(7))
LL(7) = LLg(7)
 
do x=1,4
  if (focal == x) then
    call CalcAgeLR(A,k, B, k, k, x, .TRUE., ALR(x)) 
    if (ALR(x)==777 .or. ALR(x)<5*TF) then
      LLg(x) = 777    
    else
      if (focal==1)  call PairPO(A, B, k, focal, LLg(1))
      if (focal==2)  call PairFullSib(A, B, LLg(2)) 
      if (focal==3) then
        if (Parent(A, 3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
          call PairFullSib(A, B, LLg(2)) 
          LLg(3) = 777
        else
          call PairHalfSib(A, B, k, LLg(3)) 
        endif
      endif
      if (focal==4)  call PairGP(A, B, k, focal, LLg(4))
    endif
  endif
enddo

if (focal/=7) then                  
  if ((LLg(focal) == 777 .or. LLg(focal) - LLg(7) < TA) .and. focal/=1 .and. &
   .not. (focal==3 .and. LLg(2)<0 .and. LLg(2) - LLg(7) > TA)) then
    LL = LLg
    return
  endif
endif     

do x=1,4
  if (ALR(x) /= 999) cycle
  call CalcAgeLR(A,k, B,k, k, x, .TRUE., ALR(x))
  if (ALR(x) == 777) then
    LLg(x) = 777
    LL(x) = 777
  endif
enddo

if (LLg(1)==999) then
  if (AgeDiff(A,B)>0) then
    call PairPO(A, B, k, focal, LLg(1))   
  else if (focal/=1) then
    call PairPO(B, A, k, focal, LLg(1))
  endif
endif
if (LLg(2)==999)  call PairFullSib(A, B, LLg(2)) 
if (LLg(3)==999 .and. Complx>0) then
  if (Parent(A, 3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
    LLg(3) = 777
  else
    call PairHalfSib(A, B, k, LLg(3)) 
  endif
endif
  
if (LLg(4)==999) then  ! GP?
  call PairGP(A, B, k, focal, LLg(4))
endif

if (ALR(4)/=777) then   ! no ageprior for GGP
  if (AgeDiff(A,B)>=3) then
    call PairGGP(A, B, k, LLGGP(1))
    call PairGGP(A, B, 3, LLGGP(2))  ! double GGP indistinguishable from GP
  endif
  do x=1,3  ! hf
    call PairGA(A, B, k, x, LLGGP(2+x))   ! HS of GP actually 4th degree rel, but giving false pos. 
  enddo 
endif
if (focal==3 .and. MaxLL(LLg(2:3))>LLg(4) .and. Sex(B)/=3 .and. &
  .not. (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0)) then
  cgp = 0
  if (Parent(A,3-k)>0)  cgp = Parent(Parent(A,3-k),Sex(B))
  if (Parent(A,3-k)<0)  cgp = GpID(Sex(B), -Parent(A,3-k), 3-k)
  if (cgp == 0) then
    call PairGP(A, B, 3-k, focal, LLCC)
    LLg(4) = MaxLL((/ LLg(4), LLCC /))  ! Note: wrong ageprior 
  endif
endif

! FA/HA?
LLFC = 999D0
 call CalcAgeLR(A,k, B,k, 0, 6, .TRUE., ALR(6))                                       
if (ALR(6)/=777) then  ! FA & HA have same ageprior.
  LLtmpA = 999D0
  do x=1,3  ! mat, pat, FS
    call PairUA(A, B, k, x, LLtmpA(1,x))
    call PairUA(B, A, k, x, LLtmpA(2,x))
  enddo
  LLg(5) = MaxLL(LLtmpA(:,3))
  if (Complx>0)  LLg(6) = MaxLL(RESHAPE(LLtmpA(:,1:2), (/2*2/) ))
endif     
if (Complx==2 .and. (focal==2 .or. focal==3) .and. LLg(2)<0 .and. &
 Parent(A,3-k)==Parent(B,3-k) .and. &
 (MaxLL(LLtmpA(:,3)) - MaxLL(LLg(2:3)) > -TA)) then
  call FSHC(A, B, k, LLFC)
  if (LLFC > LLg(2) .and. LLFC<0)  LLg(2) = LLFC
endif

do x=1,4
  if (LLg(x) < 0) then
    LL(x) = LLg(x) + ALR(x)
  endif
enddo
do x=5,6  ! same age prior for full & half aunts/uncles
  if (LLg(x) < 0) then
    LL(x) = LLg(x) + ALR(6)
  endif
enddo

LLCC = 999D0
 call PairCC(A, B, k, LLCC) 

 LLg(6) = MaxLL( (/LLg(6), LLGGP, LLCC/) )  ! most likely 3rd degree
  LL(6) = MaxLL( (/LL(6), LLGGP, LLCC/) ) 

LLPK = 999D0
if (focal/=1 .and. AgeDiff(A, B)>0 .and. Sex(B)/=k .and. Parent(A,3-k)<=0) then  !sex=3-k or 3
  if (AgeDiff(A,B)/=999) then
    if (AgePriorM(AgeDiff(A, B)+1, 6+3-k)==0) then
      LLPK(1) = 777
    endif
  endif
  if (LLPK(1)/=777) then
    if (Parent(A,3-k)==0) then
      call PairPO(A, B, 3-k, focal, LLPK(1))
    else if (Parent(A,3-k)<0) then
      call AddParent(B, -Parent(A,3-k), 3-k, LLPK(2))
      if (LLPK(2)<0) then
        call CalcU(B,3-k,Parent(A,3-k),3-k, LLPK(3))
        LLPK(1) = LLPK(2) - LLPK(3) + LL(7)
      endif
    endif
    if (LLPK(1)<0 .and. (LLPK(1) > LLg(1) .or. LLg(1)>0)) then
      LLg(1) = LLPK(1)
      if (AgeDiff(A,B)/=999) then
        LL(1) = LLPK(1) + LOG10(AgePriorM(AgeDiff(A, B)+1, 6+3-k))
      else
        LL(1) = LLPK(1)
      endif
    endif
  endif
endif

LLX = 999D0
LLC = 999D0   
LLHH = 999D0      
if (Complx==2 .and. LL(2)<0 .and. (focal==2 .or. focal==3)) then
  if (LL(2) - LL(7) > TA) then  ! check if inbred FS (& boost LL) 
    call PairFSHA(A, B, k, LLZ(1))
    call PairFSHA(A, B, 3-k, LLZ(2))
    if (MaxLL(LLZ(1:2)) > LLg(2) .and. ANY(LLZ(1:2)<0)) then
      LLg(2) = MaxLL(LLZ(1:2))
      LL(2) = LLg(2) + ALR(2)
    endif
  endif
  if (MaxLL(LL(2:3)) - MaxLL(LL) > -TA) then
    if (Parent(A,3-k)<=0 .or. Parent(B,3-k)<=0) then  
      do x=1,3
        call PairHSHA(A, B, k, x, LLHH(1,x))
        call PairHSHA(B, A, k, x, LLHH(2,x))
      enddo
      if (ANY(LLHH<0) .and. LLg(2) - MaxLL((/LLHH(1,:), LLHH(2,:)/)) < TA) then  ! 2*TA
        LL(2) = 222
      endif
      if (LL(2)/=222) then  
         if (AgeDiff(A,B)>0 .and. Sex(B)/=k) call PairHSPO(A,B,LLX(1))  
         if (AgeDiff(B,A)>0 .and. Sex(A)/=k) call PairHSPO(B,A,LLX(2))
        if ((LLX(1)<0 .or. LLX(2)<0) .and. &
          (LLg(2) - MaxLL(LLX(1:2))) < TA) then 
          LL(2) = 222
        endif
      endif
    endif
    if (LL(2)/=222 .and. LLg(4)/=777 .and. (Parent(A,3-k)==0 .or. &
      Parent(B,3-k)==0)) then
      call PairHSGP(A, B,k, LLX(3))
      call PairHSGP(B, A,k, LLX(4))
      if (LLX(3)<0 .or. LLX(4)<0) then
        if ((LLg(2) -MaxLL(LLX(3:4)))<TA) then
          LL(2) = 222
        endif
        do x=3,4
          if (MaxLL(LLX(3:4)) > LLg(x)) then
            LLg(x) = MaxLL(LLX(3:4))
            LL(x) = LLg(x) + ALR(x)
          endif
        enddo
      endif
    endif
    if (LL(2)/=222) then
      if (Parent(A,3-k) < 0 .and. Parent(B,3-k)==0) then  
        call CheckAdd(B, -Parent(A,3-k), 3-k, .FALSE., LLC, 3)
        if (LLC(2)<0 .and. (LLC(2) - MaxLL(LLC)) < TA) then 
          LL(2) = 222
        endif
      else if (Parent(A,3-k)==0 .and. Parent(B,3-k)<0) then
        call CheckAdd(A, -Parent(B,3-k), 3-k, .FALSE., LLC, 3)
        if (LLC(2)>0 .or. (MaxLL(LLC) - LLC(2)) > TA) then  
          LL(2) = 222
          endif
        else if (Parent(A,3-k)/=0 .or. Parent(B,3-k)/=0) then
        call PairHalfSib(A, B, 3-k, LLX(5))  
        if (LLX(5) < 0 .and. (LLg(2) - LLX(5))<TA) then
          LL(2) = 222
        endif
      endif
    endif       
    if (LL(2)/=222 .and. ALL(Parent(A,:)==0) .and. ALL(Parent(B,:)==0)) then
      call pairFAHA(A, B, .FALSE., LLZ(5))
      call pairFAHA(B, A, .FALSE., LLZ(6))
      if (ANY(LLZ(5:6) < 0)) then
        if ((LLg(2) - MaxLL(LLZ(5:6))) < TA) then
          LL(2) = 222
        endif
        if (MaxLL(LLZ(5:6)) > LLg(5)) then
          LLg(5) = MaxLL(LLZ(5:6))
          LL(5) = LLg(5) + ALR(6)
        endif
      endif
    endif
    if (LL(2)/=222) then  ! check if GG in any way. can't be FS and GP
      LLC = 999D0
      call PairGP(A, B, 3-k, focal, LLC(1))
      if (AgeDiff(A,B)==999) then
        call PairGP(B, A, k, focal, LLC(2))
        call PairGP(B, A, 3-k, focal, LLC(3))
      endif
      if (MaxLL(LLC(1:3))<0 .and. (LLg(4) - MaxLL(LLC(1:3))) <TA) then
        LLg(4) = MaxLL(LLC(1:3))
        LL(4) = MaxLL(LLC(1:3))  ! TODO: ageprior
      endif
    endif
  endif
  if (LL(2)==222)  LLg(2) = 222
endif

if (Complx==2 .and. focal==3 .and. LL(3)<0 .and. &
 (MaxLL(LL)>=LL(3) .or. MaxLL(LL)==LL(2))) then
  call pairHSHAI(A, B, k, LLZ(3)) ! HS + inbr HA
  call pairHSHAI(B, A, k, LLZ(4))
  if (MaxLL(LLZ(3:4)) < 0 .and. MaxLL(LLZ(3:4)) > LLg(3)) then
    LLg(3) = MaxLL(LLZ(3:4))
    LL(3) = MaxLL(LLZ(3:4)) + ALR(3)
  endif
endif

LLP = 999D0
if ((Sex(B)==3 .or. AgeDiff(A,B)==999)  .and. focal/=1) then
  if (AgeDiff(A,B)>0 .and. Sex(B)/=2)  call PairPO(A, B, 1, 0, LLP(1)) 
  if (AgeDiff(A,B)>0 .and. Sex(B)/=1)  call PairPO(A, B, 2, 0, LLP(2)) 
  if (AgeDiff(B,A)>0 .and. Sex(A)/=2)  call PairPO(B, A, 1, 0, LLP(3))
  if (AgeDiff(B,A)>0 .and. Sex(A)/=1)  call PairPO(B, A, 2, 0, LLP(4))
  LLg(1) = MaxLL(LLP)
  LL(1) = LLg(1)
endif
if (focal==1 .and. LL(5) < LL(7)) then  ! check if FS of other parent
  call PairUA(A, B, 3-k, 3, LLP(5))
  if (LLP(5)<0 .and. (LL(5) - LLP(5)) < TA) then
    LLg(5) = LLP(5)
    if (AgeDiff(A,B)/=999) then
      LL(5) = LLP(5) + LOG10(AgePriorM(ABS(AgeDiff(A, B))+1, 6))
    else
      LL(5) = LLP(5)
    endif
  endif
endif
 
end subroutine CheckPair

! #####################################################################

subroutine PairSelf(A, B, LL)  ! currently only called w/o parents
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: l, x
double precision :: PrX(3), PrL(nSnp)

PrL = 0D0
do l=1,nSnp
  if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
  do x=1,3
    PrX = AHWE(x, l)
    if (Genos(l,A)/=-9) PrX(x) = PrX(x) * OcA(Genos(l,A), x)
    if (Genos(l,B)/=-9) PrX(x) = PrX(x) * OcA(Genos(l,B), x)
  enddo
  PrL(l) = LOG10(SUM(PrX))
enddo
LL = SUM(PrL)

end subroutine PairSelf

! #####################################################################

subroutine PairPO(A, B, k, focal, LL)
use Global
implicit none

integer, intent(IN) :: A, B, k, focal
double precision, intent(OUT) :: LL
integer :: l, x, y,m, curPar(2), AncB(2,mxA), PAB, AncPA(2,mxA)
double precision :: PrL(nSnp,5), PrX(3,3,5), PrPA(3), PrB(3),PrPB(3,2),&
  LLtmp(5), PrPAB(3), PrPAX(3), PrG(3), LLX(4)
logical :: Maybe(5)                 

LL = 999D0
if(Parent(A,k)>0) then  ! allow dummy to be replaced (need for AddFS)
  if (Parent(A,k)==B) then
    LL = 888
  else
    LL = 777
  endif
endif

if (Sex(B)/=3 .and. Sex(B)/=k) then
  LL = 777
endif
if (LL/=999) return

Maybe = .FALSE.  ! 1: non-inbred, 2: B PO & GP; 3: B PO & HS, 4: ? , 5: B selfing
Maybe(1) = .TRUE.    
if (hermaphrodites)  Maybe(5) = .TRUE.           

if (Complx==2) then
  if (Parent(A,3-k)==0) then
    Maybe(2) = .TRUE.
  else if (Parent(A,3-k)>0) then
    if (Parent(Parent(A,3-k),k) == B .or. Parent(Parent(A,3-k),k)==0) then
      Maybe(2) = .TRUE.
    endif
  else if (Parent(A,3-k)<0) then
    if (GpID(k, -Parent(A,3-k), 3-k) == B .or. GpID(k, -Parent(A,3-k), 3-k) == 0) then
      Maybe(2) = .TRUE.
    endif
  endif
endif
 call GetAncest(Parent(A,3-k), 3-k, AncPA)                                                                                 
if (Complx==2 .and. (Parent(A,3-k)==Parent(B,3-k) .or. Parent(A,3-k)==0 &
  .or. Parent(B,3-k)==0)) then
    if (focal == 1 .and. Parent(A,3-k)==0 .and. Parent(B,3-k)/=0) then
      Maybe(3) = .FALSE.
    else if (Parent(A,3-k)/=0) then
      if (ANY(AncPA == B)) then
      Maybe(3) = .FALSE.
    else
      Maybe(3) = .TRUE. 
    endif
  else
    Maybe(3) = .TRUE.
  endif
endif
if (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0) then
  Maybe(1) = .FALSE.  ! becomes config (3)
else if (Parent(A,3-k) > 0) then
  if (any(Parent(Parent(A,3-k), :) == B)) then
    Maybe(1) = .FALSE.  ! becomes config (2). TODO: use AncPA above.
  endif
else if (Parent(A,3-k) < 0) then
  if (any (GpID(:, -Parent(A,3-k), 3-k) == B)) then
    Maybe(1) = .FALSE.
  endif
endif

if (ANY(Parent(B,:)/=0)) then
  call getAncest(B, k, AncB)
  if (ANY(AncB == A)) then  ! if agediff unknown
    LL = 777
    return
  else if (ANY(AncB(3-k, 3:4) == Parent(A,3-k))) then
    Maybe(4) = .TRUE.
  endif
endif

PAB = Parent(A,3-k)
if(Maybe(3) .and. Parent(A,3-k)==0) then  ! B PO & HS
  PAB = Parent(B,3-k)
endif            

  PrL = 0D0
  LLtmp = 999
  do l=1,nSnp
    PrX = 0       
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
    PrB = LindX(:,l,B)
    do m=1,2  
      call ParProb(l, Parent(B,m), m, B, 0, PrPB(:,m))
    enddo
    if (Maybe(3)) then
      call ParProb(l, PAB, 3-k, A, B, PrPAB)
    endif
    if (Maybe(2)) then
      call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX)
      if (Parent(A,3-k)>0) then
        call ParProb(l, Parent(Parent(A,3-k),3-k), 3-k, Parent(A,3-k), 0, PrG)
      else if (Parent(A,3-k)<0) then
        call ParProb(l, GpID(3-k,-Parent(A,3-k),3-k), 3-k, 0, 0, PrG)
      else
        PrG = AHWE(:,l)
        PrPAX = 1         
      endif
    endif
    if (Genos(l,A)==-9) then
      PrL(l,1) = LOG10(SUM(PrB))
      if (Parent(A,3-k)==0)  PrL(l,2) = LOG10(SUM(PrB))
      if (Parent(A,3-k)==Parent(B,3-k))  PrL(l,3) = LOG10(SUM(PrB))
      cycle
    endif
    do x=1,3  ! B
      do y=1,3  ! parent(A,3-k)
        PrX(x,y,1) = OKA2P(Genos(l, A), x, y) * PrB(x) * PrPA(y)
        if (Maybe(2)) then  ! B PO & GP
          PrX(x,y,2) = OKA2P(Genos(l, A),x,y)* PrB(x)* PrPAX(y) *&
            SUM(AKA2P(y,x,:) * PrG)                         
        endif
        if (Maybe(3)) then  ! B PO & HS
          PrX(x,y,3) = OKA2P(Genos(l, A), x, y) * PrPAB(y) * &
            SUM(AKA2P(x,y,:) * PrPB(:,k))
          if (Genos(l,B)/=-9) then
            PrX(x,y,3) = PrX(x,y,3) * OcA(Genos(l,B), x)
          endif
        endif
        if (hermaphrodites) then  ! B selfing
          if (x/=y)  PrX(x,y,5) = 0.0
          if (x==y)  PrX(x,y,5) = OKA2P(Genos(l, A), x, x) * PrB(x)                                                    
        endif
      enddo
    enddo
    do x=1,5
      if (Maybe(x))  PrL(l,x) = LOG10(SUM(PrX(:,:,x)))
    enddo                                                         
  enddo
  
  LLtmp = SUM(PrL, DIM=1)
  if (Parent(A,3-k) > 0) then
    LLtmp(2) = LLtmp(2) - Lind(Parent(A,3-k))
  endif
  do x=1,5
    if (.not. Maybe(x) .or. LLtmp(x)>=0)  LLtmp(x) = 777
  enddo
   if (PAB<0 .and. ANY(Maybe(2:4))) then 
    curPar = Parent(A, :)
    call NewPar(A, k, 0)     ! B vs none.    
    call CalcU(A,3-k, B,k, LLX(1))
    call CalcU(Parent(A,3-k),3-k, B,k, LLX(2))
    
    call NewPar(A, k, B)   
    call CalcU(Parent(A,3-k),3-k, B,k, LLX(3))
    LLtmp(1) = LLX(1) + (LLX(3) - LLX(2))
    
    if (Maybe(2) .and. Parent(A,3-k)<0 .and. .not. ANY(AncPA == B).and. &
      .not. ANY(AncB(3-k,:)==Parent(A,3-k))) then   ! B PO & GP
      GpID(k, -Parent(A,3-k), 3-k) = B
      call CalcCLL(-Parent(A,3-k), 3-k)
      call CalcU(Parent(A,3-k),3-k, B,k, LLX(4))
      LLtmp(2) = LLX(1) + (LLX(4) - LLX(2))
      GpID(k, -Parent(A,3-k), 3-k) = 0
      call CalcCLL(-Parent(A,3-k), 3-k)
    endif
    
    call NewPar(A, k, curPar(k))  ! restore
  endif                          

  LL = MaxLL(LLtmp)

end subroutine PairPO

! #####################################################################

subroutine PairFullSib(A, B, LL)
use Global
implicit none

integer, intent(IN) :: A, B
double precision, intent(OUT) :: LL
integer :: x, y, l,k, Par(2), AncA(2,mxA), AncB(2,mxA), ix, i
double precision :: PrL(nSnp), PrXY(3,3), Px(3,2), LUX(2), LLtmp, &
  dx(maxSibSize)            

LL = 999D0
Par = 0  ! joined parents of A & B
if (Parent(A,1)==Parent(B,1) .and. Parent(A,1)/=0 .and. &
  Parent(A,2)==Parent(B,2) .and. Parent(A,2)/=0) then ! already FS
  LL = 888
    return
else 
  do k=1,2
    if (Parent(A,k) == B .or. Parent(B,k) == A) then
      LL = 777
      return
    else if (Parent(A,k)/=Parent(B,k) .and. .not. (Parent(A,k)==0 .or. &
      Parent(B,k)==0)) then
      LL = 777
      return
    else if (Parent(A,k)/=0) then
      Par(k) = Parent(A,k)
    else
      Par(k) = Parent(B,k)
    endif        
  enddo
endif  

 call GetAncest(A,1,AncA)
 call GetAncest(B,1,AncB)
if (ANY(AncA == B) .or. ANY(AncB == A)) then
  LL = 777
  return
endif
   
PrL = 0D0 
LUX = 0
if (Par(1) < 0 .and. Par(2)<0) then  ! call AddFS
  call CalcU(A, 1, B, 1, LUX(1))
  do k=1,2 
    if (Parent(A,k)==Par(k) .and. Parent(B,k)==0) then
      call addFS(B, -Par(k), k, 0, k, LLtmp, ix, dx)
      do i=1, nS(-Par(k),k)
        if (SibID(i,-Par(k),k) == A) then
          LL = LUX(1) + dx(i)
          return
        endif
      enddo
    else if (Parent(B,k)==Par(k) .and. Parent(A,k)==0) then
      call addFS(A, -Par(k), k, 0, k, LLtmp, ix, dx)
      do i=1, nS(-Par(k),k)
        if (SibID(i,-Par(k),k) == B) then
          LL = LUX(1) + dx(i)
          return
        endif
      enddo      
    endif
  enddo     
else  
  do l=1, nSnp
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
    do k=1,2
      if (Parent(A,k)==Parent(B,k)) then  
        call ParProb(l, Par(k), k, A, B, Px(:,k))
      else if (Parent(A,k)==Par(k)) then
        call ParProb(l, Par(k), k, A, 0, Px(:,k))
      else if (Parent(B,k)==Par(k)) then
        call ParProb(l, Par(k), k, B, 0, Px(:,k))
      else
        call ParProb(l, Par(k), k, 0, 0, Px(:,k))
      endif       
    enddo 
  
    do x=1,3
      do y=1,3
        PrXY(x,y) = Px(x,1) * Px(y,2)
        if(Genos(l,A)/=-9) then
          PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
        endif
        if(Genos(l,B)/=-9) then
          PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,B), x, y)
        endif
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
  enddo
  LL = SUM(PrL) 
endif

end subroutine PairFullSib

! #####################################################################

subroutine PairHalfSib(A, B, k, LL)
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: x,y, l, Par, Inbr, AB(2), AncA(2,mxA), AncB(2,mxA)
double precision :: PrL(nSnp), PrX(3), PrPx(3,2), PrXY(3,3)

LL = 999D0
Par = 0  ! parent K
if (Parent(A,k)/=0) then
  if (Parent(A,k)/=Parent(B,k) .and. Parent(B,k)/=0) then
    LL = 777 ! mismatch
  else if (Parent(A,k)==Parent(B,k)) then
    LL = 888 ! already captured under H0
  else
    Par = Parent(A,k)
    if (Par>0) then
      if (AgeDiff(B, Par) <= 0) then  ! Par(k) younger than B
        LL = 777
      endif
    endif
  endif
else if (Parent(B,k)/=0) then
  Par = Parent(B,k)
  if (Par>0) then
    if (AgeDiff(A, Par) <= 0) then  ! Par(k) younger than A
      LL = 777
    endif
  endif
endif
if (LL/=999) return

if (Parent(A,k)/=0) then
  call GetAncest(Parent(A,k),k,AncA)
  if (ANY(AncA == B)) LL = 777 
endif
if (Parent(B,k)/=0) then
  call GetAncest(Parent(B,k),k,AncB)
  if (ANY(AncB == A)) LL = 777 
endif
if (LL/=999) return
 
AB = (/ A, B /)
Inbr = 0
if (Parent(A,3-k)==B)  Inbr = 1
if (Parent(B,3-k)==A)  Inbr = 2
PrL = 0D0
do l=1,nSnp
   if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
  if (Par==Parent(A,k) .and. Par/=0) then
    call ParProb(l, Par, k, A, 0, PrX)
  else if (Par==Parent(B,k) .and. Par/=0) then
    call ParProb(l, Par, k, B, 0, PrX)
  else
    call ParProb(l, Par, k, 0, 0, PrX)    
  endif
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPx(:,1))
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPx(:,2))
  if (Inbr==0) then
    do x=1,3
      do y=1,3 
        PrXY(x,y) = PrX(x) * PrPX(y,1)
        if (Genos(l,A)/=-9) then
          PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A),x,y)
        endif
        if (Genos(l,B)/=-9) then
          PrXY(x,y) = PrXY(x,y) * SUM(OKA2P(Genos(l,B),x,:) * PrPX(:,2))
        endif
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
  else 
    do x=1,3
      do y=1,3
        PrXY(x,y)=PrX(x) * SUM(AKA2P(y, x, :) * PrPx(:,3-Inbr))
        if (Genos(l,AB(Inbr))/=-9) then
          PrXY(x,y)= PrXY(x,y) * OKA2P(Genos(l,AB(Inbr)), x, y)                                   
        endif
        if (Genos(l,AB(3-Inbr))/=-9) then
          PrXY(x,y) = PrXY(x,y) * OcA(Genos(l, AB(3-Inbr)), y)
        endif
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXY))
  endif
enddo
LL = SUM(PrL)

end subroutine PairHalfSib

! #####################################################################

subroutine pairHSHA(A, B, k, hf, LL)  !HS via k, & parent A is HS of B via 3-k
use Global
implicit none

integer, intent(IN) :: A,B, k, hf
double precision, intent(OUT) :: LL
integer :: l, x, y, z, PAB
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrY(3), PrZ(3)

LL = 999D0
PAB = 0
if (Parent(A,3-k)/=0) then
  LL = 444
  return
else if (Parent(A,k)/=Parent(B,k)) then
  if (Parent(A,k)/=0) then
    if(Parent(B,k)/=0) then
      LL = 777
      return
    else
      PAB = Parent(A,k)
    endif
  else
    PAB = Parent(B,k)
  endif
endif    

PrL = 0D0
do l=1, nSnp
  if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrY)
  if (Parent(A,k)==Parent(B,k)) then
    call ParProb(l, PAB, k, A, B, PrZ)
  else if (Parent(A,k) == PAB) then
    call ParProb(l, Parent(A,k), k, A, 0, PrZ)
  else if (Parent(B,k) == PAB) then
    call ParProb(l, Parent(B,k), k, B, 0, PrZ)
  endif
  
  do x=1,3
    do y=1,3    
      do z=1,3
        if (hf==1) then
          PrXYZ(x,y,z) = PrY(y) * PrZ(z) * AKAP(x, y, l)
        else if (hf==2) then
          PrXYZ(x,y,z) = PrY(y) * PrZ(z) * AKA2P(x, y, z)
        else if (hf==3) then
          PrXYZ(x,y,z) = PrY(y) * PrZ(z) * AKAP(x, z, l)
        endif
        if (Genos(l,B)/=-9) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), y,z)
        endif
        if (Genos(l,A)/=-9) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,A), x,z)
        endif
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine pairHSHA

! #####################################################################

subroutine pairHSHAI(A, B, k, LL)  !HS via k, & A inbred
use Global
implicit none

integer, intent(IN) :: A, B, k
double precision, intent(OUT) :: LL
integer :: l, x, y, AncA(2, mxA)
double precision :: PrL(nSnp,2), PrXY(3,3,2), PrPA(3), PrPB(3), LLU

if (Parent(A,k)/=0 .or. Parent(B,k)/=0) then
  LL = 777
endif 
if (Parent(A, 3-k)==B) LL = 777
if (LL==777) return 
 call GetAncest(A,1,AncA)
if (ANY(AncA == B)) LL = 777
if (LL==777) return 
if (Parent(A,3-k)==Parent(B,3-k) .and. Parent(A,3-k)/=0) then
    LL = 444  ! likely picked up elsewhere
    return
endif

PrL = 0D0
do l=1, nSnp
  if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPB)
  do x=1,3
    do y=1,3    
        PrXY(x,y,:) = AKAP(x, y, l) * PrPA(y)
        if (Genos(l,A)/=-9) then
          PrXY(x,y,:) = PrXY(x,y,:) * OKA2P(Genos(l,A), x,y)
        endif
        if (Genos(l,B)/=-9) then
          PrXY(x,y,1) = PrXY(x,y,1) * SUM(OKA2P(Genos(l,B), x,:)*PrPB)
          PrXY(x,y,2) = PrXY(x,y,2) * SUM(OKAP(Genos(l,B),:,l)*PrPB)
        endif
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrXY(:,:,1)))
  PrL(l,2) = LOG10(SUM(PrXY(:,:,2)))
enddo
 call CalcU(A,k,B,k, LLU)
if (SUM(PrL(:,2)) > LLU) then
  LL = SUM(PrL(:,1)) - SUM(PrL(:,2)) + LLU
else
  LL = SUM(PrL(:,1))
endif

end subroutine pairHSHAI

! #####################################################################

subroutine pairFAHA(A, B, withFS, LL)  !B FA via k & HA via 3-k
use Global
implicit none

integer, intent(IN) :: A, B
logical, intent(IN) :: withFS
double precision, intent(OUT) :: LL
integer :: l, k, x, y, z, v, i, AA(maxSibSize), BB(maxSibSize), nA, nB
double precision :: PrL(nSnp, 2), PrXY(3,3,3,3, 2), PrP(3), LLU

LL = 999D0
if (ANY(Parent(A,:)>0) .or. ANY(Parent(B,:)>0)) then  
  LL = 444
  return
endif

AA = 0
BB = 0
if (withFS) then
  nA = nFS(FSID(maxSibSize+1, A))  
  AA(1:nA) = FSID(1:nA, FSID(maxSibSize+1, A))
  nB = nFS(FSID(maxSibSize+1, B))
  BB(1:nB) = FSID(1:nB, FSID(maxSibSize+1, B))
else
  nA = 1
  AA(1) = A
  nB = 1
  BB(1) = B
endif

PrL = 0D0
do l=1, nSnp
  if (Genos(l,A)==-9 .and. Genos(l,B)==-9 .and. .not. withFS) cycle
  call ParProb(l, 0, k, 0, 0, PrP)
  do x=1,3
    do y=1,3    
      do z=1,3
        do v=1,3
          PrXY(x,y,z,v,1) = PrP(y) * PrP(z) * AKAP(v,y,l) * AKA2P(x,y,z)
          PrXY(x,y,z,v,2) = PrP(y) * PrP(z) * PrP(v) * PrP(x)  ! A, B unrelated
          do i=1, nA
            if (Genos(l,AA(i))/=-9) then
              PrXY(x,y,z,v,:) = PrXY(x,y,z,v,:) * OKA2P(Genos(l,AA(i)), x, v)
            endif
          enddo
          do i=1, nB
            if (Genos(l,BB(i))/=-9) then
              PrXY(x,y,z,v,:) = PrXY(x,y,z,v,:) * OKA2P(Genos(l,BB(i)), y, z)
            endif
          enddo
        enddo
      enddo
    enddo                          
  enddo
  PrL(l,1) = LOG10(SUM(PrXY(:,:,:,:,1)))
  PrL(l,2) = LOG10(SUM(PrXY(:,:,:,:,2)))
enddo
if (.not. withFS) then
  LL = SUM(PrL(:,1))
else
  call CalcU(A, 0, B, 0, LLU)
  LL = SUM(PrL(:,1)) - SUM(PrL(:,2)) + LLU
endif

end subroutine pairFAHA

! #####################################################################

subroutine pairHSPO(A, B, LL)   ! HS via k, & PO via 3-k
use Global
implicit none

integer, intent(IN) :: A,B
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), PrXY(3,3)

if (ANY(Parent(A,:)/=0) .or. ANY(Parent(B,:)/=0)) then
  LL = 777
  return   ! else not necessary.
endif  

PrL = 0D0
do l=1, nSnp
  if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
  do x=1,3 
    do y=1,3    ! B
      PrXY(x,y) = AHWE(x,l) * AKAP(y,x,l)
      if (Genos(l,B)/=-9) then
        PrXY(x,y) = PrXY(x,y) * OcA(Genos(l,B),y)
      endif
      if (Genos(l,A)/=-9) then
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
      endif
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine pairHSPO

! #####################################################################

subroutine clustHSHA(SA, SB, k, LL)   ! HS via 3-k, & SB parent of SA; SA,SB FS
use Global
implicit none

integer, intent(IN) :: SA,SB, k
double precision, intent(OUT) :: LL
integer :: l, x, y, z,i, Par(2), GC(2), u
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrGA(3),PrGC(3,2),PrUZ(3,3)

! all checks done by CheckMerge.

! grandparents of opp. parent
 call getFSpar(SA, k, .TRUE., Par(1))
 call getFSpar(SB, k, .TRUE., Par(2))
GC = 0
do i=1,2
  if (Par(1)<0) then
    GC(i) = GpID(i, -Par(1),3-k)
    if (GpID(i, -Par(1),3-k)/=0) then
      if (Par(2) < 0) then
       if (GpID(i,-Par(2),3-k)/=GC(i) .and. GpID(i,-Par(2),3-k)/=0) then
          GC(i) = 0   ! shouldn't happen
        else if (GC(i)==0 .and. GpID(i, -Par(2),3-k)/=0) then
          GC(i) = GpID(i, -Par(2),3-k)
        endif
      endif
    endif
  endif
enddo

LL = 999D0
PrL = 0D0
do l=1, nSnp
  call ParProb(l, GpID(3-k,SA,k), 3-k, 0, 0, PrGA)
  do i=1,2
    call ParProb(l, GC(i), i, 0, 0, PrGC(:,i))
  enddo
  do z=1,3
    do u=1,3
      PrUZ(u,z) = SUM(AKA2P(z,u,:) * PrGC(u,1) * PrGC(:,2))
    enddo
    do x=1,3    
      do y=1,3
        PrXYZ(x,y,z) = SUM(AKA2P(x,y,:) * PrGA) * XPr(2,y,l,SB,k) *&
         SUM(PrUZ(:,z))
        do i=1,nS(SA,k)
          if (Genos(l, SibID(i,SA,k))/=-9) then
           PrXYZ(x,y,z) =PrXYZ(x,y,z) *OKA2P(Genos(l,SibID(i,SA,k)),x,z)
          endif
        enddo 
        do i=1,nS(SB,k)
          if (Genos(l, SibID(i,SB,k))/=-9) then
           PrXYZ(:,y,z) =PrXYZ(:,y,z) *OKA2P(Genos(l,SibID(i,SB,k)),y,z)
          endif
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine clustHSHA

! #####################################################################

subroutine FSHC(A, B, k, LL)  ! FS + parents are HS; B may be neg
use Global
implicit none

integer, intent(IN) :: A,B, k
double precision, intent(OUT) :: LL
integer :: l, x, y, z, Par(2), m, GG(2,2), kG, i, PM(2)
double precision :: PrL(nSnp), PrG(3,2), PrXYZ(3,3,3), PrZ(3)

LL = 999D0
if (B < 0 .and. A>0) then
    Par(k) = B
    if (Parent(A,k)/=0) then
        LL = 777
        return
    endif
    if (ALL(Parent(SibID(1:nS(-B,k), -B, k), 3-k) == 0)) then
        Par(3-k) = Parent(A, 3-k)
    else
        call getFSpar(-B, k, .TRUE., Par(3-k))
        if (Par(3-k)==0 .or. (Parent(A,3-k)/=Par(3-k) .and. &
          Parent(A, 3-k)/=0)) then
            LL = 777
            return
        endif
    endif
else if (B > 0 .and. A>0) then
    do m=1,2
        if (Parent(B,m)==0) then
            Par(m) = Parent(A,m)
        else if (Parent(B,m) /= Parent(A,m) .and. Parent(A,m)/=0) then
            LL = 777
            return
        else
            Par(m) = Parent(B,m)
        endif
    enddo
else if (B<0 .and. A<0) then
    if (ANY(GpID(:,-B,k)/=0) .or. ANY(GpID(:,-A,k)/=0)) then
        LL = 444
        return      
    else
        Par = 0
    endif
endif

GG = 0
kG = 0
PM = 0
do m=1,2
  if (Par(m)>0)  GG(:, m) = Parent(Par(m), :)
  if (Par(m)<0)  GG(:, m) = GpID(:, -Par(m), m)
enddo
do m=1,2
  if (GG(m,1)==0 .or. GG(m,2)==0) then  ! GG(m,1)==GG(m,2) not needed
    kG = m
  endif
  if (Par(m)>0) then
    PM(m) = Par(m)
  endif
enddo       
if (kG==0) then
    LL = 888
    return
endif

PrL = 0D0
do l=1, nSnp
  if (B>0 .and. A>0) then
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
  endif
  do m=1,2
    call ParProb(l, GG(3-kG, m), 3-kG, PM(m),0, PrG(:,m))
  enddo
  if (GG(kG,1)/=0) then
    call ParProb(l, GG(kG,1), kG, PM(1),0, PrZ)
  else
    call ParProb(l, GG(kG,2), kG, PM(2),0, PrZ)
  endif
  do x=1,3  ! Par(1)
    do y=1,3  !Par(2)
      do z=1,3  ! GG(kG)
        PrXYZ(x,y,z) = PrZ(z) * SUM(AKA2P(x, z, :) * PrG(:,1)) * &
          SUM(AKA2P(y, z, :) * PrG(:,2))
        if (A>0) then
            if (Genos(l,A)/=-9) then
              PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,A), x, y)
            endif
        else
          do i=1, nS(-A,k)
            if (Genos(l,SibID(i,-A,k))/=-9) then
              PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,SibID(i,-A,k)), x, y)
            endif
          enddo
        endif
        if (B>0) then
          if (Genos(l,B)/=-9) then
            PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), x, y)
          endif
        else
          do i=1, nS(-B,k)
            if (Genos(l,SibID(i,-B,k))/=-9) then
              PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,SibID(i,-B,k)), x, y)
            endif
          enddo
        endif
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine FSHC

! #####################################################################

subroutine pairFSHA(A, B, k, LL) !inbred FS: par k offspring of par 3-k
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), PrXY(3,3), PrY(3)

if (Parent(A,k)/=0 .or. Parent(B,k)/=0) then
  LL = 444   ! TODO (prob. not necessary)
  return
endif  

PrL = 0D0
do l=1, nSnp
  if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
  call ParProb(l, Parent(A,3-k), 3-k, -1,0, PrY) 
  do x=1,3
    do y=1,3    
      PrXY(x,y) = PrY(y) * AKAP(x, y, l)
      if (Genos(l,B)/=-9) then
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,B), x, y)
      endif
      if (Genos(l,A)/=-9) then
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
      endif
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine pairFSHA

! #####################################################################

 subroutine trioFS(A, B, C, LL)
use Global
implicit none

integer, intent(IN) :: A,B,C
double precision, intent(OUT) :: LL
integer :: l, x, y
double precision :: PrL(nSnp), PrXY(3,3)

!if (any(Parent(A,:)/=0) .or. any(Parent(B,:)/=0) .or. any(Parent(C,:)/=0)) then
!  LL = 444   ! not necessary
!  return
!endif

PrL = 0D0
do l=1, nSnp
  do x=1,3
    do y=1,3    
      PrXY(x,y) = AHWE(x,l) * AHWE(y,l)   ! parents are ignored
      if (Genos(l,A)/=-9) then
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
      endif
      if (Genos(l,B)/=-9) then
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,B), x, y)
      endif
      if (Genos(l,C)/=-9) then
        PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,C), x, y)
      endif
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine trioFS

! #####################################################################

subroutine trioFA(A, B, C, LL)  ! B & C both FS of par(A,k)
use Global
implicit none

integer, intent(IN) :: A,B,C
double precision, intent(OUT) :: LL
integer :: l, x, y, z
double precision :: PrL(nSnp), PrXYZ(3,3,3)

 ! parents are ignored 

PrL = 0D0
do l=1, nSnp
  do x=1,3
    do y=1,3    
      do z=1,3 
        PrXYZ(x,y,z) = AKA2P(x,y,z) * AHWE(y,l) * AHWE(z,l)  
        if (Genos(l,A)/=-9) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKAP(Genos(l,A), x, l)
        endif
        if (Genos(l,B)/=-9) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,B), y, z)
        endif
        if (Genos(l,C)/=-9) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,C), y,z)
        endif
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine trioFA

! #####################################################################

subroutine pairHSGP(A, B,k, LL)   ! HS via k, B is GP of A via 3-k
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y, z
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrPB(3)

if (Parent(A,3-k)/=0) then
  LL = 444
  return
endif 

PrL = 0D0
do l=1, nSnp
  if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
  call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPB)
  do x=1,3  ! parent 3-k of A, offspring of B
    do y=1,3  ! shared parent k 
      do z=1,3  ! B
        PrXYZ(x,y,z) =AKAP(x,z,l)*AHWE(y,l)*SUM(AKA2P(z,y,:)*PrPB)
        if (Genos(l,A)/=-9) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OKA2P(Genos(l,A),x, y)
        endif
        if (Genos(l,B)/=-9) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * OcA(Genos(l,B), z)
        endif
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo
LL = SUM(PrL)

end subroutine pairHSGP

! #####################################################################

subroutine PairGP(A, B, k, focal, LL)  
! calculates LL that B is maternal(k=1) or paternal(k=2) gp of A
use Global
implicit none

integer, intent(IN) :: A,B,K, focal
double precision, intent(OUT) :: LL
integer :: l, x, y, curGP(2), m, z, AncB(2, mxA), v, AncPA(2,mxA)
double precision :: PrL(nSnp,3), PrPA(3,2), PrG(3), LLtmp(3),&
   PrXZ(3,3,3,3), PrB(3), PrGx(3), PrPB(3), PrV(3), PrPAX(3,2), ALR
logical :: cat(3)                

LL = 999D0
 curGP = 0  
if (Parent(A,k)>0) then
  curGP = Parent(Parent(A,k),:) ! current grandparents of A (shortcut)
else if (Parent(A,k)<0) then
  curGP = GpID(:, -Parent(A,k), k)
endif

if (AgeDiff(A, B) <= 0) then  ! B younger than A
  LL = 777
  return      
else if (Sex(B)/=3) then
  if (curGP(Sex(B))>0 .or. (curGP(Sex(B))/=0 .and. focal==4)) then
    if (curGP(Sex(B))/=B) then
      LL = 777  ! conflict
    endif
  endif
  m = Sex(B)
else 
  if (curGP(1)/=0 .and. curGP(2)/=0) then
    do y=1,2
      if (curGP(y) == B) then
        LL = 888
        exit
      else if (focal==4 .or. (curGP(1)>0 .and. curGP(2)>0)) then
        LL = 777
      endif
    enddo
  endif
  if (curGP(1)==0) then
    m = 1  ! doesn't really matter.
  else
    m = 2 
  endif
endif

LLtmp = 999D0
if (Parent(A,k)>0 .and. LL==999) then
  if (AgeDiff(Parent(A,k), B) <= 0) then  ! B younger than Parent(A,k)
    LL = 777 
  else
    if (Sex(B)/=3) then
      call PairPO(Parent(A,k), B, Sex(B), focal, LLtmp(1))
    else
      call PairPO(Parent(A,k), B, 1, focal, LLtmp(1))
    endif
    if (LLtmp(1) > 0) then    ! impossible
      LL = 777
    else 
      call CalcU(Parent(A,k), k, B,k, LLtmp(2))
      if (LLtmp(1) - LLtmp(2) < TA) then
        LL = 777
      endif
    endif
  endif
endif
if (LL/=999) return

 call GetAncest(B, k, AncB)
if (ANY(AncB == A)) then
  LL = 777
endif
if (LL/=999) return

 cat = .TRUE.
if (Parent(B,3-k)==Parent(A,3-k) .and. Parent(A,3-k)/=0) then
  cat(1) = .FALSE.
endif
if (complx < 2) then
  cat(2) = .FALSE.
else if (Parent(A,3-k)/=0) then
  if (ANY(AncB(3-k,:) == Parent(A,3-k))) then
    cat(2) = .FALSE.
  else
    call CalcAgeLR(Parent(A,3-k), 3-k, B,k, k, 1, .TRUE., ALR)
    if (ALR == 777 .or. ALR < 3*TF) then
      cat(2) = .FALSE.
    endif 
  endif
endif

if (complx==2 .and. (Parent(B,3-k)==0 .or. Parent(A,3-k)==0 .or. & 
 Parent(B,3-k)==Parent(A,3-k))) then
  cat(3) = .TRUE.
  if (Parent(A,3-k)/=0 .and. Parent(B,3-k)==0) then
    call getAncest(Parent(A,3-k), 3-k, AncPA)
    if (ANY(AncPA == B)) then
      cat(3) = .FALSE.
    else
      call CalcAgeLR(B,k, Parent(A,3-k), 3-k, 3-k, 1, .TRUE., ALR)
      if (ALR == 777 .or. ALR < 3*TF) then
        cat(3) = .FALSE.
      endif
    endif
  endif
else
  cat(3) = .FALSE.
endif   

PrL = 0D0
do l=1,nSnp
  if (Genos(l,A)==-9) cycle
  call ParProb(l, curGP(3-m), 3-m, 0,0, PrG) 
  call ParProb(l, Parent(A,k), k, A, -4, PrPA(:, k)) 
  if (Parent(A,3-k)==Parent(B,3-k)) then
    call ParProb(l, Parent(A,3-k), 3-k, A, B, PrPA(:, 3-k))
  else
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA(:, 3-k))
  endif
  call ParProb(l, Parent(A,3-k), 3-k, A, -4, PrPAX(:, 3-k))                                                           
  call ParProb(l, B, 0, 0, 0, PrB)
  call ParProb(l, Parent(B,k), k, B, 0, PrPB)
  if (cat(2)) then
    if (Parent(A,3-k) > 0) then
      call ParProb(l, Parent(Parent(A,3-k),3-m), 3-m, Parent(A,3-k), 0, PrGx)
    else if (Parent(A,3-k) < 0) then
      call ParProb(l, GpID(3-m, -Parent(A,3-k),3-k), 3-m, 0, 0, PrGx)   
    else
      PrGx = AHWE(:,l)
    endif
  endif    

  PrXZ = 1
  do x=1,3  ! PA(k)
    if (Parent(A,k)/=0) then
      PrXZ(x,:,:,:) = PrPA(x, k)
    endif      
    do y=1,3  ! PA(3-k)
      do z=1,3  !  PrG(3-m)
        PrXZ(x,y,z,:)=PrXZ(x,y,z,:)*OKA2P(Genos(l,A),x,y)*PrG(z)
        if (cat(1)) then
          PrXZ(x,y,z,1) = PrXZ(x,y,z,1) * PrPA(y,3-k) *&
           SUM(AKA2P(x, :, z) * PrB)  !non-inbred
        endif
        if (cat(2)) then   !inbreeding loop; B double gp
          do v=1,3
            PrV(v) = SUM(AKA2P(y,v,:) * PrGx) * AKA2P(x,v,z) * PrB(v)
          enddo
          PrXZ(x,y,z,2) =PrXZ(x,y,z,2) * PrPAX(y,3-k) * SUM(PrV)                                               
        endif
         if (cat(3)) then  !B GP and HS of A
          do v=1,3
            PrV(v) = AKA2P(x, v,z) * SUM(AKA2P(v,y,:)*PrPB)*PrPA(y,3-k)
            if (Genos(l,B)/=-9) then
              PrV(v) = PrV(v) * OcA(Genos(l,B), v)
            endif
          enddo
          PrXZ(x,y,z,3) =  PrXZ(x,y,z,3) * SUM(PrV)
        endif
      enddo
    enddo
  enddo
  do z=1,3  ! inbred/non-inbred
    if (cat(z)) then
      PrL(l,z) = LOG10(SUM(PrXZ(:,:,:,z)))
    endif             
  enddo
enddo

LLtmp = SUM(PrL, DIM=1)
WHERE(LLtmp(1:2) <0)  LLtmp(1:2) = LLtmp(1:2) + Lind(B)                                        
 if (Parent(A,k)>0) then
  WHERE(LLtmp <0)  LLtmp = LLtmp - Lind(Parent(A,k))
endif
if (Parent(A,3-k)>0 .and. cat(2)) then
  LLtmp(2) = LLtmp(2) - Lind(Parent(A,3-k))
endif   

LL = MaxLL(LLtmp)
if (LL >= 0) then
  LL = 777
endif

end subroutine PairGP

! #####################################################################

subroutine PairGGP(A, B, kIN, LL)   
! calculates LL that B is maternal(k=1) paternal(k=2), or double(k3) ggp
use Global
implicit none

integer, intent(IN) :: A,B,kIN
double precision, intent(OUT) :: LL
integer :: l, x, y,z,w, AncB(2, mxA), m, MaybeF, k, n, AncA(2,mxA)
double precision :: PrL(nSnp,3), PrXY(3,3), PrXZ(3,3,3,3,3), LLtmp(3),&
  PrPA(3),PrB(3), PrG(3,2), PrPAX(3,2)  
  
LL = 999D0
LLtmp = 999D0
do k=1,2
  if (k/=kIN .and. kIN/=3) cycle
  if (AgeDiff(A, B) <= 0) then  ! B younger than A
    LL = 777
  else if (B==Parent(A,k)) then
    LL = 777
  else
    call GetAncest(B, k, AncB)
    if (ANY(AncB == A)) then
      LL = 777
    else if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
      LL = 444   ! not impossible, just unlikely.
    endif
  endif
  if (LL==777) return

  if (Parent(A,k)>0) then 
    if (ANY(Parent(Parent(A,k), :)/=0) .and. kIN/=3) then
      LL = 444    ! should be picked up elsewere
    else
      call PairGP(Parent(A,k), B, k, 4, LLtmp(1))
      if (LLtmp(1) > 0) then    
        LL = LLtmp(1)
      endif
    endif
  else if (Parent(A,k)<0) then
    if (ANY(GpID(:,-Parent(A,k),k)/=0)) LL = 444
  endif
  if (LL/=999) return
enddo

MaybeF = 0
call getAncest(A, k, AncA)                           
if (kIN==3) then
  MaybeF = 1
  k = 1
else 
  k = kIN
  if (Parent(A,3-k)==0) then
    MaybeF = 1
  else if (Parent(A,k)==0 .and. Complx==2) then
    if (ANY(AncA(:,5-k) == 0)) then   ! GP via parent 3-k
      MaybeF = 1  ! maybe inbreeding loop
    endif
  endif
endif                          

PrL = 0D0    
do l=1,nSnp
  PrXY = 0
  PrXZ = 0
  if (Genos(l,A)==-9) cycle
  call ParProb(l, B, 0, 0, 0, PrB)
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)  ! with GP con
  PrPAX = 1
  do n=1,2
    if (n/=k .and. MaybeF/=1) cycle
    if (Parent(A,n)/=0) then
      call ParProb(l, Parent(A,n), n, -4, 0, PrPAX(:,n))  !=OcA if >0
    endif
    PrG(:,n) = AHWE(:,l)
    do m=1,2
      if (AncA(m, n+2)/=0) then ! either GP ==0
        call ParProb(l,AncA(m, n+2), m,0,0, PrG(:,n))    
      endif
    enddo
  enddo
  
  PrXZ = 0
  do x=1,3  
    do y=1,3 
      PrXY(x,y) = SUM(OKA2P(Genos(l,A),x,:) * PrPA) * PrPAX(x,k) * &
       AKAP(x,y,l) * SUM(AKAP(y, :, l) * PrB)
      if (MaybeF==1) then
        do z=1,3  !consider double GGP (2x k, or k & 3-k)
          do w=1,3
            PrXZ(x,y,z,w,2) = OKA2P(Genos(l,A),x,z) *PrPAX(x,k) *&
             PrPAX(z,3-k) * SUM(AKA2P(z,w,:)* PrG(:,3-k))* &
             SUM(AKA2P(x,y,:)*PrG(:,k)) * &
             SUM(AKAP(y,:,l) * AKAP(w,:,l) * PrB)
          enddo
          PrXZ(x,y,z,y,3) = OKA2P(Genos(l,A),x,z) *PrPAX(x,k) *&
             PrPAX(z,3-k) * SUM(AKA2P(z,y,:)* PrG(:,3-k))* &
             SUM(AKA2P(x,y,:)*PrG(:,k)) * &
             SUM(AKAP(y,:,l) * PrB)
        enddo
      endif
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrXY))
  if (MaybeF==1) then
    PrL(l,2) = LOG10(SUM(PrXZ(:,:,:,:,2)))
    PrL(l,3) = LOG10(SUM(PrXZ(:,:,:,:,3)))
  endif
enddo

LLtmp = SUM(PrL, DIM=1)
do n=1,2
  if (Parent(A,n)>0 .and. (n==k .or. MaybeF==1)) then
    if (n==k) then
      LLtmp = LLtmp - Lind(Parent(A,n))
    else if (MaybeF==1) then
      LLtmp(2:3) = LLtmp(2:3) - Lind(Parent(A,n))
    endif
  endif
enddo

if (kIN==3) then
  LL = MaxLL(LLtmp(2:3))
else if (MaybeF==1) then
  LL = MaxLL(LLtmp)
else
  LL = LLtmp(1)
endif
LL = LL + Lind(B)

end subroutine PairGGP

! #####################################################################

 subroutine PairGA(A, B, k, hf, LL)   ! B FS/HS of GP
use Global
implicit none

integer, intent(IN) :: A,B,k, hf
double precision, intent(OUT) :: LL
integer :: l, x, y,v,w, AncB(2, mxA), m, n, i, BB(maxSibSize), nFSB
double precision :: PrL(nSnp, 2), PrX(3,3,3,3,2), PrPA(3), PrGG(3, 2) 

LL = 999D0
if (AgeDiff(A, B) <= 0) then  ! B younger than A
  LL = 777
else if (ANY(Parent(A,:)==B)) then
  LL = 444
else
  call GetAncest(B, k, AncB)
  if (ANY(AncB == A)) then
    LL = 777
  else if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
    LL = 444   ! not impossible, just unlikely.
  endif
endif
if (LL==777) return

m = 3-k  ! most neutral, doesn't matter in most cases
if (Parent(A,k)/=0) then
  LL = 444
  return
endif

nFSB = nFS(FSID(maxSibSize+1,B))
BB = FSID(1:maxSibSize, FSID(maxSibSize+1,B)) 

PrL = 0D0    
do l=1,nSnp
  PrX = 0D0
  call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA) 
  do n=1,2
    call ParProb(l, Parent(B, n), n, B, -1, PrGG(:,n))
  enddo

  do v=1,3
    do w=1,3
      do x=1,3  
        do y=1,3
          PrX(x,y,v,w,:) = AKAP(x,y,l) * PrGG(v,1) * PrGG(w,2)
          if (hf==3) then
            PrX(x,y,v,w,:) = PrX(x,y,v,w,:) * AKA2P(y, v, w)
          else if (hf==1) then
            PrX(x,y,v,w,:) = PrX(x,y,v,w,:) * AKAP(y, v, l)
          else if (hf==2) then
            PrX(x,y,v,w,:) = PrX(x,y,v,w,:) * AKAP(y, w, l)
          endif
          do i = 1, nFSB
            if (BB(i) == B)  cycle
            if (Genos(l,BB(i))/=-9) then
              PrX(x,y,v,w,:) = PrX(x,y,v,w,:) * OKA2P(Genos(l,BB(i)), v, w)
            endif  
          enddo
          if (Genos(l,A)/=-9) then
            PrX(x,y,v,w,2) = PrX(x,y,v,w,2) * SUM(OKA2P(Genos(l,A),x,:) * PrPA)
          endif          
        enddo
      enddo
      if (Genos(l,B)/=-9) then
        PrX(:,:,v,w,2) = PrX(:,:,v,w,2) * OKA2P(Genos(l,B), v, w)
      endif
    enddo
  enddo
  do n=1,2
    PrL(l,n) = LOG10(SUM(PrX(:,:,:,:,n)))
  enddo
enddo

LL = SUM(PrL(:,2)) - SUM(PrL(:,1))

end subroutine PairGA

! #####################################################################

subroutine PairUA(A, B, kA, kB, LL)
! B half sib or full sib (kB=3) of parent kA of A?
use Global
use qsort_c_module
implicit none

integer, intent(IN) :: A,B,kA, kB  ! kB=3 : full sibs
double precision, intent(OUT) :: LL
integer :: l, x, g, y, z, GG(2), GA(2), PB(2), PA, i, nA, r,u,j,e,Ei,m,&
  AncA(2,mxA), AncG(2, 2,mxA), AA(maxSibSize), cat(maxSibSize), AncB(2,mxA), &
  doneB(maxSibSize), BB(maxSibSize), nB, catG(2), GGP, catB(maxSibSize), &
  nBx(2), BBx(maxSibSize, 2), Bj, Mates(maxSibSize, 2), w, BBf(maxSibSize), &
  nBf, AB(2*maxSibSize)
double precision :: PrL(nSnp), PrG(3,2), PrXYZ(3,3,3), PrPA(3), PrA(3),&
  PrPB(3), PrGA(3), PrAB(3,3,3,2), PrE(3), PrH(3), PrGG(3), &
  PrLX(nSnp, 2), PrEW(3,3), PrW(3), PrY(3)
integer, allocatable, dimension(:) :: UseEE, MateABpar, TypeEE
double precision, allocatable, dimension(:,:,:) :: PrEE
logical :: MateLoop(maxSibSize,2), SIMPL

AA = 0
BB = 0  
if (A>0) then  
  nA = 1
  AA(1) = A
  PA = Parent(A,kA)
  if (PA<0) then
    nA = ns(-PA,kA)
    AA(1:nA) = SibID(1:nA, -PA, kA)
    GA = GpID(:,-PA,kA)
  else if (PA>0) then
    GA = Parent(PA, :)
  else
    GA = 0
  endif
else
  nA = nS(-A, kA)
  AA(1:nA) = SibID(1:nA, -A, kA)
  PA = A
  GA = GpID(:, -A, kA)
endif
  
Mates = 0
nBx = 0
BBx = 0
nBf = 0
BBf = 0
if (B > 0) then
  nB = 1
  BB(1) = B
  PB = Parent(B,:)
  do m=1,2
    if (kB/=3 .and. m/=kB)  cycle
    if (Parent(B, m) >=0) then
      nBx(m) = 1
      BBx(1, m) = B
    else 
      nBx(m) = nS(-Parent(B, m), m)
      BBx(1:nBx(m), m) = SibID(1:nBx(m), -Parent(B, m), m)
    endif
    do j=1,nBx(m)
      Mates(j,m) = Parent(BBx(j, m), 3-m)
    enddo
  enddo
  nBF = nFS(FSID(maxSibsize+1,B))
  BBf = FSID(1:maxSibsize, FSID(maxSibsize+1,B))
else if (B < 0) then
  nB = nS(-B, kB)
  BB(1:nB) = SibID(1:nB, -B, kB)
  PB(kB) = B
  PB(3-kB) = 0
  do j=1,nB
    Mates(j,kB) = Parent(BB(j), 3-kB)
  enddo
endif

LL = 999D0
if (kB < 3) then
  if (B < 0 .and. GA(kB) == B) then
    LL = 888
  else if (B > 0 .and. GA(kB) == PB(kB) .and. PB(kB)/=0) then
    LL = 888
  else if (GA(3-kB)==PB(3-kB) .and. GA(3-kB)/=0) then ! B>0; FA not HA
    LL = 777
  endif
else if (GA(1)==PB(1) .and. GA(2)==PB(2) .and. GA(1)/=0 .and. &
  GA(2)/=0) then  ! kB==3
  LL = 888
endif
if (LL /= 999D0) return

GG = 0  ! parent of B, GP of A
AncG = 0
 call GetAncest(A, kA, AncA)
do x=1,2
  if (x/=kB .and. kB/=3) cycle
  if (GA(x)==0) then
    GG(x) = PB(x)
  else if (GA(x)/=PB(x) .and. PB(x)/=0) then
    LL = 777
  else
    GG(x) = GA(x)
  endif
  if (ANY(AA(1:nA)==GG(x))) then
    LL = 777
  endif
  do j=1,nBx(x)
    if (PA<0 .and. Parent(BBx(j,x), 3-x)<0) then
      if (GpID(kA, -Parent(BBx(j,x), 3-x), 3-x) == PA) then
        LL = 444
      endif
    endif
  enddo
enddo
if (LL /= 999D0) return
do x=1,2
  call GetAncest(GG(x), x, AncG(x, :, :))
enddo

if (A > 0) then
  if (ANY(AncG == A)) then
    LL = 777
  endif
else if (A < 0) then
  if (ANY(AncG(:, kA, 2:mxA) == A)) then
    LL = 777
  endif
endif
if (B > 0) then
  if (ANY(AncG == B)) then
    LL = 777
  endif
else if (B < 0) then
  if (ANY(AncG(:, kB, 3:mxA) == B)) then
    LL = 777
  endif
endif
if (kB<3) then
  if (B<0) then
    if (ANY(AncA(kB,3:mxA)==B))  LL = 444  ! B is GGP; 
  else if (B>0) then
    if (ANY(AncA(:,3:mxA)==B))  LL = 444
  endif
endif
if (LL /= 999D0) return

do x=2,mxA
  do y=1,2
    do g=1,2
      if (AncG(g,y,x) > 0) then
        if (A > 0) then
          if (AgeDiff(A, AncG(g,y,x)) < 0) then
            LL = 777  ! A older than putative ancestor
          endif 
        else if (A<0) then
          if (ANY(AgeDiff(SibID(1:nS(-A,kA),-A,kA),AncG(g,y,x))<0)) then
            LL = 777  ! A older than putative ancestor
          endif
        endif
        if (x==2) cycle 
        if (B > 0) then
          if (AgeDiff(B, AncG(g,y,x)) < 0) then
            LL = 777  
          endif 
        else if (B<0) then
          if (ANY(AgeDiff(SibID(1:nS(-B,kB),-B,kB),AncG(g,y,x))<0)) then
            LL = 777 
          endif
        endif
      endif
    enddo
  enddo
enddo
if (LL /= 999D0) return
!==============================================
 PrL = 0D0
 
if (A>0 .and.  B>0) then  ! quicker.
  if (ALL(Parent(A,:)>=0) .and. ALL(Parent(B,:)>=0)) then
    do l=1, nSnp
      if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
      call ParProb(l, Parent(A,3-kA), 3-kA, A, 0, PrPA)
      if (kB == 3) then
        do g=1,2
          call ParProb(l, GG(g), g, 0, 0, PrG(:,g))  ! >=0
        enddo        
      else
        call ParProb(l, GG(kB), kB, 0, 0, PrG(:,kB))  ! >=0
        if (PA>0) then
          call ParProb(l, GA(3-kB), 3-kB, PA, 0, PrGA)
        else
          call ParProb(l, GA(3-kB), 3-kB, 0, 0, PrGA)
        endif
        call ParProb(l, PB(3-kB), 3-kB, B, 0, PrPB)
      endif
      
      PrXYZ = 0
      do z=1,3
        do y=1,3
          do x=1,3
            if (kB == 3) then
              PrXYZ(x,y,z) =AKA2P(x,y,z)*PrG(y,1)*PrG(z,2)
            else
              PrXYZ(x,y,z) =AKA2P(x,y,z)*PrG(y,kB)*PrGA(z)   
            endif
            if (Parent(A,3-kA)/=B .or. kB/=3) then
              if (Genos(l,A)/=-9) then
                PrXYZ(x,y,z) =PrXYZ(x,y,z) *SUM(OKA2P(Genos(l,A), x, :)&
                 * PrPA)  
              endif
              if (Genos(l,B)/=-9) then
                if (kB==3) then
                  PrXYZ(x,y,z) = PrXYZ(x,y,z) *OKA2P(Genos(l,B), y, z)   
                else
                  PrXYZ(x,y,z) =PrXYZ(x,y,z) *SUM(OKA2P(Genos(l,B),y,:)&
                   * PrPB)
                endif
              endif
            else  ! FS mating
              do u=1,3
                PrH(u) = AKA2P(u,y,z)
                if (Genos(l,A)/=-9)  PrH(u)=PrH(u)*OKA2P(Genos(l,A),x,u) 
                if (Genos(l,B)/=-9)  PrH(u) = PrH(u) *OcA(Genos(l,B), u)
              enddo
              PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(PrH)   
            endif
            if (Parent(A,kA)>0) then
              if (Genos(l,Parent(A,kA))/=-9) then
                PrXYZ(x,y,z) =PrXYZ(x,y,z) *OcA(Genos(l,Parent(A,kA)),x)
              endif
            endif
          enddo
        enddo
      enddo
      PrL(l) = LOG10(SUM(PrXYZ))
    enddo
    
    LL = SUM(PrL)
    return
  endif
endif

!==============================================

if (B<0 .and. A>0) then
  SIMPL = .TRUE.
  if(ANY(Parent(A,:)/=0)) then
    SIMPL = .FALSE.
  endif
  do j=1,nB
    call getAncest(BB(j), kA, AncB)
    if (ANY(AncB == A)) SIMPL = .FALSE.
  enddo
  if (SIMPL) then
    do l=1, nSnp
      do y=1,3
        PrY(y) = XPr(3,y,l, -B,kB)
        if (Genos(l,A)==-9) cycle
        PrY(y) = PrY(y)*SUM(OKAP(Genos(l,A),:,l) * AKAP(:,y,l))
      enddo
      PrL(l) = LOG10(SUM(PrY))
    enddo
    
    LL = SUM(PrL)
    return
  endif
endif 

!==============================================

allocate(UseEE(nA+nB))
allocate(TypeEE(nA+nB))
allocate(PrEE(3,3, nA+nB))
allocate(MateABpar(nA+nB))
UseEE = 0
AB = 0

if (kA==kB) then
  AB(1:nA) = AA(1:nA)
  AB((nA+1):(nA+nB)) = BB(1:nB)
  call FindEE(AB(1:(nA+nB)), nA, nB, kA, UseEE, MateABpar)  ! may reorder AA, BB
  AA(1:nA) = AB(1:nA)
  BB(1:nB) = AB((nA+1):(nA+nB))
  TypeEE = 3-kA
else if (kA/=kB  .and. kB/=3) then
  call FindEE(AA(1:nA), nA, 0, kA, UseEE(1:nA), MateABpar(1:nA)) 
  call FindEE(BB(1:nB), nB, 0, kB, UseEE((nA+1):(nA+nB)), MateABpar((nA+1):(nA+nB)))
  do i=1, nB
    if (UseEE(nA+i)/=0) then
      UseEE(nA+i) = nA + UseEE(nA+i)
    endif
  enddo
  TypeEE(1:nA) = 3-kA
  TypeEE((nA+1):(nA+nB)) = 3-kB
endif                          

!============================================

 cat=0
 catG = 0
 catB = 0
GGP = 0
do i = 1, nA
  if (Parent(AA(i),3-kA)==0) cycle
  if (kA/=kB .and. GG(3-kA)<0 .and. &
    Parent(AA(i), 3-kA) == GG(3-kA)) then  !incl. kB=3
    cat(i) = 1  
  else if (kA==kB .and. Parent(AA(i), 3-kA)==GA(3-kA) .and. GA(3-kA)<0) then
    cat(i) = 2
    UseEE(i) = 0
  else 
    if (Parent(AA(i), 3-kA)<0) then                               
      if (GpID(kA,-Parent(AA(i), 3-kA),3-kA) == PA .and. PA/=0) then
        cat(i) = 7  ! Ai inbred
      endif
    endif 
    do j=1, nB
      if (AA(i) == BB(j) .or. kA/=kB) cycle
      if (Parent(AA(i), 3-kA) == Parent(BB(j), 3-kA)) then
        cat(i) = 3
      else if (Parent(AA(i), 3-kA) == BB(j)) then
        cat(i) = -j
      endif
    enddo
  endif
  do g=1,2
    if (kB/=g .and. kB/=3) cycle
    if (Parent(AA(i), 3-kA) < 0) then
      if (GpID(g,-Parent(AA(i), 3-kA),3-kA) == GG(g) .and. GG(g)/=0) then
        if (g==kB .or. (kB==3 .and. g==3-kA)) then
          if (ALL(GpID(:,-Parent(AA(i), 3-kA),3-kA) == GG) .and. ALL(GG/=0)) then
            cat(i) = 10
          else
            cat(i) = 8  ! via y
          endif
        else
          cat(i) = 9  ! via z
        endif
      endif
    endif
    if (GG(g) > 0) then
      if (Parent(AA(i),3-kA) == Parent(GG(g), 3-kA)) then
        cat(i) = 5  ! TODO? 4+g when kB==3
        catG(g) = 2
        GGP = Parent(GG(g), kA)
        UseEE(i) = 0
      endif
    else if (GG(g) < 0) then
      if (Parent(AA(i),3-kA) == GpID(3-kA, -GG(g),g)) then
        cat(i) = 5  
        catG(g) = 2
        GGP = GpID(kA, -GG(g),g)
        UseEE(i) = 0
      endif
    endif
  enddo
enddo
if (kB/=3) then   ! TODO: for kB==3
  do j=1, nB
    if (Parent(BB(j),3-kB)==0) cycle
    if (Parent(BB(j), 3-kB) == GA(3-kB) .and. GA(3-kB)<0) then
      catB(j) = 2
      UseEE(nA+j) = 0
    else if (Parent(BB(j),3-kB)<0) then
      if (GpID(kB, -Parent(BB(j),3-kB),3-kB) == GG(kB) .and. GG(kB)/=0) then
        catB(j) = 7
      else if (GpID(kA, -Parent(BB(j),3-kB),3-kB) == PA .and. PA/=0) then
        catB(j) = 8
      endif
    endif
    do g=1,2
      if (GG(g) > 0) then
        if (Parent(BB(j),3-kB) == Parent(GG(g), 3-kB)) then
          if(catG(g)==0) then
            catB(j) = 5  
            catG(g) = 3
            GGP = Parent(GG(g), kB)
            UseEE(nA+j) = 0  ! ??
          endif
        endif
      else if (GG(g) < 0) then
        if (Parent(BB(j),3-kB) == GpID(3-kB, -GG(g),g)) then           
           if(catG(g)==0) then
            catB(j) = 5          
            catG(g) = 3
            GGP = GpID(kB, -GG(g),g)
            UseEE(nA+j) = 0
          endif
        endif
      endif   
    enddo
    if (ANY(cat == 8) .and. kB/=3 .and. catB(j)==0) then
      do i=1,nA
        if (PA<0 .and. NFS(AA(i))==0) cycle
        if (Parent(AA(i), 3-kA)>=0) cycle
        if (GpID(3-kB,-Parent(AA(i), 3-kA),3-kA) == Parent(BB(j),3-kB)) then
          catB(j) = -i
        endif
      enddo
    endif    
  enddo
endif

if (kB/=3) then
  if (GG(kB) > 0) then
    if (Parent(GG(kB), 3-kB) == GA(3-kB) .and. GA(3-kB)/=0) then
      catG(kB) = 1
      GGP = Parent(GG(kB), kB)
    endif
  else if (GG(kB) < 0) then
    if (GpID(3-kB, -GG(kB), kB) == GA(3-kB) .and. GA(3-kB)/=0) then
      catG(kB) = 1
      GGP = GpID(kB, -GG(kB), kB)
    endif
  endif
endif

MateLoop = .FALSE.
if (B>0) then  
  do m=1,2
    if (kB/=3 .and. m/=kB)  cycle
    do j=1, nBx(m)
      Bj = BBx(j, m)
      if (nFS(Bj)==0) cycle  !  .and. Bj/=B
      if (kB==3 .and. Parent(Bj,1)==GG(1) .and.  Parent(Bj,2)==GG(2))  cycle
      if (Parent(Bj,m)<0 .and. Parent(Bj,3-m)<0) then
        do g=1, nS(-Parent(Bj, 3-m),3-m)
          Ei = SibID(g,-Parent(Bj,3-m),3-m)
          if (Parent(Ei,m)>=0 .or. Parent(Ei,m)==Parent(Bj,m)) cycle
          if (ANY(Mates(:,3-m) == Parent(Ei, m))) then
            MateLoop(j,m) = .TRUE.
          endif
        enddo
      endif
    enddo
  enddo
endif

!==============================================

PrL = 0D0
PrLx = 0D0
DoneB = 0
SIMPL = ALL(cat==0) .and. ALL(catG==0) .and. ALL(catB==0) .and. &
  ALL(GG >=0) .and. ALL(UseEE==0) .and. .not. ANY(MateLoop)
 !   ( .or. (B>0 .and. kB<3 .and. .not. ANY(MateLoop)))
 
do l=1,nSnp
  if (A>0 .and. B>0) then
    if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
  endif
  do g=1,2
    if (g/=kB .and. kB/=3) cycle
    if (catG(g)==0) then
      if (SIMPL .and. B>0) then
        call ParProb(l, GG(g), g, B, -1, PrG(:,g)) 
      else
        call ParProb(l, GG(g), g, -1, 0, PrG(:,g))
      endif
    else
      if (GG(g) > 0) then 
        call ParProb(l, GG(g), g, 0, 0, PrG(:,g)) 
        if (catG(g)==1) then
          call ParProb(l, GGP, kB, GG(g), 0, PrGG)  !ALL(GG >=0)
        else if (catG(g)==2) then
          if (Genos(l,GG(g))/=-9) then
            PrG(:,g) = OcA(Genos(l,GG(g)),:)
          endif
          call ParProb(l, GGP, kA, GG(g), 0, PrGG)
        else if (catG(g)==3) then
          call ParProb(l, GGP, kB, GG(g), 0, PrGG)
        endif
      else if (GG(g) < 0) then
        PrG(:,g) = 1
        if (catG(g)==1) then
          call ParProb(l, GGP, kB, 0, 0, PrGG)
        else if (catG(g)==2) then
          call ParProb(l, GGP, kA, 0, 0, PrGG)
        else if (catG(g)==3) then
          call ParProb(l, GGP, kB, 0, 0, PrGG)
        endif
      else 
          PrG(:,g) = AHWE(:,l)
      endif
    endif
  enddo
  if (kB/=3) then 
    if (ANY(cat==2)  .or. ANY(catB==2)) then
      call ParProb(l, GA(3-kB), 3-kB, -1, 0, PrGA)
    else if (PA>0) then
      call ParProb(l, GA(3-kB), 3-kB, PA, 0, PrGA)
    else if (catG(kB)==1 .and. GG(kB)>0) then
      call ParProb(l, GA(3-kB), 3-kB, GG(kB), 0, PrGA)
    else
      call ParProb(l, GA(3-kB), 3-kB, 0, 0, PrGA)
    endif
    if (B>0) then
      call ParProb(l, PB(3-kB), 3-kB, B, 0, PrPB)   ! -1?
    endif
  endif

  ! === 
  
  PrA = 1
  if (SIMPL) then
    if (A < 0) then
      PrA = Xpr(1,:,l, -A,kA)
    else if (A>0) then
      call ParProb(l, Parent(A,3-kA), 3-kA, A, 0, PrPA)
      if (Genos(l,A)/=-9) then
        do x=1,3 
          PrA(x) = SUM(OKA2P(Genos(l,A), x, :) * PrPA)
        enddo
      else
        PrA = 1 ! ? AHWE(:,l)  
      endif
      if (Parent(A,kA)>0) then
        if (Genos(l,Parent(A,kA))/=-9) then
          PrA = PrA * OcA(Genos(l,Parent(A,kA)), :)  
        endif
      else if (Parent(A,kA)<0 .and. Genos(l,A)/=-9) then   
        do x=1,3
          PrH(x) = Xpr(1,x,l, -Parent(A,kA),kA) /&
            SUM(OKA2P(Genos(l,A),x,:)*PrPA)
        enddo
        PrH = PrH / SUM(PrH) 
        do x=1,3 
          PrA(x) = SUM(OKA2P(Genos(l,A),x,:) * PrPA * PrH(x))
        enddo
      endif
    endif
  
    do x=1,3  ! PA, kA
      do y=1,3  ! PrG, kB
        do z=1,3  ! PrGA, 3-kB / PrG, 3-kB
          if (kB==3 .and. B>0) then  ! SA/PA FS of B; 
            PrXYZ(x,y,z) = PrA(x) * AKA2P(x,y,z)*PrG(y,3-kA) * PrG(z,kA)    
            if (Genos(l,B)/=-9) then
              PrXYZ(x,y,z) = PrXYZ(x,y,z)*OKA2P(Genos(l,B), y, z)
            endif
          else 
            PrXYZ(x,y,z) = PrA(x) * AKA2P(x,y,z) * PrGA(z)
            if (B>0) then
              PrXYZ(x,y,z) = PrXYZ(x,y,z) * PrG(y, kB)
              if (Genos(l,B)/=-9) then
                PrXYZ(x,y,z) = PrXYZ(x,y,z)*&
                 SUM(OKA2P(Genos(l,B),y,:) * PrPB)
              endif
            else if (B<0) then
              PrXYZ(x,y,z) =PrXYZ(x,y,z) *XPr(3,y,l,-B,kB)
            endif
          endif     
        enddo
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))   

  else 

    PrAB = 0   
    do y=1,3  ! PrG, kB
      PrEE = 0
      do x=1,3  ! PA, kA
        do z=1,3
          if (kB==3) then
            PrAB(x,y,z,:) = AKA2P(x,y,z) *PrG(z,kA) *PrG(y,3-kA)
          else if (catG(kB)==1) then
            PrAB(x,y,z,:) = AKA2P(x,y,z)*PrGA(z)*&
              SUM(AKA2P(y,z,:) * PrGG) * PrG(y, kB)  
          else if (catG(kB)==2) then
            PrAB(x,y,z,:) = AKA2P(x,y,z)*PrGA(z)
          else
            PrAB(x,y,z,:) = AKA2P(x,y,z) *PrGA(z) * PrG(y,kB) 
          endif
        enddo
        if (PA>0) then
          if (Genos(l,PA)==-9) cycle
          PrAB(x,y,:,:) = PrAB(x,y,:,:) * OcA(Genos(l,PA), x) 
        endif
      enddo   
        
      do x=1,3
       doneB = 0
       do z=1,3       
        do r=1, nA
          if (PA<0 .and. NFS(AA(r))==0) cycle
          if (cat(r)>2 .and. cat(r)<7) then
            call ParProb(l,Parent(AA(r),3-kA),3-kA,-1,0,PrE)
            PrEW = 0
          else if (cat(r)==7) then
            call ParProb(l, GpID(3-kA,-Parent(AA(r),3-kA),3-kA), 3-kA, 0,0,PrH)
            do e=1,3
              PrE(e) = SUM(AKA2P(e,x,:) * PrH)
            enddo
          else if (cat(r)>=8 .and. cat(r)<=10) then
            if (kB < 3) then
              if (cat(r)==8) then
                g=kB
              else if (cat(r)==9) then
                g=3-kB
              endif
            else 
              if (cat(r)==8) then
                g=3-kA
              else if (cat(r)==9) then
                g=kA
              endif
            endif
            if (cat(r) < 10) then
              call ParProb(l, GpID(3-g,-Parent(AA(r),3-kA),3-kA), 3-g, 0,0,PrH)
            endif
            do e=1,3
              if (cat(r)==8) then
                PrE(e) = SUM(AKA2P(e,y,:) * PrH)
              else if (cat(r)==9) then
                PrE(e) = SUM(AKA2P(e,:,z) * PrH)
              else if (cat(r)==10) then
                PrE(e) = AKA2P(e,y,z)
              endif
            enddo  
            PrE = PrE/SUM(PrE)        
          else if (cat(r)==42) then
            cycle ! do with B  (catB(j) = -i)                               
          else if (UseEE(r)/=0) then
            call ParProb(l, MateABpar(r), 3-TypeEE(r), 0,0,PrH)
            do e=1,3
              do u=1, 3
                PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,x,UseEE(r)) * PrH)
              enddo
              PrE(e) = SUM(PrW)
            enddo
            PrE = PrE/SUM(PrE)
          else if (cat(r) < 0) then
            if (kB<3) then
              call ParProb(l, Parent(BB(-cat(r)),3-kB), 3-kB, BB(-cat(r)),0,PrH)
            else
              PrH = PrG(:,kA)
            endif
            do e=1,3
              PrE(e) = SUM(AKA2P(e,y,:) * PrH)
            enddo
            if (Genos(l,BB(-cat(r)))/=-9) PrE = PrE * OcA(Genos(l,BB(-cat(r))), :)
            PrE = PrE/SUM(PrE)
          else if (cat(r)==0) then
            call ParProb(l,Parent(AA(r),3-kA),3-kA,-1,0,PrE)
          else
            PrE = 1
          endif

          if (Parent(AA(r), 3-kA) < 0 .and. Parent(AA(r), 3-kA)/=GG(3-kA)) then
            if (A>0) then
              do i=1, MAX(nFS(AA(r)),1)
                if (FSID(i, AA(r))==A) cycle
                if (ANY(GG == FSID(i, AA(r)))) cycle
                if (Genos(l,FSID(i, AA(r)))==-9) cycle
                 PrE=PrE*OKA2P(Genos(l,FSID(i,AA(r))),x,:)  ! FS of A
              enddo
            endif
            
            do e=1,3
              if (cat(r)==1 .and. e/=y)  cycle
              if (cat(r)==2 .and. e/=z)  cycle                                                             
              do g=1, nS(-Parent(AA(r), 3-kA), 3-kA)
                Ei = SibID(g, -Parent(AA(r), 3-kA),3-kA)
                if (nFS(Ei) == 0) cycle 
                if (Parent(Ei,kA)==PA .and. PA/=0) cycle
                if (kA==kB .and. Parent(Ei,kA)==GG(kA) .and. GG(kA)/=0)  cycle                                              
                if (kB<3) then
                  if (Parent(Ei, kB)== PB(kB) .and. PB(kB)/=0) cycle
                endif                
                if (cat(r)==5 .and. Parent(Ei,kA)==GGP .and. GGP/=0) then
                  PrH = PrGG
                  do i=1, nFS(Ei)
                    if (GG(kA) == FSID(i, Ei)) cycle
                    if (Genos(l,FSID(i, Ei))==-9) cycle
                    PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)),:,e)
                  enddo
                  if (catG(kA)==2 .and. kB==3) then 
                      PrEW(e,z) = PrEW(e,z) * SUM(AKA2P(z,e,:) * PrH)  
                  else if (ANY(catG==2)) then
                    PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrH)
                  endif                
                else
                  call ParProb(l,Parent(Ei,kA),kA,Ei,-1,PrH) 
                  do i=1, nFS(Ei)
                    if (A>0 .and. FSID(i, Ei)==A) cycle
                    if (B>0 .and. FSID(i, Ei)==B) cycle
                    if (Genos(l,FSID(i, Ei))==-9) cycle
                    PrH=PrH*OKA2P(Genos(l,FSID(i,Ei)),:,e)
                  enddo
                  if (SUM(PrH)<3) then
                    PrE(e) = PrE(e) * SUM(PrH)
                  endif
                endif
              enddo  ! g
            enddo  ! e
            if (cat(r)==3 .and. B>0) then   ! TODO: nBx?
              do j=1,nB
                if (Parent(BB(j), 3-kA) /= Parent(AA(r), 3-kA)) cycle
                do i=1, MAX(nFS(BB(j)),1)
                  if (FSID(i,BB(j))==B) cycle
                  if (Genos(l,FSID(i,BB(j)))==-9) cycle
                  PrE = PrE * OKA2P(Genos(l,FSID(i,BB(j))), y, :)
                enddo
              enddo
            endif   
          endif
          if (cat(r)==5 .and. (Parent(AA(r), 3-kA)>0 .or. GGP==0)) then
            do e=1,3
              PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)
            enddo
          endif
          
          if (cat(r)==5) then
            if (catG(kA)==2 .and. SUM(PrEW)>0) then      
              PrAB(x,y,z,1)=PrAB(x,y,z,1)*SUM(PrEW(:,z))  
            else if (SUM(PrE)<3) then   
              PrAB(x,y,z,1) = PrAB(x,y,z,1) * SUM(PrE)
            endif
          else if (cat(r)==1) then 
            PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(y)
          else if (cat(r)==2) then 
            PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(z)
          else if (SUM(PrE)<3) then   
            PrAB(x,y,z,1)=PrAB(x,y,z,1)*SUM(PrE)
          endif       

          do i=1, MAX(nFS(AA(r)),1)
            if (A>0 .and. FSID(i, AA(r))/=A .and. .not. &
             ANY(BB==FSID(i, AA(r)))) cycle
            if (Genos(l,FSID(i, AA(r)))/=-9) then
              PrE = PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)  ! <- A
              if (SUM(PrEW)>0) then
                do e=1,3
                  PrEW(e,:) = PrEW(e,:) * OKA2P(Genos(l,FSID(i,AA(r))), x, e)
                enddo
              endif
            endif
          enddo

          if (cat(r)==3 .or. (cat(r)==5 .and. ANY(catB==5)) .or. &
            (cat(r)==2 .and. ANY(catB==2))) then 
            do j=1,nB
              if (Parent(BB(j),3-kA) /= Parent(AA(r),3-kA)) cycle
              if (ANY(AA == BB(j)))  cycle
              if (B>0 .and. BB(j)/=B) cycle
                if (Genos(l,BB(j))==-9) cycle
                PrE = PrE * OKA2P(Genos(l,BB(j)), y, :)
              DoneB(j) = 1
            enddo
          endif
          
          if (cat(r)==1) then
            PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(y)
          else if (cat(r)==2) then
            PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(z)
          else if (cat(r)==5) then
            if (catG(kA)==2 .and. SUM(PrEW)>0) then    
              PrAB(x,y,z,2)=PrAB(x,y,z,2)*SUM(PrEW(:,z))  
            else if (SUM(PrE)<3) then   
              PrAB(x,y,z,2) = PrAB(x,y,z,2) * SUM(PrE)
            endif
          else if (SUM(PrE)<3) then
            PrAB(x,y,z,2)=PrAB(x,y,z,2)*SUM(PrE)
          endif
          PrEE(:,x,r) = PrE
        enddo  ! r
       enddo  ! z
      enddo  ! x  
      
      do x=1,3
        if (x>1 .and. ALL(UseEE==0) .and. all(catB>=0))  cycle
      if (B<0) then            
        do j=1,nB
          if (nFS(BB(j))==0) cycle
          if (DoneB(j)==1) cycle
          if (kA/=kB .and. PA<0 .and. Parent(BB(j), 3-kB)==PA) cycle
          DoneB(j) = 2  ! for output check only
          if (catB(j)==2) then
            PrE = 1
          else if (catB(j)==7) then
            call ParProb(l, GpID(3-kB,-Parent(BB(j),3-kB),3-kB), 3-kB, 0,0,PrH)
            do e=1,3
              PrE(e) = SUM(AKA2P(e,y,:) * PrH)
            enddo
          else if (catB(j)==8) then
            call ParProb(l, GpID(3-kA,-Parent(BB(j),3-kB),3-kB), 3-kA, 0,0,PrH)
            do e=1,3
              PrE(e) = SUM(AKA2P(e,x,:) * PrH)
            enddo
          else if (UseEE(nA+j)/=0) then
            call ParProb(l, MateABpar(nA+j), 3-TypeEE(nA+j), 0,0,PrH)
            do e=1,3
              do u=1, 3
                PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,x,UseEE(nA+j)) * PrH)
              enddo
              PrE(e) = SUM(PrW)
            enddo
            PrE = PrE/SUM(PrE)
          else
            call ParProb(l, Parent(BB(j), 3-kB), 3-kB, -1,0,PrE)
          endif

          if (Parent(BB(j), 3-kB) < 0) then  
            do e=1,3
              do g=1, nS(-Parent(BB(j), 3-kB), 3-kB)
                Ei = SibID(g, -Parent(BB(j), 3-kB),3-kB)
                if (nFS(Ei) == 0) cycle
                if (Parent(Ei, kB) == PB(kB) .and. PB(kB)/=0) cycle  
                if (Parent(Ei, kA)== PA .and. PA/=0) cycle
                call ParProb(l,Parent(Ei,kB),kB,Ei,-1,PrH) 
                do i=1, nFS(Ei)
                  if (A>0 .and. FSID(i, Ei)==A) cycle
                  if (Genos(l,FSID(i, Ei))/=-9) then
                    PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                  endif
                enddo
                PrE(e) = PrE(e) * SUM(PrH) 
              enddo
            enddo                   
          endif
          
          if (ANY(UseEE/=0) .or. ANY(catB<0)) then
            if (catB(j)==2) then
              PrAB(x,y,:,1) = PrAB(x,y,:,1) * PrE
            else if (SUM(PrE)<3) then
              PrAB(x,y,:,1) = PrAB(x,y,:,1) * SUM(PrE)
            endif
          else if (catB(j)==2) then
            do z=1,3
              PrAB(:,y,z,1) = PrAB(:,y,z,1) * PrE(z)
            enddo
          else if (SUM(PrE)<3) then
            PrAB(:,y,:,1) = PrAB(:,y,:,1) * SUM(PrE)
          endif

          if (catB(j)==5) then
            do e=1,3
              PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)
            enddo
          endif                    
          do u=1, nFS(BB(j))
            if (Genos(l,FSID(u, BB(j)))==-9) cycle
            PrE = PrE * OKA2P(Genos(l,FSID(u,BB(j))), y, :)
          enddo
          
          if (ANY(UseEE/=0) .or. ANY(catB<0)) then
            if (catB(j)==2) then
              PrAB(x,y,:,2) = PrAB(x,y,:,2) * PrE
            else  if (SUM(PrE)<3) then
              PrAB(x,y,:,2) = PrAB(x,y,:,2) * SUM(PrE)
            endif
            PrEE(:,x,nA+j) = PrE
          else if (catB(j)==2) then 
            do z=1,3
              PrAB(:,y,z,2) = PrAB(:,y,z,2) * PrE(z)
            enddo
          else if (SUM(PrE)<3) then
            PrAB(:,y,:,2) = PrAB(:,y,:,2) * SUM(PrE) 
          endif
        enddo   ! j

      else if (B>0) then 
        do z=1,3  
          do m=1,2
            if (kB/=3 .and. m/=kB)  cycle
            do j=1, nBx(m)
              Bj = BBx(j, m)
              if (nFS(Bj)==0 .and. Bj/=B) cycle    
              if (B>0 .and. ANY(FSID(:,Bj)==B) .and. DoneB(1)==1)  cycle                                                          
              if (kA/=kB .and. PA<0 .and. Parent(Bj, kA)==PA) cycle
              if (kB==3 .and. Parent(Bj, 3-m)==GG(3-m) .and. GG(3-m)/=0) then  ! FS of B
                if (Parent(Bj,1)<0 .and. Parent(Bj,2)<0 .and. m==2) cycle
                do u=1, nFS(Bj)
                  if (FSID(u,Bj)==B) cycle
                  if (ANY(AA(1:nA)==FSID(u,Bj))) cycle
                  if (Genos(l,FSID(u, Bj))==-9) cycle
                  if (ALL(UseEE==0)) then
                    PrAB(:,y,z,:) = PrAB(:,y,z,:) * OKA2P(Genos(l,FSID(u,Bj)),y,z) 
                  else
                    PrAB(x,y,z,:) = PrAB(x,y,z,:) * OKA2P(Genos(l,FSID(u,Bj)),y,z)
                  endif  
                enddo
              else if (kB==3 .and. Parent(Bj,1)<0 .and. Parent(Bj,2)<0 &
               .and. MateLoop(j,m)) then  
                if (m==2) cycle
                call ParProb(l,Parent(Bj,m),m,-1,0,PrE)
                call ParProb(l,Parent(Bj,3-m),3-m,-1,0,PrW)
                
                do g=1, nS(-Parent(Bj, 3-m),3-m)
                  Ei = SibID(g,-Parent(Bj,3-m),3-m)
                  if (nFS(Ei) == 0) cycle
                  do e=1,3
                    do w=1,3
                      PrEW(e,w) = PrE(e) * PrW(w)
                      if (Parent(Ei,m)==Parent(Bj,m) .and. &
                       Parent(Ei,3-m)==Parent(Bj,3-m)) then
                        do i=1, nFS(Ei)
                          if (FSID(i, Ei)==B) cycle
                          if (Genos(l,FSID(i, Ei))==-9) cycle
                          PrEW(e,w) = PrEW(e,w) * OKA2P(Genos(l,FSID(i,Ei)), e, w)
                        enddo
                      else if (Parent(Ei,m)==Parent(Bj,m)) then
                        call ParProb(l, Parent(Ei,3-m),3-m, Ei, -1, PrH)
                        do i=1, nFS(Ei)
                          if (FSID(i, Ei)==B) cycle
                          if (Genos(l,FSID(i, Ei))==-9) cycle
                          PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                        enddo
                        PrEW(e,w) = PrEW(e,w) * SUM(PrH)
                      else if (Parent(Ei,3-m)==Parent(Bj,3-m)) then
                        call ParProb(l, Parent(Ei,m),m, Ei, -1, PrH)
                        do i=1, nFS(Ei)
                          if (FSID(i, Ei)==B) cycle
                          if (Genos(l,FSID(i, Ei))==-9) cycle
                          PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, w)
                        enddo
                        PrEW(e,w) = PrEW(e,w) * SUM(PrH)
                      endif
                    enddo  ! w
                  enddo  ! e
                enddo  ! sib g
                if (ALL(UseEE==0)) then
                  PrAB(:,y,z,:) = PrAB(:,y,z,:) * SUM(PrEW)
                else
                  do e=1,3
                    PrEE(:,x,nA+1) = SUM(PrEW(e,:))
                  enddo
                  PrAB(x,y,z,:) = PrAB(x,y,z,:) * SUM(PrEW)
                endif

              else
                if (catB(j)==2) then
                  PrE = 1
                else if (catB(j)==7) then
                  call ParProb(l, GpID(3-kB,-Parent(BB(j),3-kB),3-kB), 3-kB, 0,0,PrH)
                  do e=1,3
                    PrE(e) = SUM(AKA2P(e,y,:) * PrH)
                  enddo              
                else if (catB(j)==8) then
                  call ParProb(l, GpID(3-kA,-Parent(BB(j),3-kB),3-kB), 3-kA, 0,0,PrH)
                  do e=1,3
                    PrE(e) = SUM(AKA2P(e,x,:) * PrH)
                  enddo
                else if (UseEE(nA+1)/=0 .and. Bj==B) then  
                  call ParProb(l, MateABpar(nA+1), 3-TypeEE(nA+1), 0,0,PrH)
                  do e=1,3
                    do u=1, 3
                      PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,x,UseEE(nA+1)) * PrH)
                    enddo
                    PrE(e) = SUM(PrW)
                  enddo
                  PrE = PrE/SUM(PrE)
                else
                  call ParProb(l,Parent(Bj,3-m),3-m,-1,0,PrE)
                endif
                if (Parent(Bj,3-m)<0 .and. Parent(Bj,3-m)/=GG(3-m)) then
                  do g=1, nS(-Parent(Bj, 3-m),3-m)
                    Ei = SibID(g,-Parent(Bj,3-m),3-m)
                    if (nFS(Ei) == 0) cycle
                    if (ANY(AA(1:nA)==Ei)) cycle
                    if (Parent(Ei,m)==GG(m) .and. GG(m)/=0) cycle
                    do e=1,3
                      call ParProb(l, Parent(Ei, m), m, Ei, -1, PrH)
                      do i=1, nFS(Ei)
                        if (ANY(AA(1:nA)==FSID(i, Ei))) cycle
                        if (FSID(i, Ei)==B) cycle
                        if (FSID(i, Ei)==PA) cycle
                        if (Genos(l,FSID(i, Ei))/=-9) then
                          PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                        endif
                      enddo
                      PrE(e) = PrE(e) * SUM(PrH)
                    enddo
                  enddo 
                endif
                do u=1, MAX(nFS(Bj),1)
                  if (FSID(u,Bj)==B) cycle
                  if (FSID(u, Bj)==PA) cycle
                  if (ANY(AA(1:nA)==FSID(u,Bj))) cycle
                  if (Genos(l,FSID(u, Bj))==-9) cycle
                  if (kB==3 .and. m==kA) then
                    PrE = PrE * OKA2P(Genos(l,FSID(u,Bj)), z, :) 
                  else
                    PrE = PrE * OKA2P(Genos(l,FSID(u,Bj)), y, :) 
                  endif
                enddo
                if (catB(j)==5) then
                  do e=1,3
                    PrE(e) = PrE(e) * SUM(AKA2P(y,e,:) * PrGG)
                  enddo
                endif
                
                if (kB<3) then
                  if (ANY(UseEE/=0) .or. ANY(catB<0)) then
                    if (catB(j)==2) then 
                      PrAB(x,y,z,1) = PrAB(x,y,z,1) * PrE(z)
                    else if (SUM(PrE)<3) then
                      PrAB(x,y,z,1) = PrAB(x,y,z,1) * SUM(PrE)
                    endif
                  else if (catB(j)==2) then 
                    PrAB(:,y,z,1) = PrAB(:,y,z,1) * PrE(z)
                  else if (SUM(PrE)<3) then
                    PrAB(:,y,z,1) = PrAB(:,y,z,1) * SUM(PrE)
                  endif
                  do u=1, MAX(nFS(Bj),1)
                    if (FSID(u,Bj)==B .and. Genos(l,B)/=-9) then
                      PrE = PrE * OKA2P(Genos(l,B), y, :) 
                    endif
                  enddo
                  
                  if (ANY(UseEE/=0) .or. ANY(catB<0)) then
                    if (catB(j)==2) then 
                      PrAB(x,y,z,2) = PrAB(x,y,z,2) * PrE(z)
                    else if (SUM(PrE)<3) then
                      PrAB(x,y,z,2) = PrAB(x,y,z,2) * SUM(PrE)
                    endif
                    if (Bj==B)  PrEE(:,x,nA+1) = PrE
                  else if (catB(j)==2) then 
                    PrAB(:,y,z,2) = PrAB(:,y,z,2) * PrE(z)
                  else if (SUM(PrE)<3) then
                    PrAB(:,y,z,2) = PrAB(:,y,z,2) * SUM(PrE)
                  endif
                else if (SUM(PrE)<3) then
                  PrAB(:,y,z,:) = PrAB(:,y,z,:) * SUM(PrE)
                endif
              endif
            enddo  ! j_m
          enddo  ! m
          if (kB==3 .and. Genos(l,B)/=-9) then
            PrAB(:,y,z,2) = PrAB(:,y,z,2) * OKA2P(Genos(l,B), y, z)
          endif
        enddo  ! z
      endif  ! kB==3
      enddo  ! x (ANY(UseEE/=0)  only)
    enddo  ! y

    PrL(l) = LOG10(SUM(PrAB(:,:,:,2))) - LOG10(SUM(PrAB(:,:,:,1)))
    PrLx(l,1) = LOG10(SUM(PrAB(:,:,:,1)))
    PrLx(l,2) = LOG10(SUM(PrAB(:,:,:,2)))
  endif
enddo
LL = SUM(PrL)

deallocate(UseEE)
deallocate(PrEE)
deallocate(MateABpar)
deallocate(TypeEE)

end subroutine PairUA 

! #####################################################################

subroutine addFA(A, SB, kB, LL)
use Global
implicit none

integer, intent(IN) :: A,SB,kB
double precision, intent(OUT) :: LL
integer :: x, y, z, Par, i, l, ParChk
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrY(3), PrZ(3)

LL = 999D0
call getFSpar(SB, kB, .TRUE., Par)  ! TODO: strict=FALSE ? 
if (Par/=0) then
  LL = 777
  return
endif
if (Par <0) then
  call getFSPar(-Par, 3-kB, .TRUE., ParChk)
  if (ParChk /= -SB) then
    LL = 444
    return
  endif
endif
if (Parent(A,kB)/=0) then
  LL = 444
  return
else if (Parent(A, 3-kB)/=0 .and. Parent(A, 3-kB)/=Par) then
  LL = 777
  return
endif

PrL = 0D0
do l=1, nSnp
  call ParProb(l, SB, kB, -1, 0, PrY)
  call ParProb(l, Par, 3-kB, -1, 0, PrZ)
  do x=1,3
    do y=1,3
      do z=1,3
        PrXYZ(x,y,z) = PrY(y) * PrZ(z) * AKA2P(x,y,z)
        if (Genos(l,A)/=-9) then
          PrXYZ(x,y,z) = PrXYZ(x,y,z) *  OKA2P(Genos(l,A), x, z)
        endif
        do i=1, ns(SB, kB)
          if (Genos(l,SibID(i,SB,kB))/=-9) then
            PrXYZ(x,y,z) = PrXYZ(x,y,z) *  OKA2P(Genos(l,SibID(i,SB,kB)), y, z)
          endif
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))
enddo

LL = SUM(PrL)

end subroutine addFA

! #####################################################################

subroutine pairCC(A,B,k, LL)  ! full 1st cousins
use Global
implicit none

integer, intent(IN) :: A,B,k
double precision, intent(OUT) :: LL
integer :: l, x, y, u, v, AreHS, z, AncA(2,mxA), AncB(2,mxA)
double precision :: PrL(nSnp, 3), PrXY(3,3), PrUV, PrPA(3), PrPB(3), &
  PrC(3,3,3), PrZ(3), PrXYf(3,3,2), LLself(2)
logical :: MaybeInbr(2)

LL = 999D0
if (Parent(A,k)/=0 .and. Parent(B,k)/=0) then
  if (Parent(A,k)==Parent(B,k)) then
    LL = 777
  endif
endif

if (Parent(A,k)/=0 .or. Parent(B,k)/=0) then  ! TODO? See ParentHFS
  LL = 444
endif
if (LL/=999) return

 call GetAncest(A,1,AncA)
 call GetAncest(B,1,AncB)
if (ANY(AncA == B) .or. ANY(AncB == A)) then
  LL = 777
  return
endif

AreHS = 0
MaybeInbr = .TRUE.
if (Parent(A,3-k)/=0 .and. Parent(A,3-k)==Parent(B,3-k)) then
    AreHS = 1
endif
if (hermaphrodites) then
  if (Parent(A,3-k)>0) then
    call PairSelf(B, Parent(A,3-k), LLself(1))
    call PairFullSib(B, Parent(A,3-k), LLself(2))
    if (LLself(1) > LLself(2)) then
      MaybeInbr(1) = .FALSE.
    endif
  endif
  if (Parent(B,3-k)>0) then
    call PairSelf(A, Parent(B,3-k), LLself(1))
    call PairFullSib(A, Parent(B,3-k), LLself(2))
    if (LLself(1) > LLself(2)) then
      MaybeInbr(2) = .FALSE.
    endif
  endif
endif

PrL = 0D0
do l=1, nSnp
  if (Genos(l,A)==-9 .and. Genos(l,B)==-9) cycle
  if (AreHS==0) then
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrPA)
    call ParProb(l, Parent(B,3-k), 3-k, B, 0, PrPB)
  else
    call ParProb(l, Parent(A,3-k), 3-k, A, B, PrPA)
  endif
  
  PrC = 0
  do u=1,3  ! GG 3-k
    do v=1,3  ! GG k
      PrUV = AHWE(u,l) * AHWE(v,l)
      do x=1,3  !PA
        do y=1,3    !PB
          PrXY(x,y) = AKA2P(x,u,v) * AKA2P(y,u,v) * PrUV
          if (AreHS==0) then
            PrXYf(x,y,1) = PrXY(x,y) * PrPA(u) / AHWE(u,l) 
            PrXYf(x,y,2) = PrXY(x,y) * PrPB(u) / AHWE(u,l)  
            if (Genos(l,A)/=-9) then
              PrXY(x,y) =PrXY(x,y) * SUM(OKA2P(Genos(l,A), x, :) * PrPA)
              PrXYf(x,y,1) = PrXYf(x,y,1) * OKA2P(Genos(l,A), x, u)  ! A inbred
              PrXYf(x,y,2) = PrXYf(x,y,2) * SUM(OKA2P(Genos(l,A), x, :) * PrPA)
            endif   
            if (Genos(l,B)/=-9) then
              PrXY(x,y) = PrXY(x,y) * SUM(OKA2P(Genos(l,B), y, :) * PrPB)
              PrXYf(x,y,1) = PrXYf(x,y,1) * SUM(OKA2P(Genos(l,B), y, :) * PrPB)
              PrXYf(x,y,2) = PrXYf(x,y,2) * OKA2P(Genos(l,B), y, u) 
            endif  ! not considered: both A & B inbred. 
          else if (AreHS==1) then
            PrZ = PrPA
            do z=1,3
              if (Genos(l,A)/=-9) PrZ(z) =PrZ(z) * OKA2P(Genos(l,A),x,z)
              if (Genos(l,B)/=-9) PrZ(z) =PrZ(z) * OKA2P(Genos(l,B),y,z)
            enddo
            PrXY(x,y) = PrXY(x,y) * SUM(PrZ)
          endif
        enddo
      enddo
      PrC(u,v,1) = SUM(PrXY)
      if (AreHS==0) then
        if (MaybeInbr(1))  PrC(u,v,2) = SUM(PrXYf(:,:,1))
        if (MaybeInbr(2))  PrC(u,v,3) = SUM(PrXYf(:,:,2))
      endif
    enddo
  enddo 
  do z=1,3
    PrL(l,z) = LOG10(SUM(PrC(:,:,z)))
  enddo
enddo
LL = MaxLL(SUM(PrL, DIM=1))

end subroutine pairCC

! #####################################################################

subroutine Clustering
use Global
implicit none

integer :: k, x, n, m, ij(2), sx(2), topX, u, fcl, par, topFS
double precision :: LL(7), dLL, LLx(7, 2), dLLtmp(maxSibSize)
logical :: FSM, UseAge(2), IsPair, OK  !, IsFSpair

FSM = .FALSE.
UseAge = (/.TRUE., .FALSE./)
do x=1, nPairs
  if (Error/=0) return
  if (MODULO(x,500)==0) call rchkusr()
  LL = 999D0
  ij = PairID(x,:)
  
  do k=1,2
    if (Parent(ij(1),k)>0 .or. Parent(ij(2),k)>0) cycle                                                  
    if (hermaphrodites .and. ANY(parent(ij,3-k)>0) .and. ANY(parent(ij,k)==0))  cycle                                                   
    sx(1) = -Parent(ij(1),k)  
    sx(2) = -Parent(ij(2),k)
    fcl = 3
    if (Complx==0) then
      fcl = 2
    else if (ALL(Parent(ij(1),:)==0) .and. ALL(Parent(ij(2),:)==0)) then
      fcl = 2
    else
      do m=1,2
        if (Parent(ij(m),k) < 0) then
          call getFSpar(sx(m), k, .TRUE., Par)  
          if (Par < 0) then
            if (ALL(Parent(ij(3-m),:)==0)) then
              fcl = 2
            else if (Parent(ij(3-m),k)<0) then
              call getFSpar(sx(3-m), k, .TRUE., Par) 
              if (Par < 0) then
                fcl = 2
              endif
            endif
          endif
        endif
      enddo
    endif
    if (fcl==2 .and. Complx>0 .and. AgeDiff(ij(1),ij(2))/=999) then  ! exception
      if (AgePriorM(ABS(AgeDiff(ij(1),ij(2)))+1, 9) == 0.0 .or. &
       AgePriorM(ABS(AgeDiff(ij(1),ij(2)))+1, 3-k) == 0.0) then
        fcl = 3       
      endif
    endif

    if (sx(1)==0 .and. sx(2)==0) then
      IsPair = .TRUE.                                                                   
      do n=1,2
        if (AgePhase==0 .and. n==1)  cycle
        topX = 0
        dLL = 0
        LLx = 999D0             
        if (AgeDiff(ij(1),ij(2))==999) then
          call CheckRel(ij(1), k, ij(2), k, UseAge(n), fcl, LLx(:,1))  
          call CheckRel(ij(2), k, ij(1), k, UseAge(n), fcl, LLx(:,2))
          do u=1,7
            LL(u) = MaxLL(LLx(u,:))
          enddo
        else if(AgeDiff(ij(1),ij(2))>=0) then
          call CheckRel(ij(1), k, ij(2), k, UseAge(n), fcl, LL)
        else
          call CheckRel(ij(2), k, ij(1), k, UseAge(n), fcl, LL)
        endif
        if (LL(2)>0 .and. LL(3)>0) then
          IsPair = .FALSE.
          exit
        endif
        call BestRel(LL, fcl, topX, dLL)
        if (topX==fcl .or. (fcl==3 .and. topX==2)) then
          IsPair = .TRUE.
        else if (AgePhase==2 .and. n==2) then
          if (MaxLL(LL((/2,3/))) - MaxLL(LL((/1,6,7/)))>TA .and. &
            MaxLL(LL((/4,5/))) - MaxLL(LL((/2,3/))) < TA) then
              IsPair = .TRUE.
          else
            IsPair = .FALSE.
            exit
          endif
        else
          IsPair = .FALSE.
          exit
        endif             
      enddo
      if (.not. IsPair) cycle
      
      nC(k) = nC(k)+1  ! new sibship (pair)              
      nS(nC(k), k) = 2
      SibID(1:2, nC(k), k) = ij
      Parent(ij(1), k) = -nC(k)
      call NewPar(ij(2), k, -nC(k))
      call CalcCLL(nC(k), k)
      if (fcl==2 .or. (topX==2 .and.  dLL>2*TA)) then
        if (ALL(Parent(ij, 3-k)==0)) then  ! another new sibship
          nC(3-k) = nC(3-k)+1  
          nS(nC(3-k), 3-k) = 2
          SibID(1:2, nC(3-k), 3-k) = ij
          Parent(ij(1), 3-k) = -nC(3-k)
          call NewPar(ij(2), 3-k, -nC(3-k))
          call CalcCLL(nC(3-k), 3-k)
        else
          do m=1,2
            if (Parent(ij(m),3-k)/=0) then
              call NewPar(ij(3-m), 3-k, Parent(ij(m),3-k))
            endif
          enddo
        endif
      endif
      if (Parent(ij(1), 3-k) /=0 .and. &
        Parent(ij(1), 3-k)==Parent(ij(2), 3-k)) then
        call MakeFS(ij(1), ij(2))
      endif
      call CalcCLL(nC(k), k)
      do n=1, 2
        u = SibID(n, nC(k), k)
        if (Parent(u,3-k) < 0) then
          call CalcCLL(-Parent(u,3-k), 3-k)
        endif                
        call CalcLind(SibID(n, nC(k), k))
      enddo
      call CalcCLL(nC(k),k)
      
    else if (sx(1)>0 .and. sx(2)>0 .and. sx(1) /= sx(2)) then
      IsPair = .TRUE.
      do n=1,2
        if (AgePhase==0 .and. n==1)  cycle                                  
        if (AgePhase==2 .and. n==2)  cycle
        call CheckMerge(sx(1), sx(2), k,k, UseAge(n), LL,1, FSM)  
        call BestRel(LL, 1, topX, dLL)
        if (topX /=1 .or. dLL < TA * MIN(nS(sx(1),k), nS(sx(2),k)) &
         .or. (fcl==2 .and. (.not. FSM .or. &
          dLL < 2*TA * MIN(nS(sx(1),k), nS(sx(2),k))))) then
          IsPair = .FALSE.
          exit
        endif
      enddo
      if (.not. IsPair) cycle      
      
      if (FSM .and. fcl==2) then
        call DoFSmerge(sx(1), sx(2), k)
      else 
        call DoMerge(sx(1), sx(2), k, OK)
      endif
      
    else
      do m=1,2
        if (sx(m)>0 .and. sx(3-m)==0) then
          IsPair = .TRUE.
          topX = -9                   
          do n=1,2                                                                       
            if (AgePhase==0 .and. n==1)  cycle                                 
            if (AgePhase==2 .and. n==2)  cycle
            call CheckRel(ij(3-m), 0, -sx(m), k, UseAge(n), fcl, LL)
            call BestRel(LL, fcl, topX, dLL)
            if (.not. (topX==fcl .or. (fcl==3 .and. topX==2))) then
              IsPair = .FALSE.
              exit                               
            endif 
          enddo
          if (.not. IsPair) cycle
          if (ANY(SibID(1:ns(sx(m),k), sx(m), k) == Parent(ij(3-m),3-k))) then  ! inbreeding
            call CheckPair(ij(3-m), Parent(ij(3-m),3-k), k, 3, LLx(:,1), LLx(:,2))
            call BestRel(LL, 3, topX, dLL)
            if (topX /= 3) then
              IsPair = .FALSE. 
              cycle
            endif
          endif
          dLLtmp = 999
          FSM = .FALSE.
          if (topX==2 .and. fcl/=2) then  
            call BestRel(LL, 2, topX, dLL)  
            if (topX==2 .and. dll > 2*TA) then 
              FSM = .TRUE.
            endif
          endif
          if (fcl==2 .or. FSM) then                         
            call getFSpar(sx(m), k, .TRUE., Par)
            if (Par/=0) then
              call NewPar(ij(3-m), 3-k, Par)
            else
              call AddFS(ij(3-m), sx(m), k,0,k, LL(2), topFS, dLLtmp)
              if (topFS>0 .and. MAXVAL(dLLtmp, mask=dLLtmp<777)>2*TA) then
                call CheckPair(ij(3-m), topFS, k, 2, LLx(:,1), LLx(:,2))
                call BestRel(LLx(:,2), 2, topX, dLL)
                if (topX==2 .and. dll > 2*TA) then                
                  if (Parent(topFS, 3-k)/=0) then                
                    call NewPar(ij(3-m), 3-k, Parent(topFS, 3-k))
                  else if (Parent(ij(3-m), 3-k)/=0) then
                    call NewPar(topFS, 3-k, Parent(ij(3-m), 3-k))                                     
                  else
                    nC(3-k) = nC(3-k)+1  ! new sibship (pair)      
                    nS(nC(3-k), 3-k) = 2
                    SibID(1:2, nC(3-k), 3-k) = (/ ij(3-m), topFS /)
                    Parent(ij(3-m), 3-k) = -nC(3-k)
                    call NewPar(topFS, 3-k, -nC(3-k))
                    call CalcCLL(nC(3-k), 3-k)
                  endif
                endif
              endif
            endif
          endif
          call DoAdd(ij(3-m), sx(m), k)
        endif
      enddo
    endif
  enddo
enddo

end subroutine Clustering

! #####################################################################

subroutine Merging ! check if any of the existing clusters can be merged
use Global
implicit none

integer :: k, s, r, topX, xr, ParOpp(2), n
double precision :: LLm(7), dLL
logical :: FSM, OK, UseAge(2), MergeOK

FSM = .FALSE.
UseAge = (/.FALSE., .TRUE./)
do k=1,2
  if (Complx==0 .and. k==2) cycle
  do s=1,nC(k)-1
    r = s
    if (Error/=0) return
    if (MODULO(s,20)==0) call rchkusr()                       
    do xr=s+1, nC(k)
      r = r + 1
      if (r > nC(k)) exit   ! possible due to merged sibships
      topX = 0
      OK = .TRUE.
      do n=1,2
        if (AgePhase==2 .and. n==1)  cycle                                 
        call CheckMerge(s, r, k, k, UseAge(n), LLm, 1, FSM)
        call BestRel(LLm, 1, topX, dLL)
        if ((n==1 .and. MaxLL(LLm(2:7)) - LLm(1) > TA) .or. &
         (n==2 .and. (topX/=1 .or. dLL < TA * MIN(nS(s,k), nS(r,k))))) then
          OK = .FALSE.
          exit
        endif
      enddo
      if (.not. OK)  cycle
        
      ParOpp(1) = Parent(SibID(1,s,k), 3-k)
      ParOpp(2) = Parent(SibID(1,r,k), 3-k)
      MergeOK = .TRUE.
      if (FSM .and. dLL > 2*TA * MIN(nS(s,k), nS(r,k))) then
         call DoFSmerge(s, r, k)
      else
          call DoMerge(s, r, k, MergeOK)
      endif
      call CalcCLL(s,k)
      if(MergeOK)  r = r-1  ! otherwise a cluster is skipped
    enddo
  enddo
enddo

end subroutine Merging

! #####################################################################

subroutine GrowClusters
! for each individual, check if they can be added to any sibship cluster
use Global
implicit none

integer :: k, s, i, nMaybe, ClM(100), topX, n, fcl, j, ctop(100), opPar, topFS
double precision :: LLtmp(7), dLL, dLLM(100), ALR, dLLtmp(maxSibSize), LR, LLtmpG(7)
logical :: UseAge(2), Addi

UseAge = (/.TRUE., .FALSE./)
do k=1,2   
  if (Complx==0 .and. k==2) exit
  do i=1, nInd
    if (Error/=0) return
    if (MODULO(i,500)==0) call rchkusr()
    if (Parent(i,k)/=0) cycle                                                                                       
    nMaybe = 0
    Clm = 0
    dLLM = 999D0
    do s=1,nC(k)  ! sibships
      if (ANY(GpID(:,s,k)==i)) cycle
      if (BY(i)>=0 .and. ALL(BY(SibID(1:ns(s,k),s,k))>=0)) then
        if(ANY(AgePriorM(ABS(AgeDiff(SibID(1:ns(s,k),s,k),i))+1, k)==0.0)) cycle
      endif
      call CalcAgeLR(i,0,-s,k,k,1, .TRUE., ALR)
      if (ALR==777 .or. ALR < 3*TF)  cycle
      call Qadd(i, s, k, LR)
      if (LR < ns(s,k)*TF) cycle                     
      LLtmp = 999D0
      Addi=.TRUE.
      if (Complx==0 .or. (Hermaphrodites .and. Parent(i,3-k)==0)) then
        fcl = 2 ! FS
      else if (Parent(i,3-k)==0 .and. ALL(Parent(SibID(1:ns(s,k),s,k),3-k)<=0)) then
         call getFSpar(s, k, .FALSE., opPar)
        if (opPar < 0) then
          fcl = 2
        else
          fcl = 3
        endif  
        do j=1,ns(s,k)
          call CalcAgeLR(i,0,j,0,0,2, .TRUE.,ALR)
          if (ALR == 777) then
            fcl = 3   ! not FS with all
            exit
          else if (AgePhase==2 .and. AgeDiff(i,j)/=999) then
            if (AgePriorM(ABS(AgeDiff(i,j))+1, k) - &
              AgePriorM(ABS(AgeDiff(i,j))+1, 3-k) > 3*TF) then
              fcl = 3
              exit
            endif  
          endif
        enddo  
      else
        fcl = 3                                   
      endif
      
      do n=1,2
        if (AgePhase==0 .and. n==1)  cycle
        if (AgePhase==2 .and. n==2)  cycle
        topX = 0
        dLL = 0
        call CheckRel(i, 0, -s, k, UseAge(n), fcl, LLtmp)
        call BestRel(LLtmp, fcl, topX, dLL)
        if (.not. ((topX==2 .and. dll>2*TA) .or. &
          fcl==topX .or. (fcl==3 .and. topX==2))) then
          Addi = .FALSE.
          exit
        endif
      enddo
      if (ANY(SibID(1:ns(s,k), s, k) == Parent(i,3-k))) then  ! inbreeding
        call CheckPair(i, Parent(i,3-k), k, 3, LLtmpG, LLtmp)
        call BestRel(LLtmp, 3, topX, dLL)
        if (topX /= 3) then
          Addi = .FALSE. 
        endif
      endif
      if (Addi) then
        nMaybe = nMaybe+1
        Clm(nMaybe) = s 
        dLLM(nMaybe) = dLL
        ctop(nMaybe) = topX
      endif
    enddo 
    
    if (nMaybe==1) then
      if (fcl==2 .and. ctop(1)==2 .and. dLLM(1)>2*TA) then
        call DoAdd(i, Clm(1), k)                                                  
        call getFSpar(Clm(1), k, .TRUE., opPar)
        if (opPar/=0) then
          call NewPar(i, 3-k, opPar)
        else
          call AddFS(i, Clm(1), k,0,k, LLtmp(2), topFS, dLLtmp)
          if (Complx==0 .or. (topFS>0 .and. MAXVAL(dLLtmp, mask=dLLtmp<777)>2*TA)) then
            if (Parent(topFS, 3-k)/=0) then
              call NewPar(i, 3-k, Parent(topFS, 3-k))                  
            else if (Parent(i, 3-k)/=0) then
              call NewPar(topFS, 3-k, Parent(i,3-k))
            else
              nC(3-k) = nC(3-k)+1  ! new sibship (pair)      
              nS(nC(3-k), 3-k) = 2
              SibID(1, nC(3-k), 3-k) = i
              SibID(2, nC(3-k), 3-k) = topFS
              Parent(i, 3-k) = -nC(3-k)
              call NewPar(topFS, 3-k, -nC(3-k))
              call CalcCLL(nC(3-k), 3-k)              
            endif
          endif
        endif  
      else if (fcl==3) then
        call DoAdd(i, Clm(1), k) 
      endif
    endif                       
  enddo
enddo
  
end subroutine GrowClusters

! #####################################################################

subroutine SibParent  
! for each sibship, check if a parent can be assigned
use Global
implicit none

integer :: k, s, xs, i, n,maybe, topX, CurNumC(2), OH, Par, MaybeOpp, &
  j, nCandPar, CandPar(20), h, SibTmp(maxSibSize), nSib, sib1
double precision :: LL(7), dLL, LLtmp(7), ALR, LLO, LR
logical :: OK             

 CurNumC = nC 
maxOppHom = MaxMismatch - FLOOR(-nSNP * Er)                                                                        
do k=1,2
  s = 0
  do xs=1, CurNumC(k)
    s = s+1
    if (Error/=0) return
    if (MODULO(s,20)==0) call rchkusr()
    if (s > nC(k)) exit   
    nCandPar = 0
    if (hermaphrodites) then
      call getFSpar(s, k, .TRUE., Par)
      if (Par < 0)  cycle
    endif  
    do i=1,nInd                                                                 
      if (nCandPar == 20) exit  !unlikely
      if (Sex(i)/=k .and. Sex(i)/=3) cycle
      if (Parent(i,k)==-s) cycle
      if (GpID(1,s,k)==i .or. GpID(2,s,k)==i) cycle
      if (ANY(AgeDiff(SibID(1:ns(s,k), s, k), i) <= 0))  cycle                                                        
      maybe=1
      call CalcAgeLR(-s, k, i, k, k, -1, .TRUE., ALR)
      if (ALR==777 .or. ALR < 3*TF)  cycle                             
      do n=1,nS(s,k)
        call CalcOH(i, SibID(n,s,k), OH)
        if (OH > maxOppHom) then
          maybe = 0
          exit
        endif   
      enddo                                                                         
      if (maybe==0) cycle
      do n=1,2
        if (GpID(n,s,k) <= 0) cycle
        call CalcOH(i, GpID(n,s,k), OH)
        if (OH > maxOppHom) then
          maybe = 0
          exit
        endif   
      enddo                                                                         
      if (maybe==0) cycle
      call QPO(i, s, k, LR)
      if (LR < TF*nS(s,k))  cycle                          
       LL = 999D0                                                                 
      call CheckRel(-s, k, i, k, AgePhase > 1, 1, LL)
      call BestRel(LL, 1, topX, dLL)
      if (topX/=1) then
        maybe = 0
        cycle
      endif
      if (Sex(i)==3 .or. hermaphrodites) then  ! check if parent of opposite sex instead
        MaybeOpp = 1
        call getFSpar(s, k, .TRUE., Par)                       
        if (Par > 0) then
          MaybeOpp = 0
        else if (Par/=0 .and. Parent(i, 3-k) == Par) then
          MaybeOpp = 0  ! are HS
        else if (Par==0) then
          if (ANY(Parent(SibID(1:nS(s,k), s,k),3-k)>0)) then
            MaybeOpp = 0
          else  ! check if could all be FS
            do j=1, ns(s,k)-1
              do h=j+1, ns(s,k)
                call CalcAgeLR(sibID(j,s,k), 0, SibID(h,s,k), 0, 0, 2, .TRUE., ALR)
                if (ALR == 777 .or. ALR < 3*TF) then
                  MaybeOpp = 0
                  exit
                endif                
              enddo
              if (MaybeOpp == 0) exit
            enddo
          endif
          if (MaybeOpp == 1) then
            call OppMerge(s,k,LLO)
            if (LLO>444 .or. (CLL(s,k) - LLO) > ns(s,k)*TF) then
              MaybeOpp = 0
            endif
          endif           
          else if (Par < 0) then
          do n=1,nS(-Par,3-k)
            call CalcOH(i, SibID(n,-par,3-k), OH)
            if (OH > maxOppHom) then
              maybeOpp = 0
              exit
            endif   
          enddo
          if (maybeOpp == 1) then
            do n=1,2
              if (GpID(n,-Par,3-k) <= 0) cycle
              call CalcOH(i, GpID(n,-Par,3-k), OH)
              if (OH > maxOppHom) then
                maybeOpp = 0
                exit
              endif   
            enddo
          endif
        endif 
        if (MaybeOpp == 1) then
          if (Par < 0) then  ! may have more/fewer sibs
            call CheckRel(Par, 3-k, i, 3-k, AgePhase > 1, 1, LLtmp)
            call BestRel(LLtmp, 1, topX, dLL)
            if (topX/=1)  MaybeOpp = 0
          else if (Par == 0 .and. ns(s,k)>0) then
            sib1 = SibID(1,s,k)
            call PairPO(sib1, i, 3-k, 0, LLtmp(1))   
            call CalcU(sib1, k, i, 3-k, LLtmp(2))
            if (LLtmp(1)>0 .or. (LL(1)-LL(7)) - (LLtmp(1)-LLtmp(2)) > &
             TA*ns(s,k))  MaybeOpp = 0
          endif
        endif
        if (MaybeOpp==1) cycle
      endif     
      if (maybe==1) then               
        nCandPar = nCandPar + 1
        CandPar(nCandPar) = i
      endif
    enddo  ! i
    
    if (nCandPar == 1) then        
      SibTmp = SibID(:,s,k)       
      nSib = ns(s,k)
      do n=1,nSib 
        call NewPar(SibTmp(n), k, CandPar(1))                                         
      enddo  
      if (Sex(CandPar(1))==3)  Sex(CandPar(1)) = k
      call DoMerge(0, s, k, OK)  !removes cluster s, 
      s = s-1  ! otherwise a cluster is skipped
    endif
  enddo ! s
enddo ! k
  
end subroutine SibParent

! #####################################################################

subroutine MoreParent(pre)  
! for each individual, check if a parent can be assigned now.
use Global
implicit none

logical, intent(IN) :: pre
integer :: i, j, l, OH, Lboth, k, curPar(2), AncJ(nInd,2,mxA), topX
double precision :: LRP, LLP(7), curLL, ALR(2)
logical :: dropped, OK

do j=1,nInd
  call GetAncest(j,1,AncJ(j,:,:))  ! store globally? how time consuming?
enddo

do i=1, nInd
  if (Error/=0) return
  if (MODULO(i,500)==0) call rchkusr()
  call CalcLind(i)
  CurPar = Parent(i,:)
  curLL = Lind(i)
  dropped = .FALSE.        
    if (ANY(Parent(i,:)>0) .and. ANY(Parent(i,:)<=0)) then  ! double check
      do k=1,2
        if (Parent(i,k)<=0)  cycle
        if (AgeDiff(i,Parent(i,k)) == 999)  cycle
        if (Parent(i,3-k) < 0) then
          if (ns(-Parent(i,3-k),3-k)>2) then
            call NewPar(i, k, 0)      
            call CheckRel(i, k, CurPar(k), k, .TRUE., 1, LLP)
            call BestRel(LLP, 1, topX, LRP)
            if (topX==1 .and. LRP>TA) then  ! restore                          
              call NewPar(i, k, curPar(k))
          endif
        endif        
      endif
      if (Parent(i,k)==0) then
        if (Complx == 0) then
          call NewPar(i, 3-k, 0)
        else if (Parent(i,3-k) < 0) then
          call NewPar(i, 3-k, 0)
          dropped = .TRUE.
          call CheckRel(i, 0, curPar(3-k), 3-k, AgePhase>0, 1, LLP)
          call BestRel(LLP, 3, topX, LRP)
          if ((topX==3 .or. topX==2) .and. LRP>TA) then
            call NewPar(i, 3-k, curPar(3-k))
            dropped = .FALSE.
          endif
        endif
        call CalcLind(i)
      endif
      enddo
    endif

  if (ALL(Parent(i,:)/=0)) cycle                                                   
  if (pre)  cycle                                                    
  do j=1, nInd  ! candidate parent.
    if (i==j) cycle
    if (ANY(Parent(j,:)==i) .or. ANY(Parent(i,:)==j)) cycle
    if (AgeDiff(i,j) <= 0)  cycle  ! note: unknown = 999D0 > 0   
    if (ANY(AncJ(j,:,:) == i)) cycle
    if (Sex(j) == 3) then
      if (ALL(Parent(i,:)==0)) cycle 
    else
      if (Parent(i, Sex(j)) < 0) cycle 
    endif
    OH = 0
    do l=1,nSnp
      if ((Genos(l,i)==1).and.(Genos(l,j)==3)) then
        OH = OH+1
        if (OH > maxOppHom) exit
      endif                       
      if ((Genos(l,i)==3).and.(Genos(l,j)==1)) then
        OH = OH+1
        if (OH > maxOppHom) exit
      endif                       
    enddo
    if (OH > maxOppHom) cycle  
    if (OH <= maxOppHom) then
      Lboth = COUNT(Genos(:,i)/=-9 .and. Genos(:,j)/=-9)
      if (Lboth < nSnp/4.0)  cycle   ! >3/4th of markers missing
    endif               
    call CalcAgeLR(i,sex(i), j, sex(j), 0, 1, .TRUE., ALR(1))
    if (ALR(1) == 777 .or. ALR(1)<= 3*TF)  cycle
    call CalcAgeLR(j, sex(j), i,sex(i), 0, 1, .TRUE., ALR(2))
    if (ALR(2) /= 777 .and. (ALR(1)-ALR(2)) < TF)  cycle
    call CalcPO(i, j, LRP)
    if (LRP < 5*TF) cycle 
     if (Sex(j)/=3) then
      call PairPO(i, j, sex(j), 0, LLP(1))
    else
      call PairPO(i, j, 1, 0, LLP(1))
    endif  
    if (LLP(1) > 0) cycle
    call CalcU(i,sex(i), j, sex(j), LLP(7))
    if ((LLP(1) - LLP(7)) < TF) cycle
    call CalcPOZ(i,j, AgePhase>0)  ! assigns parent as side effect;
    if (CurPar(1)/=Parent(i,1) .or. CurPar(2)/=Parent(i,2)) then
      call GetAncest(j,1,AncJ(j,:,:))
      call GetAncest(i,1,AncJ(i,:,:))
      do k=1,2
        if(CurPar(k)>0)  call GetAncest(CurPar(k),k,AncJ(curPar(k),:,:))
      enddo
    endif
  enddo  ! j
  
  ! restore original parent?
  do k=1,2
    if (curPar(k)<0 .and. Parent(i,k)==0) then
      call CheckRel(i, 0, curPar(k), k, AgePhase>0, 1, LLP)
      call BestRel(LLP, 3, topX, LRP)
      if (topX==3 .or. topX==2) then
        call NewPar(i, k, curPar(k))
      else if (ns(-curPar(k), k)==0 .or. (ns(-curPar(k),k)==1 .and. &
         ALL(GpID(:,-curPar(k), k)==0))) then
        call DoMerge(0, -curPar(k), k, OK)  !remove cluster
      endif
    endif
  enddo

  if ((Parent(i,1)/=0 .and. curPar(1)==0) .or. (Parent(i,2)/=0 .and. &
   curPar(2)==0)) then
    do k=1,2
      if (Parent(i,k)<=0) cycle
      if (AgeDiff(i, Parent(i,k)) == 999 .or. Sex(Parent(i,k))==3) then
        if (Parent(i,3-k) == 0) then
          call NewPar(i, k, 0) ! reset
        else if (Parent(i,3-k)>0) then
          if (Sex(Parent(i,k))==3 .and. Sex(Parent(i,3-k))==3) then
            call NewPar(i, k, 0) 
            call NewPar(i, 3-k, 0) 
          endif
        endif
      endif
      if (Parent(i,k)/=0) then
        Sex(Parent(i,k)) = k
      endif
    enddo
  endif 
enddo

end subroutine MoreParent

! #####################################################################

subroutine SibGrandparents 
! for each sibship, find the most likely male and/or female grandparent
use Global
implicit none

integer :: k, s, i,j, r, m, par, xs, candGP(20, 2), nCG(2), curGP(2),&
   u, v, AncI(nInd,2,mxa), AncR(maxval(nC),2,2,mxA), GGP(2), GpIN(2), ix, &
   Best(2), x, sib1
double precision :: LRG, ALRtmp(2), LLX(3,2), curCLL, &
 LLGGP(2,2), dx(maxSibSize), LL(7),LLG(21,21),gLL
logical :: SexUnk(20,2), OK, sharing, skip(nInd/2,2) 

skip = .FALSE.
do k=1,2
  do s=1, nC(k)
    if (ALL(Parent(SibID(1:nS(s,k), s, k), 3-k) < 0)) then
      call getFSpar(s, k, .TRUE.,par)
      if (par < 0) then
        if (nS(-par, 3-k) == nS(s,k)) then  !cannot tell if mat or pat
          skip(s,k) = .TRUE.
        endif
      endif          
    endif
    call GetAncest(-s, k, AncR(s,k,:,:))
  enddo
enddo

do i=1, nInd
  call GetAncest(i,1,AncI(i,:,:))
enddo
  
do k=1,2
  s = 0
  do xs=1, nC(k)
    s = s+1
    if (Error/=0) return
    if (MODULO(s,5)==0) call rchkusr()
    if (s > nC(k)) exit
    if (ALL(GpID(:,s,k)/=0) .and. ns(s,k)>2) cycle   
    if (skip(s,k))  cycle
    nCG = 1  ! first pair = 0,0
    CandGP = 0
    SexUnk = .FALSE.
    CurGP = GpID(:,s,k)
    call calcCLL(s,k)                 
    curCLL = CLL(s,k)
    do m=1,2
      if (GpID(m,s,k)/=0) then
        nCG(m) = 2
        CandGP(2,m) = GpID(m,s,k)
      endif
    enddo
    GpID(:,s,k) = 0
    call CalcCLL(s,k)
    
    do i=1,nInd
      if (Parent(i,k)==-s) cycle
      if (CandGP(2,1)==i .or. CandGP(2,2)==i) cycle
      if (Sex(i)/=3) then
        if (nCG(Sex(i))==20)  cycle
        if (ANY(curGP/=0) .and. .not. hermaphrodites) then
          if (curGP(Sex(i)) > 0) cycle
          if (curGP(Sex(i))<0) then
            if (ns(-curGP(Sex(i)), Sex(i)) > 1) cycle
          endif  
        endif
      else
        if (ANY(nCG==20))  cycle
      endif   
      if (ANY(AgeDiff(SibID(1:ns(s,k), s, k), i) <= 1))  cycle
      call CalcAgeLR(-s,k, i,Sex(i), 0,1, .FALSE., ALRtmp(1))
      if (ALRtmp(1) == 777 .or. ALRtmp(1) < 3*TF) cycle
      call CalcAgeLR(i,Sex(i), -s,k, 0,1, .FALSE., ALRtmp(2))
      if (ALRtmp(2)/=777 .and. (ALRtmp(1)-ALRtmp(2)) < TF)  cycle
      if (ALL(ALRtmp==0.0))  cycle  ! no age info      
      if (nS(s,k)>1) then
        call QGP(i, 0, s, k, LRG) 
        if (LRG < TF*nS(s,k))  cycle    !  
      else if (ns(s,k)==1) then
        sib1 = SibID(1,s,k)
        call PairQHS(i, sib1, LRG)
        if (LRG < 2*TF)  cycle
      endif
      if (ANY(AncI(i,k,:) == -s)) cycle
      LLX = 999D0
      call AddGP(i, s, k, LLX(1,1))
      if (LLX(1,1)>0) cycle
      call CalcU(i, k, -s, k, LLX(1,2))
      if ((LLX(1,1) - LLX(1,2)) < TA)  cycle   
      call AddGGP(i, s, k, LLX(2,2))    ! i GP of dummy?
      if ((LLX(1,1) - MaxLL(LLX(:,2))) < TF)  cycle 
      call pairUA(-s, i, k, 3, LLX(3,2))    ! i FS of dummy?  (todo?: PHS, MHS)
      if ((LLX(1,1) - MaxLL(LLX(:,2))) < TF)  cycle 
      do m=1,2
        if (Sex(i)/=3 .and. Sex(i)/=m) cycle
        if (ncG(m) < 20) then  ! arbitrary threshold to limit comp. time
          ncG(m) = nCG(m) + 1
          CandGP(nCG(m), m) = i
          if (Sex(i)==3) SexUnk(ncG(m),m) = .TRUE.
        else
          exit  
        endif
      enddo
    enddo
    
    do m=1,2
      if (curGP(m)>0 .and. ns(s,k)>2) cycle
      if (curGP(m)<0 .and. ns(s,k)>2) then
        if (ns(-curGP(m),m) > 2) cycle
      endif
      do r=1, nC(m) 
        if (m==k .and. s==r) cycle
        if (CandGP(2,m) == -r) cycle  ! current GP
        if (ANY(AncR(r,m,k,:)==-s)) cycle
        do j=1, nS(s,k)
          if (ANY(AncR(r,m,:,:) == SibID(j,s,k))) cycle ! covers 1 extra gen
        enddo
        if (nS(r,m)==1 .and. ANY(SibID(1:nS(s,k),s,k) == SibID(1,r,m))) cycle
        if (m/=k .and. complx<2) then
          if (ALL(Parent(SibID(1:ns(s,k),s,k),m) == -r))  cycle
          if (ALL(Parent(SibID(1:ns(r,m),r,m),k) == -s))  cycle
        endif
        if (ANY(GpID(:,r,m)>0)) then
          do u=1,2
            if (GpID(u,r,m)>0) then
              if (ANY(AgeDiff(SibID(1:ns(s,k), s, k), GpID(u,r,m)) <= 2))  cycle
            endif
          enddo
        endif
        sharing = .FALSE.
        do i=1, ns(s,k)
          if (m/=k .and. ANY(SibID(:,r,m) == SibID(i,s,k))) then
            sharing = .TRUE.
            exit
          else if (m==k .and. Parent(SibID(i,s,k),3-k)<0) then
            if (ANY(Parent(SibID(1:ns(r,m),r,m), 3-m) == Parent(SibID(i,s,k),3-k))) then
              sharing = .TRUE.
              exit
            endif
          endif
        enddo  
        if (.not. sharing) then  ! QGP doesn't work when inbreeding loops
          call QGP(-r, m, s, k,  LRG) 
          if (LRG < TF*MIN(nS(s,k), nS(r,m))) cycle  ! conservative.
        endif
        call CalcAgeLR(-s,k, -r,m, 0,1, .FALSE., ALRtmp(1))
        if (ALRtmp(1) == 777 .or. ALRtmp(1) < 2*TF) cycle
        call CalcAgeLR(-r,m, -s,k, 0,1, .FALSE., ALRtmp(2))
        if (ALRtmp(2)/=777 .and. (ALRtmp(1)-ALRtmp(2)) < TF)  cycle  ! what threshold?
        if (ALL(ALRtmp==0.0))  cycle  ! no age info
        LLX = 999D0
        call PairUA(-s, -r, k, m, LLX(1,1))
        if (LLX(1,1)>0) cycle
        call CalcU(-s,k, -r,m, LLX(1,2)) 
        if ((LLX(1,1) - LLX(1,2)) < nS(s,k)*TF)  cycle   
        call addFS(0, r, m, s, k, LLX(2,1), ix, dx) 
        if ((MaxLL(LLX(:,1)) - LLX(1,2)) < TA)  cycle
        call dummyGP(s, r, k, m, LLX(2,2))  ! SB GGP of A's
        if (MaxLL(LLX(:,1)) - MaxLL(LLX(:,2)) < TA)  cycle   ! nS(s,k)*TF 
        if (ncG(m) < 20) then
          nCG(m) = nCG(m) + 1
          CandGP(nCG(m), m) = -r
          if (COUNT(nFS(SibID(1:nS(r,m),r,m))>0)==1 .and. &
            ALL(Parent(SibID(1:nS(r,m),r,m), 3-m)<=0)) then
              SexUnk(nCG(m),m) = .TRUE.
          endif
        else
          exit 
        endif                
      enddo
    enddo
    
    LLG = 999D0
    LLG(1,1) = 0D0  ! no GP                       
    if (ALL(nCG==1)) then
      cycle  ! no candidate GPs
    else if (ANY(nCG==1)) then  ! cand GP of 1 sex only
      do m=1,2
        if (nCG(m)>1) then
          if (nCG(m)==2 .and. CandGP(2,m)==curGP(m)) then
            GpID(m,s,k) = curGP(m)  ! restore
          else 
            do u=2, nCG(m)  ! first pair = 0,0
              call CalcGPz(s, k, CandGP(u,m), m, gLL)
            enddo
          endif
        endif
      enddo
    else if (ALL(nCG==2)) then
      GpID(1,s,k) = CandGP(2,1)
      call CalcGPz(s, k, CandGP(2,2), 2, gLL)
      if (GpID(2,s,k) /= CandGP(2,2)) then
        GpID(1,s,k) = curGP(1)  ! undo temp. assignment
        call CalcGPz(s, k, CandGP(2,1), 1, gLL)
      endif   
    else ! if (all(nCG>20)) then   ! test all combo's
      if (.not. ANY(SexUnk)) then
        do x=2, MAXVAL(nCG)
          do m=1,2
            if (CandGP(x,m) == 0) cycle
            GpIN = GPID(:,s,k)
            u = which(CandGP(:,1), GpIN(1))
            v = which(CandGP(:,2), GpIN(2))
            if (m==1)  LLG(x, v) = 505
            if (m==2)  LLG(u, x) = 505
            call CalcGPz(s, k, CandGP(x, m), m, gLL)
            if (GpID(m,s,k) == CandGP(x,m)) then
              if (GpID(3-m,s,k) == GpIN(3-m)) then
                if (m==1)  LLG(x, v) = gLL
                if (m==2)  LLG(u, x) = gLL
              else if (GpID(3-m,s,k) == 0) then
                if (m==1 .and. LLG(x,1)>500)  LLG(x, 1) = gLL
                if (m==2 .and. LLG(1,x)>500)  LLG(1, x) = gLL
              endif
            else if (GpID(1,s,k) == GpIN(1) .and. GpID(2,s,k)==0) then
              if (LLG(u,1)>500)  LLG(u, 1) = gLL
            else if (GpID(2,s,k) == GpIN(2) .and. GpID(1,s,k)==0) then 
              if (LLG(1,v)>500)  LLG(1, v) = gLL
            endif
          enddo
        enddo
      endif    
      do u = 2, nCG(1)    
        do v = 2, nCG(2)
          if (SexUnk(u,1) .and. SexUnk(v,2))  cycle
          if (CandGP(u, 1)==CandGP(v, 2) .and. CandGP(u,1)>0) cycle 
          if (GpID(1,s,k) == CandGP(u,1)) then
            call CalcGPz(s, k, CandGP(v, 2), 2, gLL)
          else if (GpID(2,s,k) == CandGP(v,2)) then
            call CalcGPz(s, k, CandGP(u, 1), 1, gLL)
          else
            GPIN = GpID(:,s,k)
            GpID(1,s,k) = CandGP(u,1)
            call CalcGPz(s, k, CandGP(v, 2), 2, gLL)
            if (GPIN(1)/=0 .and. GpID(1,s,k) == 0) then  ! restore GP1?
              if (GpID(2,s,k) == GPIN(2)) then
                GPID(1,s,k) = GPIN(1)  
              else if (GPIN(2)/=0 .and. GpID(2,s,k)==0) then
                GPID(:,s,k) = GPIN
              else
                call CalcGPz(s, k, GPIN(1), 1, gLL)
              endif
            endif
          endif
          if (GpID(1,s,k)==CandGP(u,1) .and. GpID(2,s,k)==CandGP(v,2)) then
            LLG(u,v) = gLL
          else if (GpID(1,s,k)==CandGP(u,1) .and. GpID(2,s,k)==0) then
            if (LLG(u,1)>500)  LLG(u,1) = gLL
          else if (GpID(2,s,k)==CandGP(v,2) .and. GpID(1,s,k)==0) then
            if (LLG(1,v)>500)  LLG(1,v) = gLL
          endif
        enddo
      enddo
      Best = MAXLOC(LLG, MASK=LLG<505)
      if (ALL(Best==1)) then
        GpID(:,s,k) = 0
      else 
        if (Best(1) >1 .and. GpID(1,s,k)/=CandGP(Best(1),1) .and. &
         Best(2) >1 .and. GpID(2,s,k)/=CandGP(Best(2),2)) then
          GpID(1,s,k) = CandGP(Best(1),1)
          call CalcGpz(s,k, CandGP(Best(2),2), 2, gLL)
          if (GpID(2,s,k) /= CandGP(Best(2),2) .or. &
           GpID(1,s,k) /= CandGP(Best(1),1)) then
            GpID(:,s,k) = 0
            call CalcCLl(s,k)
          endif
        endif
        do m=1,2
          if (GpID(m,s,k)/=CandGP(Best(m),m)) then
            if ((SexUnk(Best(m),m) .and. GpID(3-m,s,k)==0) .or. Best(m)==1) then
             GpID(m,s,k) = 0 
            else
              call CalcGpz(s,k, CandGP(Best(m),m), m,gLL)
            endif
          endif
        enddo        
      endif
    endif
    ALRtmp = 999.0D0          
    
    if (ANY(GpID(:,s,k)==0) .and. ANY(GpID(:,s,k)<0)) then   ! safety net
      do m=1,2
        if (GpID(m,s,k)<0) then
          if (skip(-GpID(m,s,k), m)) then
            GpID(:,s,k) = 0  ! can't be sure if grandmother or grandfather
          endif
        endif
      enddo
    endif
    
    if (ANY(GpID(:,s,k)==0) .and. ns(s,k)>1) then   ! double check that not offspring assigned as parent
      GPIN = GpID(:,s,k)
      GpID(:,s,k) = 0
      do m=1,2
        if (GPIN(m)/=0) then
          call CalcAgeLR(-s,k, GPIN(m),m, 0,1, .FALSE., ALRtmp(1))
          call CalcAgeLR(GPIN(m),m, -s,k, 0,1, .FALSE., ALRtmp(2))     
          if (ALRtmp(1)/=777 .and. (ALRtmp(2)==777 .or. &
           (ALRtmp(1)-ALRtmp(2)) > 3*TA)) then
            GpID(m,s,k) = GPIN(m)
          ! else drop - might be other way around
          endif
          call CalcCLL(s,k)
        endif
      enddo
      
      if (nS(s,k)==1 .and. ANY(GpID(:,s,k)>0)) then
        do m=1,2
          if (GPIN(m)>0) then
            call CheckRel(-s, k, GPIN(m), m, AgePhase>0, 4, LL)
            call BestRel(LL, 4, ix, LRG)
            if (ix == 4) then
              GpID(m,s,k) = GPIN(m)  
            endif
            call CalcCLL(s,k)
          endif
        enddo
      endif
    endif
    
    if (nS(s, k)==1 .and. ALL(GpID(:,s,k)==0)) then  ! single sib left; remove sibship  
      u = SibID(1,s, k)
      call RemoveSib(u, s, k)
      call DoMerge(0, s, k, OK)
      cycle
    endif
    
    call CalcCLl(s,k)
    do m=1,2  ! Shouldn't be necessary, as GGP among candGP?
      if (GpID(m,s,k)==curGP(m)) cycle
      if (GpID(m,s,k) < 0) then ! check if non-sampled non-dummy sib-of-dummy is GP
        GGP = GpID(:, -GpID(m,s,k), m)
        if (ALL(GGP <= 0))  cycle
        GPIN = GpID(:,s,k)
        do u=1,2
          if (GGP(u) > 0) then
            call CalcU(-s, k, GPIN(m), m, LLGGP(1,2))
            GpID(m,s,k) = 0
            call CalcU(-s, k, GPIN(m), m, LLGGP(1,1))
            call CalcU(-s, k, GGP(u), u, LLGGP(2,1))
            call AddGGP(GGP(u), s, k, LLGGP(2,2))
            if ((LLGGP(1,2)-LLGGP(1,1)) - (LLGGP(2,2)-LLGGP(2,1)) > TA) then
              GpID(m,s,k) = GPIN(m)  ! restore; else don't trust.
            else
              exit
            endif
            call CalcCLl(s,k)
          endif
        enddo
      else if (GpID(m,s,k)>0) then
        if (Sex(GpID(m,s,k))==3) then
          Sex(GpID(m,s,k)) = m  
        endif
      endif
    enddo
    
     if (nS(s, k)==1 .and. ALL(GpID(:,s,k)==0)) then  ! single sib left; remove sibship  
      u = SibID(1,s, k)
      call RemoveSib(u, s, k)
      call DoMerge(0, s, k, OK)
      cycle
    endif 
    !update ancestors & update likelihoods
    call CalcCLl(s,k)
    if (GPID(1,s,k)/=curGP(1) .or. GPID(2,s,k)/=curGP(2)) then
      call GetAncest(-s, k, AncR(s,k,:,:))
      do i=1, ns(s,k)
        u = SibID(i,s, k)
        call GetAncest(u,1,AncI(u,:,:))
        call calcLind(u)
        if (Parent(u, 3-k) < 0) then
          call CalcCLL(-Parent(u,3-k), 3-k)
          call calcLind(u)
          call CalcCLL(s,k)
        endif
      enddo
    endif
  enddo  ! s
enddo  ! k

end subroutine SibGrandparents

! #####################################################################

subroutine CalcGPz(SA, kA, B, kB, gLL) ! B or SB gp of sibship SA?
use Global
implicit none

integer, intent(IN) :: SA, kA, B, kB
double precision, intent(OUT) :: gLL                                    
integer :: m,n, CurGP(2), topX, x, mid(5), y, CY(4), kY(4), v
double precision :: LLA(2,7,7), dLL, PairLL(2,2,2), LLcp(3,2), LLU(3), &
 LLtmp(3), ALR(3), TopLL(2), LLAA(2,7,7)
logical :: ConPar(4,4), UseAge

curGP = GpID(:,SA,kA)
dLL = 999D0
PairLL = 999D0
LLA = 999D0
gLL = 999D0           
UseAge = (AgePhase>0)                     
 call CalcCLL(SA,kA)

if (B < 0) then
  n = kB
else if (B>0) then
  if (Sex(B)/=3) then
    n = Sex(B)
  else if (GpID(1,SA,kA)==0) then
    n = 1
  else !if (GpID(2,SA,kA)==0) then
    n = 2  ! TODO: check in both configs (as in POZ)
  endif
endif

if (curGP(1)==0 .and. curGP(2)==0) then
   call CheckRel(-SA, kA, B, kB, UseAge, 4, LLA(1,:,7)) 
  call BestRel(LLA(1,:,7), 4, topX, dLL)
  if (topX==4) then   !   .and. dLL>TA*nS(SA,kA)
    GpID(n,SA,kA) = B
  endif       
  if (LLA(1,4,7)<0 .and. LLA(1,7,7)<0) then
    gLL = LLA(1,4,7) - LLA(1,7,7)
  else
    gLL = 909D0
  endif
else
  ConPar = .FALSE.
  GpID(:,SA,kA) = 0
  call CalcCLL(SA,kA)
  if (ANY(curGP<0) .or. B<0) then
    do m=1,2
      if (curGP(m)==0) cycle
      call Connected(CurGP(m), m, -SA, kA, ConPar(4,m))
      call Connected(CurGP(m), m, B, n, ConPar(3,m)) 
      call Connected(curGP(3-m),3-m,CurGP(m), m, ConPar(2,1))
    enddo
    call Connected(B,n, -SA, kA, ConPar(4,3))
  endif

  LLU = 0 
  do m=1,2
    if (curGP(m)==0) cycle
    call CalcU(curGP(m),m, -SA,kA, LLtmp(1))
    LLU(m) = LLtmp(1) - CLL(SA,kA)
  enddo  
  call CalcU(B,n, -SA,kA, LLtmp(2)) 
  LLU(3) = LLtmp(2) - CLL(SA,kA)
  do m=1,2
    LLcp(m,m) = LLU(3-m) + LLU(3) 
  enddo
  LLcp(3,:) = LLU(1) + LLU(2)  ! B  
  
  CY = (/ CurGP(1), CurGP(2), B, -SA /)
  kY = (/ 1, 2, n, kA /)
  if (ANY(ConPar)) then 
    do m=1,2        
      do y=1,3  ! focal: CurGP, B 
        if (y/=m .and. y/=3) cycle           
        if (y==1) then
          call CalcU(CY(2),kY(2), CY(3),kY(3), LLCP(1,m))
        else if (y==2) then
          call CalcU(CY(1),kY(1), CY(3),kY(3), LLCP(2,m))    
        else if (y==3) then
          call CalcU(CY(1),kY(1), CY(2),kY(2), LLCP(3,m))
        endif
        do x=1,3
          if (ConPar(4,x) .and. x/=y) then
            call CalcU(CY(x),kY(x), CY(4),kY(4), LLtmp(1))
            call CalcU(CY(x),kY(x), 0,0, LLtmp(2))
            call CalcU(CY(4),kY(4), 0,0, LLtmp(3))
            LLCP(y,m) = LLCP(y,m) + (LLtmp(1)-LLtmp(2)-LLtmp(3))
          endif
          do v=1,2
            if (ConPar(x,v) .and. (x==y .or. v==y)) then
              call CalcU(CY(x),kY(x), curGP(v),v, LLtmp(1))
              call CalcU(CY(x),kY(x), 0,0, LLtmp(2))
              call CalcU(curGP(v),v, 0,0, LLtmp(3))
              LLCP(y,m) = LLCP(y,m) + (LLtmp(1)-LLtmp(2)-LLtmp(3))
            endif
          enddo
        enddo
      enddo
    enddo
  endif
  ! age based prior prob.
  ALR = 0
  do y=1,3
    call CalcAgeLR(-SA,kA, CY(y),kY(y), 0,1, .FALSE., ALR(y))
  enddo 
  
  mid = (/1,2,3,5,6/)
  GpID(:,SA,kA) = CurGP 
  call CalcCLL(SA,kA)
  dLL = 0
  LLtmp = 999D0
  do m=1,2 ! sex currently assigned GP
    if (CurGP(m)==0) cycle        
    call checkRel(-SA, kA, B, kB, UseAge, 4, LLA(m,:,4))   ! CurGP(m)=GP + A_7 
    GpID(m,SA,kA) = 0
    call CalcCLL(SA,kA)
    call checkRel(-SA, kA, B, kB, UseAge, 4, LLA(m,:,7))   ! A_7                                                            
    if (B > 0) then
      if (curGP(m)/=0 .and. Parent(B,m)==curGP(m)) then 
        LLA(m, 6, 7) = 888
        if (curGP(3-m)/=0 .and. Parent(B,3-m)==curGP(3-m)) &
          LLA(m, 5, 7) = 888  !FA
      endif
    else if (B < 0) then 
      if (curGP(m)/=0 .and. (GpID(m,-B,kB)==curGP(m) .or. &  
        ANY(SibID(1:nS(-B,kB),-B,kB)==curGP(m)))) then  
        LLA(m, 5:6, 7) = 888   
      endif
    endif
     if (curGP(m)>0) then
      if (Parent(CurGP(m),n)==B)  LLA(m, 6, 7) = 888   
    else if (curGP(m) < 0) then
      if (GpID(n,-CurGP(m),m)==B)  LLA(m, 6, 7) = 888
    endif

    call checkRel(-SA,kA, CurGP(m), m, UseAge, 4, LLA(m,7,:))  ! CurGP(m)_7
    if (curGP(m) < 0) then                  
      if (m/=n .and. ANY(Parent(SibID(1:nS(-curGP(m),m),&
        -curGP(m),m),3-m)==B)) then
        call PairUA(-SA, CurGP(m), kA, m, LLA(m,7,4))   
      else if (GpID(kB,-curGP(m),m)==B) then
        LLA(m, 7, 5:6) = 888                    
      endif
    endif
    
    if (LLA(m,4,4)<0 .or. LLA(m,4,7)<0) then  ! Else not possible
      GpID(n,SA,kA) = B
      call CalcCLL(SA,kA)
      call checkRel(-SA,kA, CurGP(m),m, UseAge, 4, LLA(m,4,:)) 
    endif

    WHERE (LLA(m,mid,4)<0) LLA(m,mid,4) =LLA(m,mid,4) +LLcp(3,m) +ALR(m)
    WHERE (LLA(m,mid,7)<0) LLA(m,mid,7) = LLA(m,mid,7) + LLcp(3,m)
    WHERE (LLA(m,4,:)<0) LLA(m,4,:) = LLA(m,4,:) + LLcp(m,m) + ALR(3)
    WHERE (LLA(m,7,:)<0) LLA(m,7,:) = LLA(m,7,:) + LLcp(m,m)
    
    if (n==m) then
      if (B>0) then
        if (Sex(B) == 3) then
          LLA(m,4,4) = dLL  ! from checkAdd(B, SA,..)
        else
          LLA(m,4,4) = 777
        endif
      else   
        LLA(m,4,4) = 777 ! cannot have 2 same-sex GPs
      endif
    endif
    
    GpID(:,SA,kA) = CurGP  ! restore
    call CalcCLL(SA,kA)  
  enddo    
  
  LLAA = LLA
  PairLL(1,1,2) = MaxLL((/LLA(3-n,4,4), LLA(n,4,:)/))  ! B + curGP(3-n)
  LLAA(3-n,4,4) = 555
  PairLL(1,2,2) = MaxLL(LLAA(3-n,4,:))   ! only B
  PairLL(2,1,1) = MaxLL(RESHAPE(LLA(:,:,4), (/2*7/)))  ! both curGP's
  LLAA(:,4,:) = 555
  LLAA(:,:,4) = 555
  PairLL(2,1,2) = MaxLL(RESHAPE(LLAA(n,:,:), (/7*7/)))   ! only curGP(3-n)
  PairLL(2,2,1) = MaxLL(RESHAPE(LLAA(3-n,:,:), (/7*7/))) ! only curGP(n)
  PairLL(2,2,2) = 555  ! TODO
  
  TopLL(1) = MaxLL(RESHAPE(PairLL, (/8/)))
  TopLL(2) = MAXVAL(PairLL, MASK=PairLL < TopLL(1))
  GpID(:,SA,kA) = 0
  if (PairLL(1,1,2) == TopLL(1) .and. TopLL(1) - TopLL(2) > TA) then
    GpID(n,SA,kA) = B
    GPID(3-n, SA, kA) = CurGP(3-n)
    gLL = LLA(3-n,4,4)
  else if (PairLL(2,1,1) == TopLL(1) .and. TopLL(1) - TopLL(2) >TA) then
    GpID(:,SA,kA) = curGP
    gLL = MaxLL(LLA(:,7,4))
  else
    PairLL(1,1,2) = 555
    PairLL(2,1,1) = 555
    TopLL(1) = MaxLL(RESHAPE(PairLL, (/8/)))
    TopLL(2) = MAXVAL(PairLL, MASK=PairLL < TopLL(1))
    if ( TopLL(1) - TopLL(2) > TA) then
      if (PairLL(1,2,2) == TopLL(1)) then
        GpID(n,SA,kA) = B  
        gLL = LLA(3-n,4,7)
      else if (PairLL(2,1,2) == TopLL(1)) then
        GpID(3-n,SA,kA) = curGP(3-n)
        gLL = LLA(n,7,7)
      else if (PairLL(2,1,2) == TopLL(1)) then
        GpID(n,SA,kA) = curGP(n)
        gLL = LLA(3-n,7,7)
      endif
    endif
  endif
  if (gLL>0) then
    gLL = 909D0
  else if (LLA(3-n,7,7) < 0) then
    gLL = gLL - LLA(3-n,7,7) !- LLU(n) 
  else if (LLA(n,7,7) < 0) then
    gLL = gLL - LLA(n,7,7) !- LLU(3-n) 
  else
    gLL = 909D0 
  endif
endif

 call CalcCLL(SA, kA)
do n=1, nS(SA,kA)
  if (Parent(SibID(n,sA,kA),3-kA) < 0) then
    call CalcCLL(-Parent(SibID(n,sA,kA),3-kA), 3-kA)
  endif         
  call calcLind(SibID(n, SA, kA))
enddo
 call CalcCLL(SA, kA)
     
end subroutine CalcGPz

! #####################################################################

subroutine GGpairs(ExtraAge)  ! find & assign grandparents of singletons
use Global
implicit none

integer, intent(IN) :: ExtraAge                               
integer :: i, j, k, x, TopX, nCG(2,2), CandG(2,2,50), m, n, curGP(2), y, s
double precision :: LLa(7), LLg(7), dLL, LRS, gLL, LL(7)
logical :: IsPair, OK

do i=1, nInd
  if (Error/=0) return
  if (MODULO(i,500)==0) call rchkusr()
  if (Parent(i,1)==0 .and. Parent(i,2)==0) cycle
  if (Parent(i,1)/=0 .and. Parent(i,2)/=0) cycle
  if (BY(i) < 0) cycle
  nCG = 0
  CandG = 0
  do k=1,2
    if (Parent(i,k)/=0) cycle
    do j=1, nInd
      if (ANY(nCG(k,:)>=50)) cycle
      if (ALL(Parent(j,:)==0 .and. AgePhase==0)) cycle
      if (AgeDiff(i,j) <= 0)  cycle
      if (AgeDiff(i,j) == 999D0) cycle  ! TODO?
      if (Sex(j)==3) cycle
      x = 0
      if (k==1) then
        if (Sex(j)==1) then
          x = 3 ! mat. grandmother
        else 
          x = 5  ! mat. grandfather
        endif
      else if (k==2) then
        if (Sex(j)==2) then
          x = 4  ! pat. grandfather
        else 
          x = 5  ! pat. grandmother
        endif
      endif
      if (AgeDiff(i,j)/=999) then
        if (LOG10(AgePriorM(AgeDiff(i, j)+1, x)) < TF) cycle  ! unlikely
      endif
      call PairQHS(i, j, LRS) 
      if (LRS < TF)  cycle
      IsPair = .TRUE.
      call CheckPair(i,j,k, 4, Llg, LLa)                                  
      do n=1,2
        if (AgePhase==0 .and. n==1)  cycle
        if (n==1)  LL = LLa
        if (n==2)  LL = LLg                   
        call BestRel(LL, 4, topX, dLL)
        if (topX == 4 .and. dLL>TA) then
          IsPair = .TRUE.
        else if (LL(4)- MaxLL(LL((/1,2,6,7/))) > TA .and. n==2) then
          IsPair = .TRUE.
        else
          IsPair = .FALSE.
          exit
        endif
      enddo
      if (.not. IsPair) cycle
      nCG(k, Sex(j)) = nCG(k, Sex(j)) +1
      CandG(k, Sex(j), nCG(k,Sex(j))) = j
    enddo
  enddo
      
  if (SUM(nCG)>=1) then
    do k=1,2
      if (SUM(nCG(k,:))>=1) then
        nC(k) = nC(k) + 1
        nS(nC(k),k) = 1
        SibID(1, nC(k), k) = i
        Parent(i, k) = -nC(k)
        s = nC(k)
        do x=1, MAXVAL(nCG(k,:))
          do m=1,2
            if (candG(k,m,x) == 0)  cycle
            call calcGPZ(s, k, CandG(k,m,x),m,gLL)   !!!!
          enddo
        enddo
        if (ALL(GpID(:,s,k)==0) .and. ALL(nCG(k,:)>=1)) then  ! check all pairs of cand gps
          do x=1, nCG(k,1)
            do y=1, nCG(k,2) 
              if (GpID(1,s,k) == CandG(k,1,x)) then
                call CalcGPz(s, k, CandG(k,2,y), 2, gLL)
              else if (GpID(2,s,k) == CandG(k,2,y)) then
                call CalcGPz(s, k, CandG(k,1,x), 1, gLL)
              else
                curGP = GpID(:,s,k)                                    
                GpID(1,s,k) = CandG(k,1,x)
                call CalcGPz(s, k, CandG(k,2,y), 2, gLL)
                if (curGP(1)/=0 .and. GpID(1,s,k) == 0) then  ! restore GP1?
                  if (GpID(2,s,k) == curGP(2)) then
                    GPID(1,s,k) = curGP(1)  
                  else if (curGP(2)/=0 .and. GpID(2,s,k)==0) then
                    GPID(:,s,k) = curGP
                  else
                    call CalcGPz(s, k, curGP(1), 1, gLL)
                  endif
                endif                                                         
              endif
            enddo
          enddo
        endif
        
        if (ANY(GpID(:,s,k)==0)) then  ! double check doesn't depend only on age
          curGP = GpID(:,s,k)
          GpID(:,s,k) = 0
          call CalcCLL(s, k)
          do m=1,2
            if (curGP(m)==0)  cycle
            call CheckPair(i,curGP(m),k, 4, Llg, LLa)
            call BestRel(LLg, 4, topX, dLL)
            if (topX == 4) then
              GpID(m,s,k) = curGP(m)
            else 
              call BestRel(LLa, 4, topX, dLL)
              if (topX == 4 .and. (dLL > 3*TA .or. ExtraAge>0)) then
                GpID(m,s,k) = curGP(m) 
              endif
            endif
          enddo
          call CalcCLL(s, k)
        endif

        if (ALL(GpID(:,s,k)==0)) then
          call RemoveSib(i, s, k)
          call DoMerge(0, s, k, OK)
        endif
      endif
    enddo
  endif 
enddo

end subroutine GGpairs

! #####################################################################

subroutine Qadd(A, SB, kB, LR)
use Global
implicit none

integer, intent(IN) :: A, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x
double precision :: PrL(nSnp), PrX(3)

PrL = 0D0
do l=1,nSnp
  if (Genos(l,A)==-9) cycle
  do x=1,3
    PrX(x) = OKAP(Genos(l,A), x, l) * DumP(x,l,SB,kB) / AHWE(x,l)
  enddo   ! simple LL identical for HS and GP
  PrL(l) = LOG10(SUM(PrX))
enddo
LR = SUM(PrL)

end subroutine Qadd

! #####################################################################

subroutine QGP(A, kA, SB, kB, LR)  ! A indiv or dummy
use Global
implicit none

integer, intent(IN) :: A, kA, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x, sib1
double precision :: PrL(nSnp), PrX(3), PrXY(3,3), LL(2), PrA(3)

if (ns(SB,kB)==1 .and. A>0) then
  sib1 = SibID(1,SB,kB)
  call CalcU(sib1,kB,A,kA, LL(1))
  call PairGP(Sib1, A, kA, 4, LL(2))
  LR = LL(2) - LL(1)
else
  PrL = 0D0
  do l=1,nSnp
    call ParProb(l, A, kA, 0, 0, PrA)
    do x=1,3
      PrX(x) = XPr(1,x,l,SB,kB) * AHWE(x,l)          
      PrXY(x,:) = XPr(1,x,l,SB,kB) * AKAP(x,:,l) * PrA    
    enddo
    PrL(l) = LOG10(SUM(PrXY)) - LOG10(SUM(PrX))
  enddo
  LR = SUM(PrL)
endif

end subroutine QGP

! #####################################################################

subroutine QPO(A, SB, kB, LR)  ! A replaces dummy SB?
use Global
implicit none

integer, intent(IN) :: A, SB, kB
double precision, intent(OUT) :: LR
integer :: l, x, sib1
double precision :: PrL(nSnp), PrX(3,2), LL(2), PrA(3)

if (ns(SB,kB)==1) then
  sib1 = SibID(1,SB,kB)
  call CalcU(sib1,kB,A,kB, LL(1))
  call PairPO(sib1, A, kB, 1, LL(2))
  LR = LL(2) - LL(1)
else
  PrL = 0D0
  do l=1,nSnp
    call ParProb(l, A, kB, 0, 0, PrA)
    do x=1,3
      PrX(x,1) = XPr(1,x,l,SB,kB) * XPr(2,x,l,SB,kB)
      PrX(x,2) = XPr(1,x,l,SB,kB) * PrA(x)
    enddo
    PrL(l) = LOG10(SUM(PrX(:,2))) - LOG10(SUM(PrX(:,1)))
  enddo
  LR = SUM(PrL)
endif

end subroutine QPO

! #####################################################################

subroutine CheckRel(A, kA, B, kB, InclAge, focalIN, LL)  ! , LLg
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB, focalIN
logical, intent(IN) :: InclAge
double precision, intent(OUT) :: LL(7)  ! LLg(7), 
double precision :: LLg(7)
logical:: FSJ  !do separately?
integer :: k, focal

focal = focalIN               
FSJ = .FALSE.
if (A > 0 .and. B > 0) then
  if (kA /= 0) then
    k = kA
  else if (kB /= 0) then
    k = kB
  endif
  call CheckPair(A, B, k, focal, LLg, LL)  
  if (InclAge) then
    LL = LL
  else 
    LL = LLg   
  endif
else if (A > 0 .and. B < 0) then
  if (focal==1)  focal =  3  ! -B parent of A -> B's HS of A                                                            
  call CheckAdd(A, -B, kB, InclAge, LL, focal) 
else if (A < 0 .and. B > 0) then
  call CheckAdd(B, -A, kA, InclAge, LL, focal)
else if (A < 0 .and. B < 0) then
  call CheckMerge(-A, -B, kA, kB, InclAge, LL, focal, FSJ)
endif

end subroutine CheckRel

! #####################################################################

subroutine CheckAdd(A, SB, k, InclAge, LL, focal)
use Global
implicit none

integer, intent(IN) :: A, SB, k, focal
logical, intent(IN) :: InclAge
double precision, intent(OUT) :: LL(7)
double precision :: LLg(7),  LLtmp(2,3), ALR(7), LLz(6), LRHS, LHH(3), &
 ALRx(3,2), LLM(3), LLPX(2,2), LLp(7), LLx(3), LLC, dx(maxSibSize), ALRq, &  
 dLL(2), LLPg(7), LLPO, LLFH(3), ALRz(6), LLi(ns(SB,k), 3), LLHH(4), LLxp(7)                                                                     
integer :: x, y, Par, MaybeOpp, i, ParTmp(2), npt, fsi, ix, topX(2), &
  AncA(2,mxA), m, Bi, AncBi(2, mxA), AncB(2,mxA), sib1, curpar(2)                                                                               
logical :: maybeFA                       

LL = 999D0
LLg = 999D0
LLz = 999D0
LRHS = 999D0
ALR = 999D0
LLtmp = 999D0
ALRz = 999D0                        
 call Qadd(A, SB, k, LRHS)  ! 2nd degree relatives vs unrelated
if (LRHS < TF*nS(SB,k) .and. (focal/=4 .and. focal/=7)) return   
if (focal==1) then
  if (Sex(A)/=3 .and. Sex(A)/=k .and. .not. hermaphrodites) then
    LL(1) = 777
    return
  endif
  call CalcAgeLR(-SB,k, A,k, 0,-1, .TRUE., ALR(1))
  if (ALR(1)==777 .or. ALR(1)<5*TF) then
    LL(1) = 777
    return
  endif
else if (focal==2 .or. focal==3) then
  call CalcAgeLR(A,Sex(A), -SB,k, k,1, .TRUE., ALR(3))
  if (ALR(3)==777 .or. ALR(3)<5*TF) then
    LL(2:3) = 777
    return
  endif
endif

if (ANY(Parent(A,:)/=0)) then
  call GetAncest(A, k, AncA)
  if (ANY(AncA(k, :) == -SB)) then
    LL(1) = 777
    LL(4) = 777
    if (focal==1 .or. focal==4)  return
  endif
endif

 call CalcU(A,k, -SB, k, LLg(7))   ! unrelated
LL(7) = LLg(7)

fsi=0
if (focal <4) then             
  if (focal==1)  call AddParent(A, SB, k, LLg(1))
  if (focal==2)  call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
  if (focal==3)  call AddSib(A, SB, k, LLg(3))
endif 
do x=1,3
  if (focal==x) then              
    if ((LLg(focal) > 0 .or. LLg(focal) - LL(7) < TA)) then
      if (focal < 3) then
        return
      else
        call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
        if ((LLg(2) > 0 .and. LLg(3)>0) .or. &
         ((LLg(2) - LL(7) < TA) .and. (LLg(3) - LL(7) < TA))) return
      endif
    endif
  endif
enddo    
 !=======
if (LL(1)/=777 .and. LLg(1)==999) then
  call AddParent(A, SB, k, LLg(1))  ! A parent of SB
endif
if (InclAge .and. ALR(1)==999)  call CalcAgeLR(-SB,k, A,k, 0,-1, .TRUE., ALR(1))

if (LLg(2)==999)  call AddFS(A, SB, k,0,k, LLg(2), fsi, dx)
if (Complx>0 .and. LLg(3)==999)  call AddSib(A, SB, k, LLg(3))
 call CalcAgeLR(A,Sex(A), -SB,k, 0,1, .TRUE., ALR(3))  ! SB parent of A
 call CalcAgeLR(A,Sex(A), fsi,k, 0,2, .TRUE., ALR(2))

do x=1,3
  if (LLg(x) < 0 .and. ALR(x) /= 777 .and. (InclAge .or. ALR(x)>5*TF)) then
    if (InclAge)  LL(x) = LLg(x) + ALR(x)
    if (.not. InclAge)  LL(x) = LLg(x)
  else
    LL(x) = 777
  endif
enddo

if (nAgeClasses>2 .and. LL(4)/=777) then
  call AddGP(A, SB, k, LLg(4))
  if (Sex(A)<3) then
    call CalcAgeLR(-SB,k, A,Sex(A), 0,1, .TRUE., ALR(4))  ! A parent of SB
  else
    call CalcAgeLR(-SB,k, A,3-k, 0,1, .TRUE., ALR(4))
  endif
endif
if (LLg(4) < 0) then 
  if (ALR(4) /= 777 .and. (InclAge .or. ALR(4)>5*TF)) then
    if (InclAge)  LL(4) = LLg(4) + ALR(4)
    if (.not. InclAge)  LL(4) = LLg(4)
  else
    LL(4) = 777
  endif
else
  LL(4) = LLg(4)
endif

! FAU
LLC = 999D0 
Par = 0
maybeFA = .TRUE.                
call GetAncest(-SB, k, AncB)
call CalcAgeLR(-SB,k, A,3-k, 0,2, .TRUE., ALR(5))
call CalcAgeLR(A,3-k, -SB,k, k,4, .TRUE., ALR(6))   !B's FA of A?
if (ALR(5)/=777 .and. (InclAge .or. ALR(5)>5*TF)) then
   if (.not. (focal==4 .and. ALL(Parent(A,:)/=0))) then
    call pairUA(-SB, A, k, 3, LLtmp(1,3))    ! A FS with SB?
  endif
endif                                                                                        
if (ALR(6)/=777 .and. (InclAge .or. ALR(6)>5*TF) .and. .not. ANY(AncB == A)) then  
  call getFSpar(SB, k, .TRUE., Par)  ! TODO non-strict: needs check if Parent(B1,3-k)==Par                                                                                        
  if (Par /= 0) then
    call CalcAgeLR(A,3-k, Par,3-k, 0,1, .TRUE., ALRq) 
    if (ALRq < 2*TF .or. ALRq==777)  maybeFA = .FALSE. 
    if (maybeFA) then
      call getAncest(Par, 3-k, AncBi)
      if (ANY(AncBi == A))   maybeFA = .FALSE. 
    endif
  else
    maybeFA = .FALSE.
  endif
  if (maybeFA) then
    call addFA(A, SB, k, LLtmp(2,3))
  endif
endif

LLg(5) = MaxLL(LLtmp(:,3))
if (InclAge) then
  if (ALL(LLtmp(:,3)<0)) then
    LL(5) = MaxLL((/LLtmp(1,3)+ALR(5), LLtmp(2,3)+ALR(6)/))
  else if (LLtmp(1,3) < 0) then
    LL(5) = LLtmp(1,3)+ALR(5)
  else if (LLtmp(2,3) < 0) then
    LL(5) = LLtmp(2,3)+ALR(6)
  else
    LL(5) = LLg(5)    ! >0
  endif
endif

if (complx==2 .and. (focal==2 .or. focal==3 .or. focal==7) .and. LL(2)<0 .and. &
 Parent(A,3-k)==Par .and. (MaxLL(LLtmp(:,3)) - MaxLL(LL(2:3)) > -TA)) then
    call FSHC(A, -SB, k, LLC)
    if (LLC >LLg(2) .and. LLC<0 .and. ALR(2)/=777) then
      LL(2) = LLC + ALR(2)
    endif
endif
  
! HAU
ALRx = 0
do x=1,2
  call pairUA(A, -SB, x, k, LLtmp(1,x))  
  if (Parent(A,x)>0) then
    call CalcAgeLR(Parent(A,x),x, -SB,k, 0,1, .TRUE., ALRx(1,x))
  else
    call CalcAgeLR(A,Sex(A), -SB,k, x,4, .TRUE., ALRx(1,x))
  endif
  if (.not. (focal==4 .and. Parent(A,x)/=0)) then                                       
    call pairUA(-SB, A, k, x, LLtmp(2,x))
  endif                                                       
  call CalcAgeLR(-SB,k, A,Sex(A), x,3, .TRUE., ALRx(2,x))
enddo
if (Complx>0)  LLg(6) = MaxLL(RESHAPE(LLtmp(:,1:2), (/2*2/) ))
do x=1,2
  do y=1,2
    if (LLtmp(y,x) < 0 .and. ALRx(y,x) /= 777 .and. &
      (InclAge .or. ALRx(y,x)>5*TF)) then
        if (InclAge)  LLtmp(y,x) = LLtmp(y,x) + ALRx(y,x)
    else
      LLtmp(y,x) = 777
    endif
  enddo
enddo
if (Complx>0)  LL(6) = MaxLL(RESHAPE(LLtmp(:,1:2), (/2*2/) ))    
if (LLg(5)<0 .and. ALRx(2,1)/=777 .and. ALRx(2,2)/=777) then
  if (InclAge)  LL(5) = LLg(5) + ALRx(2,1) + ALRx(2,2)  ! A FS of SB
  if (.not. InclAge)  LL(5) = LLg(5)
else if (LLg(5)<0 .and. ALRx(1,1)/=777 .and. ALRx(1,2)/=777) then
  if (InclAge)  LL(5) = LLg(5) + ALRx(1,1) + ALRx(1,2)  ! one or all in SB FA of A
  if (.not. InclAge)  LL(5) = LLg(5)
else
  LL(5) = 777
endif                                                              

if (focal>0 .and. focal<8) then                               
if ((LL(focal)<0 .and. LL(focal)>=LL(7)) .or. focal==4 .or. LL(6)>0) then                                                           
  if (nAgeClasses>3 .and. ALR(4)/=777 .and. LLg(4)<0) then 
    call CalcAgeLR(-SB,k, A,Sex(A), 1,4, .TRUE., ALRz(1))
    call CalcAgeLR(-SB,k, A,Sex(A), 2,4, .TRUE., ALRz(2))
    if ((ALRz(1)/=777 .and. ALRz(1)>3*TF) .or. (ALRz(2)/=777 .and. ALRz(2)>3*TF)) then
      call AddGGP(A, SB, k, LLz(1))
    endif
  endif               
  if (nS(SB,k)>0) then
    call CalcAgeLR(A,k, -SB,k, 0,5, .TRUE., ALRx(3,1))
    if (ALRx(3,1)==777 .or. ALRx(3,1)<5*TF) then
      LLz(2:3) = 777
    else
      call ParentHFS(A, 0,1, SB, k,3, LLz(2))
      call ParentHFS(A, 0,2, SB, k,3, LLz(3))
    endif
  endif
  if (Complx==2) then
    do x=1,2   ! as checkmerge: full great-uncle  (2x 1/4)
      call CalcAgeLR(-SB,k, A,Sex(A), 1,5, .TRUE., ALRz(3+x))
      if (ALRz(3+x) /= 777) then
        if (GpID(x,SB,k) <0) then 
          call PairUA(GpID(x,SB,k), A, x, 3, LLz(3+x))
          if (LLz(3+x) < 0) then
            LLz(3+x) = LLz(3+x) - CLL(-GpID(x,SB,k), x) + CLL(SB,k)  
          endif
        else if (GpID(x,SB,k)==0) then   ! else cond. indep.                                                          
          call addGAU(A, SB, k, x, LLz(3+x))    
        endif
      endif
    enddo
  endif
  if (ALR(3)/=777 .and. ns(SB,k)>0) then
    sib1 = SibID(1,SB,k)
    call PairCC(A, sib1, k, LLz(6))
    if (LLz(6) < 0) then
      call CalcU(A, k, sib1, k, LLxp(1))
      LLz(6) = LLz(6) -LLxp(1) + CLL(SB,k)
    endif
    ALRz(6) = ALR(3)  ! no ALR for cousins yet
  endif
  LLg(6) = MaxLL((/LLg(6), LLz/))
  if (InclAge) then
    ALRz(1) = MaxLL(ALRz(1:2))
    ALRz(2:3) = ALRx(3,1)
    do x=1,5
      if (LLz(x) < 0 .and. ALRz(x) /= 777) then
        LLz(x) = LLz(x) + ALRz(x)
      else
        LLz(x) = 777
      endif
    enddo
  endif
  LL(6) = MaxLL((/LL(6), LLz/))
endif
endif     
      
 LLi = 999D0 
if ((MaxLL(LL)==LL(3) .or. MaxLL(LL)==LL(2)) .and. ANY(Parent(A,:)<=0) .and. &
  ANY(AgeDiff(A, SibID(1:ns(SB,k),SB,k)) > 1) .and. focal/=7) then
  do i=1, ns(SB,k)  ! an B GP of A
    if (AgeDiff(A, SibID(i,SB,k)) <=0) cycle
    if (AgeDiff(A, SibID(i,SB,k))/=999) then
      if (ALL(LOG10(AgePriorM(AgeDiff(A, SibID(i,SB,k))+1, 3:5)) < TF)) cycle
    endif
    do x=1,2
      Sib1 = SibID(i,SB,k)
      call PairGP(A, Sib1, x, 4, LLi(i,x))
    enddo
    call CalcU(A,k, Sib1,k, LLi(i,3))
    if (MaxLL(LLi(i,1:2)) < 0 .and. (MaxLL(LLi(i,1:2)) - LLi(i,3)) > TA) then
      if ((MaxLL(LLi(i,1:2)) - LLi(i,3) + LLg(7)) > LL(6)) then 
        LLg(6) = MaxLL(LLi(i,1:2)) - LLi(i,3) + LLg(7)
        LL(6) = LLg(6)   ! TODO: ageprior (GP)
      endif
    endif
  enddo
endif

LLM = 999D0    
LLp = 999D0
LLFH = 999D0       
LLPX = 999D0                          
MaybeOpp = 0  
Par = 0                 
call getFSpar(SB, k, .FALSE., Par)                                                             
if (complx>0 .and. (MaxLL(LL)==LL(3) .or. MaxLL(LL)==LL(2)) .and. &
  (focal==2 .or. focal==3) .and. Parent(A,3-k)==0) then
  if (MaxLL(LL)==LL(2)) then
    call CalcU(A, k, fsi, k, LLM(1))
    call PairHalfSib(A, fsi, 3-k, LLM(2)) 
    if ((LLM(2) - LLM(1)) - (LLg(2) - LLg(7)) > TA) then
      LL(2) = 222   ! more likely to be HS via 3-k
    endif
  endif

  MaybeOpp = 1
  if (Par > 0) then
    if (par/=A) then
      call CheckPair(A, par, k, 1, LLxp, LLp)    !! DANGER !!!
      if (LLp(1)<0 .and. (LLp(1) - MaxLL(LLp)) > TF) then  
        LL(2:3) = 222  ! par plausible parent of A
      endif
    else if (par==A) then ! e.g. when BY of A unknown
      do i=1, nS(SB,k)
        if (Parent(SibID(i,SB,k), 3-k)==0) then   ! todo: <=0
          call CheckPair(SibID(i,SB,k), A, 3-k, 1, LLp, LLPg)  !! DANGER !!!
          if (LLp(1)<0 .and. LLp(1) - MaxLL(LLp) > TF) then
            LL(7) = LL(7) + LLp(1) - LLp(7)
          endif
        endif
      enddo
    endif
  else if (Par==0) then
    if (ANY(Parent(SibID(1:nS(SB,k), SB,k),3-k)>0)) then
      MaybeOpp = 0
    else   ! check if opp. parent possibly to be merged
      npt = 0  ! number of unique opposite-sex dummy parents
      ParTmp = 0
      do i=1, nS(SB,k)
        if (Parent(SibID(i,SB,k), 3-k)<0 .and. &
          .not. ANY(ParTmp == Parent(SibID(i,SB,k), 3-k))) then
          npt = npt + 1
          if (npt > 2) then
            MaybeOpp = 0
            exit
          else
            ParTmp(npt) = Parent(SibID(i,SB,k), 3-k)
          endif
        endif
      enddo
      if (MaybeOpp == 1 .and. npt==2) then
        call CalcU(ParTmp(1), 3-k, ParTmp(2), 3-k, LLM(1))
        call MergeSibs(-ParTmp(1), -ParTmp(2), 3-k, LLM(2))
        if ((LLM(2) - LLM(1)) < TF*nS(SB,k))  MaybeOpp = 0
      endif
    endif
    if (ANY(Parent(SibID(1:ns(SB,k), SB, k), 3-k) < 0) .and. MaybeOpp==0) then
      do i=1, nS(SB,k)
        if (Parent(SibID(i,SB,k), 3-k)<0) then
          call CalcU(A, 3-k, Parent(SibID(i,SB,k), 3-k), 3-k, LLM(1))
          call AddSib(A, -Parent(SibID(i,SB,k), 3-k), 3-k, LLM(2))
          if ((LLM(2) - LLM(1)) - (LLg(2) - LLg(7)) > TA*nS(SB,k)) then
            LL(2) = 222   ! more likely to be added to opposing sibship only. 
          endif
          if ((LLM(2) - LLM(1)) - (LLg(3) - LLg(7)) > TA*nS(SB,k)) then
            LL(3) = 222  
          endif
          if (LL(2)==222 .and. LL(3)==222)  exit 
        endif
      enddo
    endif     
  else if (Par < 0) then
    call Qadd(A, -Par, 3-k, LLM(1))  ! 2nd degree relatives vs unrelated    
    if (LLM(1) < TF*nS(-Par,3-k))  MaybeOpp = 0
    call CalcAgeLR(A,Sex(A), Par,3-k, 0,1, .TRUE., ALRq)
    if (ALRq==777)  MaybeOpp = 0                                        
  endif
  if (MaybeOpp == 1 .and. Par < 0) then
    LLM = 999D0
    if (Par < 0) then  ! may have more/fewer sibs
      call AddFS(A, -Par, 3-k,0,3-k, LLM(1), ix, dx)
      call AddSib(A, -Par, 3-k, LLM(2))
      call CalcU(A, 3-k, Par, 3-k, LLM(3))
    else if (Par == 0  .and. nS(SB,k)>0) then
      sib1 =  SibID(1, SB, k)
      call PairFullSib(A, sib1, LLM(1))  
      call PairHalfSib(A, sib1, 3-k, LLM(2))
      call CalcU(A, 3-k, sib1, 3-k, LLM(3))
    endif
    if (LLM(2) < 0) then
      if (Par < 0 .and. complx>0) then
        if ((LLM(2) - LLM(3)) - (LLg(3) - LLg(7)) > TA*MAX(nS(SB,k),nS(-par,3-k))) then
          LL(2:3) = 222  
        endif
      endif
      if (LLM(1) < 0 .and. ((LLM(1) - LLM(2)) > 2*TA .or. Complx==0)) then
        if (Par<0 .and. Complx==2) then  ! HS + parents FS/PO?
          curPar = Parent(A,:)                                        
          call NewPar(A,k,-SB)
          call PairUA(A, Par, 3-k, 3-k, LLPX(1,1))  ! HS + HA
          call ParentHFS(A, 0,3-k,-Par, 3-k,3, LLPX(1,2))  ! HS + FC
          call NewPar(A, k, curPar(k))
          call NewPar(A, 3-k, par)
          call PairUA(A, -SB, k, k, LLPX(2,1))  ! HA + HS
          call ParentHFS(A, 0,k,SB, k,3, LLPX(2,2))  ! FC + HS
          call NewPar(A, 3-k, curPar(3-k))
          if ((LLg(2) - LLg(7)) - (MaxLL(LLPX(1,:)) - LLM(3)) < TA) then 
            LL(2) = 222
          endif
          if ((MaxLL(LLPX(1,:)) - LLM(3)) > (LLg(3) - LLg(7)) .and. &
           MaxLL(LLPX(1,:)) - MaxLL(LLPX(2,:)) > TA) then
            LLg(3) = MaxLL(LLPX(1,:)) - LLM(3) + LLg(7)
            LL(3) = LLg(3) + ALR(3)
          else if (((MaxLL(LLPX(2,:)) - LLM(3)) - (LLg(3) - LLg(7))) > TA) then  ! MAX(nS(SB,k),nS(-par,3-k))
            LL(3) = 222
          endif        
        else
          LL(2) = LL(2)
        endif
      else if (LLM(3)<0 .and. LLM(2) -LLM(3) >2*TA .and. complx>0) then
        LL(2:3) = 222  ! as likely to be added to opp. parent
      endif
      endif
  endif
  if (Par <= 0 .and. LL(2)<0 .and. Complx==2 .and. ns(SB,k)>0) then  
    sib1 = SibID(1,SB,k)
    call calcU(A,k,sib1, k, LLFH(1))
    call pairFAHA(A, sib1, .TRUE., LLFH(2))
    call pairFAHA(sib1, A, .TRUE., LLFH(3))
    WHERE(LLFH(2:3)<0)  LLFH(2:3) = LLFH(2:3) - LLFH(1) + LLg(7) 
    if (ANY(LLFH(2:3)<0) .and. MaxLL(LLFH(2:3)) > LLg(5)) then
      LLg(5) = MaxLL(LLFH(2:3))
      if (InclAge .and. ANY(ALRx < 777)) then
        LL(5) = LLg(5) + MaxLL((/ALRx(:,1), ALRx(:,2)/))
      else
        LL(5) = LLg(5)
      endif          
    endif
  endif
endif

LLHH = 999D0
if (MaxLL(LL)==LL(2) .and. (focal==2 .or. focal==3) .and. complx==2 .and. &
  Parent(A,3-k)==0) then
  do x=1,3
    call PairHSHA(A, fsi, k, x, LLHH(x))
  enddo 
  call CalcU(A, k, fsi, k, LLHH(4))
  if (MaxLL(LLHH) <0 .and. (LLg(2) - LLg(7)) - (MaxLL(LLHH(1:3)) - LLHH(4)) < TA) then
    LLg(2) = 222
    LL(2) = 222
  endif
endif

LLX = 999D0
if ((MaxLL(LL)==LL(3) .or. MaxLL(LL)==LL(2)) .and. &
  (focal==2 .or. focal==3) .and. ALL(Parent(A,:)==0)) then  
  do i=1, nS(SB,k)
    if (AgeDiff(A,SibID(i,SB,k))==999) then  ! else already assigned
      call calcU(A,k,SibID(i,SB,k), k, LLX(1))
      call PairPO(A,SibID(i,SB,k), k, 0, LLX(2))
      if (LLX(2)<0 .and. LLX(2) - LLX(1) > TA) then
        call NewPar(A, k, SibID(i,SB,k)) 
        call CalcU(A,k, -SB,k, LLX(3))
        call NewPar(A, k, 0)
        call CalcLind(A)
        if (LLX(3) > MaxLL(LL)) then
          LL(7) = LLX(3)
          exit
        endif
      endif
    endif
  enddo
endif

LHH = 999D0
if (complx==2 .and. nAgeClasses>1 .and. (focal==2 .or. focal==3 .or. focal==7) &
 .and. LL(3)<0 .and. LL(3)>=LL(7)) then                                      
  call AddSibInbr(A, SB, k, LHH)  ! 1: Par(Parent(A,3-k),k)=SB, 2: Parent(A,3-k)=GpID(3-k,SB,k), 3: as 1, A FS of B's (PA == DB)
  if (MaxLL(LHH(1:2)) - LLg(3) > 2*TA .and. MaxLL(LHH(1:2))<0) then
    if (InclAge .and. ALR(3)<777) then
      LL(3) = MaxLL(LHH(1:2)) + ALR(3)
    else
      LL(3) = MaxLL(LHH(1:2))
    endif
  endif
  if (LHH(3) - LLg(2) > 2*TA .and. LHH(3)<0) then  ! MAX(LLg(3), LLg(2))
    if (InclAge .and. ALR(2)<777) then
      LL(2) = LHH(3) + ALR(2)
    else
      LL(2) = LHH(3)
    endif
  endif
endif 

if (focal==4 .and. MaxLL(LL)==LL(4) .and. BY(A)>=0 .and. InclAge) then  ! don't rely on weak age prior only
  call BestRel(LL, focal, topX(1), dLL(1))
  call BestRel(LLg, focal, topX(2), dLL(2))
  if (topX(1) == focal .and. topX(2) /= focal) then
    if (LLg(4)-LLg(3)<.1 .and. LLg(3)<0 .and. (ALR(4) - ALR(3)) < 2*TA) then  !  .and. AgePhase<2
      LL(4) = 222 
    else if (LLg(4)-LLg(5)<.1 .and. LLg(5)<0 .and. (ALR(4) - ALR(5)) < 2*TA) then
      LL(4) = 222  ! relies on age prior only. 
    else if (LL(6) == LLz(1) .and. LLz(1)<0 .and. &
     (LLg(4) - (LLz(1)-ALRz(1)) < .1 .and. ALR(4) - ALRz(1) < 2*TA) .or. &
     LLg(4) - (LLz(1)-ALRz(1))<2*TA) then  ! GGP  
      LL(4) = 222
    endif
  endif
endif  ! TODO?: use age prior if last round. 

if (ANY(Parent(A,:)<=0)) then   ! one of Bi parent of A?
  ParTmp = Parent(A,:)
  do m=1,2
    if (Parent(A,m)<=0) then
      do i=1, ns(SB,k)
        Bi = SibID(i,SB,k)
        if (AgeDiff(A, Bi) < 0) cycle
        if (Sex(Bi)/=m .and. Sex(Bi)/=3)  cycle
        call CalcAgeLR(A,Sex(A), Bi,m, m,1, .TRUE., ALRq)
        if (ALRq < TF .or. ALRq==777)  cycle
        call getAncest(Bi,sex(Bi), AncBi)
        if (ANY(AncBi == A))  cycle
        call NewPar(A, m, Bi)                                         
        call CalcCLL(SB, k)                                                
        call CalcU(A, m, -SB, k, LLPO)
        call NewPar(A, m, ParTmp(m))
        call CalcCLL(SB, k)                   
        call CalcLind(A)
        if (LLPO < 0 .and. LLPO > MaxLL(LL)) then
          if (focal==1) then
            LL(1) = 222
          else
            LL(1) = LLPO
          endif
        endif
      enddo
    endif
  enddo
endif

end subroutine CheckAdd 

! #####################################################################

subroutine Qmerge(SA, SB, k, LR)
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LR
integer :: l, x, y
double precision :: PrL(nSnp), PrX(3), PrXY(3,3)

PrL = 0D0
do l=1,nSnp
  do x=1,3
    PrX(x) = XPr(1,x,l,SA,k)*XPr(1,x,l,SB,k)* AHWE(x,l)
    do y=1,3
      PrXY(x,y) = XPr(1,x,l,SA,k)*XPr(1,y,l,SB,k)* AHWE(x,l) * AHWE(y,l)
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrX)) - LOG10(SUM(PrXY))
enddo

LR = SUM(PrL)

end subroutine Qmerge

! #####################################################################

subroutine CheckMerge(SA, SB, kA, kB, InclAge, LL, focal, FSM) 
use Global
implicit none

integer, intent(IN) :: SA, SB, kA, kB, focal
logical, intent(IN) :: InclAge
double precision, intent(OUT) :: LL(7)
logical, intent(OUT) :: FSM
double precision :: LLg(7), LLtmp(2), ALR(7), LLx(6), LLz(2,2), LRHS, &
  LLM(5), LLMo(5), LLHA(3), dLH(nS(SB,kB)), ALRx(6), LLC, ALRtmp, &
  dx(maxSibSize), LLHHA(2), LLP
integer :: i, j, x, Par(2), AncA(2,mxA), AncB(2,mxA), ix, tmpGP
logical :: ShareOpp, ShareSib, con

LL = 999D0
LLg = 999D0
ALR = 999D0
ShareOpp = .FALSE.
ShareSib = .FALSE.                 
FSM = .FALSE.  ! merge both k & 3-k
if (kA /= kB) then
  LL(1) = 777
  if (focal==1)  return
endif
do i=1, nS(SA, kA)
  do x=1,2
    if (SibID(i, SA, kA)==GpID(x,SB,kB)) then
      LL(1) = 777
      exit
    endif
  enddo
  do j=1, nS(SB, kB)
    if (SibID(j, SB, kB)==GpID(1,SA,kA) .or. &
      SibID(j, SB, kB)==GpID(2,SA,kA)) then
      LL(1) = 777
      exit
    endif
    if (AgeDiff(SibID(i,SA,kA), SibID(j,SB,kB))==999) cycle
    if (AgePriorM(ABS(AgeDiff(SibID(i,SA,kA), &
      SibID(j,SB,kB)))+1,kA) == 0.0) then   
      LL(1) = 777
      exit
    endif
    if (LL(1)==777) exit
    if (kA/=kB) then
      if (SibID(i, SA, kA)==SibID(j,SB,kB)) then
        ShareSib = .TRUE.
      endif
    else if (kA==kB) then    
     if (Parent(SibID(i,SA,kA),3-kA)==Parent(SibID(j,SB,kB),3-kB) .and.&
        Parent(SibID(i, SA, kA), 3-kA) /=0) then
        ShareOpp = .TRUE.
      endif
    endif
  enddo
enddo 
if (LL(1) == 777 .and. focal==1) return

if (focal==1) then
  do x=1,2
    if (GpID(x,SB,kB)/=0) then
      call CalcAgeLR(-SA,kA, GpID(x,SB,kB),x, 0,1, .TRUE., ALRtmp)
      if (ALRtmp == 777) then
        LL(1) = 777
        exit
      endif
    endif
    if (GpID(x,SA,kA)/=0) then
      call CalcAgeLR(-SB,kB, GpID(x,SA,kA),x, 0,1, .TRUE., ALRtmp)
      if (ALRtmp == 777) then
        LL(1) = 777
        exit
      endif
    endif
  enddo
  if (LL(1) == 777) return
endif

AncA = 0
AncB = 0
if (focal==1) then
  call GetAncest(-SA, kA, AncA)
  call GetAncest(-SB, kB, AncB)           
  if (ANY(AncB(kA,:)==-SA) .or. ANY(AncA(kB,:)==-SB)) then
    LL(1) = 777
    return
  endif
endif       
  
if (.not. ShareOpp .and. .not. ShareSib) then
  call Qmerge(SA, SB, kB,  LRHS)
  if (LRHS < 2*TF*MAX(nS(SA,kA), nS(SB,kB))) then
    LL(1) = 777
  endif
   if (LL(1) == 777 .and. focal==1) return
endif
LLg(1) = LL(1)

 call CalcU(-SA,kA, -SB,kB, LL(7))
 
if (LL(1)/=777 .and. kA==kB) then
  call MergeSibs(SA, SB, kA, LLg(1))   ! SB parent of A's
  if (focal==1 .and. (LLg(1) > 0 .or. LLg(1) - LL(7) < TA)) then
    LL(1) = LLg(1)
    return
  endif      
  call CalcALRmerge(SA, SB, kA, ALR(1))
  if (ALR(1) == 777) then
    LL(1) = 777
  else if (LLg(1) < 0 .and. InclAge) then  
    LL(1) = LLg(1) + ALR(1)
  else
    LL(1) = LLg(1)
  endif
else
  LL(1) = 777
endif
if (focal==1 .and. (LL(1)==777 .or. (LL(1) - LL(7)) < TA)) return

  call CalcAgeLR(-SB,kB, -SA,kA, 0,1, .TRUE., ALR(2))
if (ALR(2) /= 777) then
   call addFS(0, SA, kA, SB, kB, LLg(2), ix, dx)  ! SB FS with an A
  if(complx>0)  call PairUA(-SB, -SA, kB, kA, LLg(3))  ! SB HS with an A                      
  do x=2,3
    if (LLg(x) < 0 .and. InclAge) then
      LL(x) = LLg(x) + ALR(2)
    else
      LL(x) = LLg(x)
    endif
  enddo
else
  LL(2:3) = 777
endif
 
LLtmp = 999D0
tmpGP = 0                   
call CalcAgeLR(-SA,kA, -SB,kB, 0,1, .TRUE., ALR(4))
if (ALR(4)/=777) then
  if (focal/=4 .or. complx==0) call addFS(0, SB, kB, SA, kA, LLtmp(1), ix, dx)   ! SB GP of A's
  if (focal==4) then  ! allow for replacement
    tmpGP = GpID(kB,SA,kA)
    GpID(kB,SA,kA) = 0
    call CalcCLL(SA,kA)
  endif
  if(complx>0)  call PairUA(-SA, -SB, kA, kB, LLtmp(2))  ! SB GP of A's
  if (focal==4) then  ! allow for replacement
    GpID(kB,SA,kA) = tmpGP
    call CalcCLL(SA,kA)
  endif   
  LLg(4) = MaxLL(LLtmp)                                                                       
  if (LLg(4) < 0 .and. InclAge) then
    LL(4) = LLg(4) + ALR(4)
  else
    LL(4) = LLg(4)
  endif
else
  LL(4) = 777
endif
 
call CalcAgeLR(-SA,kA, -SB,kB, 0,2, .TRUE., ALR(5))
if (ALR(5) /= 777) then 
  if(complx>0)  call ParentHFS(0, SA, kA, SB, kB,3, LLg(5))  ! SB FA of A's 
  if (InclAge .and. LLg(5)<0) then                                     
    LL(5) = LLg(5) + ALR(5)
  else
    LL(5) = LLg(5)
  endif
else
  LL(5) = 777
endif

LLx = 999D0
ALRx = 0
do x=1,4               
  if (x==1 .or. x==2) then
    if (complx==0) cycle
    call CalcAgeLR(-SA,kA, -SB,kB, x,3, .TRUE., ALRx(x))
    if (ALRx(x) /= 777) then
      call ParentHFS(0, SA, kA, SB, kB, x, LLx(x))
    endif
  else if (x==3) then
    call CalcAgeLR(-SA,kA, -SB,kB, 3-kB,4, .TRUE., ALRx(x))
    if (ALRx(x) /= 777) then
      call dummyGP(SA, SB, kA, kB, LLx(3))  ! SB GGP of A's
    endif      
  else if (x==4) then
    call CalcAgeLR(-SB,kB, -SA,kA, 3-kA,4, .TRUE., ALRx(x))
    if (ALRx(x) /= 777) then
      call dummyGP(SB, SA, kB, kA, LLx(4))  ! SA GGP of B's
    endif 
  endif
  
  if (ALRx(x) /= 777) then
    if (LLx(x) < 0 .and. InclAge) then
      LLX(x) = LLX(x) + ALRx(x)           
    endif
  else
    LLX(x) = 777
  endif
enddo

LLz = 999D0
do x=1,2
  if (GpID(x, SA, kA) > 0) then   ! TODO: more general
    do i=1,2
      if (GpID(i,SB,kB)/=0 .and. Parent(GpID(x, SA, kA), i)/=0 .and. &
       GpID(i,SB,kB) /= Parent(GpID(x, SA, kA),i)) then
        LLz(x,2) = 777
      endif
    enddo
    if (Parent(GpID(x, SA, kA), kB)==-SB) then
      LLz(x,2) = 777
    endif
    if (LLz(x,2)/=777) then    
      call CalcU(-SB, kB, GpID(x,SA,kA), x, LLz(x,1))
      call PairUA(-SB, GpID(x,SA,kA), kB, 3, LLz(x,2))
      call CalcAgeLR(-SB,kB, GpID(x,SA,kA),x,0,2, .TRUE., ALRx(4+x))
    endif
    if (LLz(x,2) < 0 .and. ALRx(4+x)/=777) then
      LLx(4+x) = LL(7) + LLz(x,2) - LLz(x,1)
      if (InclAge) then
        LLx(4+x) = LLx(4+x) + ALRx(4+x)
      endif
    endif
  endif
enddo
LL(6) = MaxLL(LLx)  ! most likely 3rd degree relative

dLH = 999D0                      
if (complx>0 .and. LL(4)<0 .and. focal/=4 .and. LLtmp(1)<0 .and. &
  LLtmp(1)>=LLtmp(2) .and. LLtmp(1) > MaxLL((/LL(1:3), LL(5:7)/))) then 
  do j=1, nS(SB, kB)
    if (GpID(3-kB,SA,kB)==0) then
      shareOpp = .TRUE.
      par(1) = Parent(SibID(j, SB, kB), 3-kB)
    else if (Parent(SibID(j, SB, kB), 3-kB)==GpID(3-kB,SA,kB) .or. &
      Parent(SibID(j, SB, kB), 3-kB)==0) then
      shareOpp = .TRUE.
      par(1) = GpID(3-kB,SA,kB)
    else
      shareOpp = .FALSE.
    endif
    if (shareOpp) then
      call PairUA(-SA, SibID(j, SB, kB), kA, 3, LLHA(1))
      call PairUA(-SA, SibID(j, SB, kB), kA, 3-kB, LLHA(2))
      call CalcU(-SA, kA, SibID(j, SB, kB), kB, LLHA(3))
      if (LLHA(1)<0)  dLH(j) = LLHA(1) - MaxLL(LLHA(2:3))
    endif
  enddo
  if (MAXVAL(dLH, MASK=dLH<999) < TA) then
    LL(4) = LLtmp(2)
    if (InclAge) LL(4) = LL(4) + ALR(4)
  endif
endif

LLM = 999D0
LLMo = 999D0
Par = 0
if (kA == kB .and. focal==1 .and. (ABS(MaxLL(LL) - LL(1))<TA .or. &
  Complx==0)) then     
  call FSMerge(SA,SB,kA, LLM)
  if (Complx/=2)  LLM(5) = 555   ! merge via 3-k + par HS                                                       
  LLM(1) = MaxLL((/LLM(1), LL(7)/))  ! do not merge  
  LLM(2) = MaxLL((/LLM(2), LL(1)-ALR(1), LLM(5)/))  ! merge via k
  call getFSpar(SA, kA, .FALSE., Par(1))
  call getFSpar(SB, kB, .FALSE., Par(2))
  if (par(1)<0 .and. par(2)<0) then
    if(.not. (ns(SA,kA)==ns(-par(1),3-kA) .and. ns(SB,kB)==ns(-par(2),3-kB))) then                                   
      call FSMerge(-par(1),-par(2),3-kA, LLMo)
      if (Complx/=2)  LLMo(5) = 555                             
      call CalcU(Par(1), 3-kA, Par(2), 3-kA, LLMo(1))  ! more accurate
      call MergeSibs(-par(1), -par(2), 3-kA, LLMo(2)) 
    endif
  endif
  if (Complx==0) then   ! FS merge only
    LL(1) = LLM(4)
    if (InclAge) LL(1) = LL(1) + ALR(1) 
    FSM = .TRUE.
  endif
  if (Complx==2 .and. LL(1) /= 222 .and.  par(1)<0 .and. par(2)<0) then
    if (nS(SA,kA) + nS(SB,kB) == nS(-Par(1),3-kA)+ns(-Par(2),3-kA)) then     
      if (LLM(4)<0) then
        call FSHC(-SA,-SB,kA,LLC)
        if (LLC > LLM(4) .and. LLC<0)  LLM(4) = LLC
      endif
      call clustHSHA(SA, SB, kA, LLHHA(1))
      if (LLHHA(1) - LLM(4) > TA)  LL(1) = 222
      call clustHSHA(SB, SA, kA, LLHHA(2))
      if (LLHHA(2) - LLM(4) > TA)  LL(1) = 222
    endif
  endif
  if (MaxLL(LLM(1:4))==LLM(4) .and. LLM(4)-LLM(2) >TA*MIN(nS(SA,kA),nS(SB,kB)) .and. &
    (ALL(LLMo==999) .or. LLMo(4)-LLMo(3) >TA*MIN(nS(SA,kA),nS(SB,kB)))) then
!    (LLM(4)-LLM(2)) - (LLMo(3) - LLMo(1)) > 2*TA) then  
    LL(1) = LLM(4)  ! FS merge most likely - go ahead.
    if (InclAge) LL(1) = LL(1) + ALR(1) 
    FSM = .TRUE.
  else if (Complx>0 .and. LLM(3)<0 .and. LLM(3)-LLM(1) > TA) then
    if (LLM(3)-LLM(4) > TA) then
      LL(1) = 222 ! likely that opp. parent need to be merged, 
    else if (Par(1) < 0 .and. Par(2)<0) then
      if (nS(SA,kA)==nS(-Par(1),3-kA) .and. nS(SB,kB)==ns(-Par(2),3-kA)) then   ! 2 FS groups
        LL(1) = 222
      else if (nS(-Par(1),3-kA)+ns(-Par(2),3-kA) > nS(SA,kA)+nS(SB,kB)) then
        if (LLMo(2)>0 .or. LLMo(1) - LLMo(2) > TA*MIN(ns(-Par(1),3-kA),ns(-Par(2),3-kA))) then
          LL(1) = LL(1)  ! opp. merge unlikely
        else                                                                       
          LL(1) = 222   ! LLM(3) and (4) not comparable
        endif
      endif
    endif  
  else if (Complx>0 .and. LLMo(2)<0) then
    if ((LLM(2)-LLM(1)) - (LLMo(2)-LLMo(1)) >TA*MIN(nS(SA,kA),nS(SB,kB))) then
      LL(1) = LL(1)
    else
      LL(1) = 222 
    endif
  endif
endif

if (focal==1 .and. (ABS(MaxLL(LL) - LL(1))<TA)) then  ! one of Bi is SA, or vv?
  do i=1, ns(SA,kA)
    if (AgeDiff(SibID(i,SA,kA), SibID(1,SB,kB)) < 0) cycle
    call AddParent(SibID(i,SA,kA), SB, kB, LLP)
    if (LLP < 0 .and. (LLP + CLL(SA,kA) - Lind(SibID(i,SA,kA))) > LL(1)) then
      LL(1) = 222
      exit
    endif
  enddo
  if (LL(1)<0) then
    do j=1, ns(SB,kB)
      if (AgeDiff(SibID(j,SB,kB), SibID(1,SA,kA)) < 0) cycle
      call AddParent(SibID(j,SB,kB), SA, kA, LLP)
      if (LLP < 0 .and. (LLP + CLL(SB,kB) - Lind(SibID(j,SB,kB))) > LL(1)) then
        LL(1) = 222
      endif
    enddo
  endif 
endif

end subroutine CheckMerge 

! #####################################################################
subroutine getFSpar(SA, kA, strict, par)  
! all individuals in SA are FS to eachother
use Global
implicit none

integer, intent(IN) :: SA,  kA
logical, intent(IN) :: strict                            
integer, intent(OUT) :: Par
integer :: i, j, ParV(ns(SA,kA))

Par = 0
ParV = 0
do i=1, nS(SA,kA)
  if (Parent(SibID(i,SA,kA), 3-kA)/=0) then
    Par = Parent(SibID(i,SA,kA), 3-kA)
    if (strict) then
      do j= i+1, nS(SA, kA)
        if (Parent(SibID(j,SA,kA), 3-kA) /= Par .and. &
         Parent(SibID(j,SA,kA), 3-kA)/=0) then
          Par = 0
          return
        endif
      enddo
    else 
      ParV(i) = Par
    endif
  endif
enddo

if (.not. strict) then ! > half by same opp. parent?
  Par = 0
  do i=1, nS(SA,kA)
    if (COUNT(ParV == ParV(i)) > nS(SA,kA)/2.d0) then
      Par = ParV(i)
      return
    else if (COUNT(ParV == ParV(i)) == nS(SA,kA)/2.d0 .and. ParV(i)<0) then
      Par = ParV(i)
      return
    endif
  enddo
endif

end subroutine getFSpar

! #####################################################################

subroutine OppMerge(SA, k, LL)  ! could opposing parents of SA all be the same dummy parent?
use Global
implicit none

integer, intent(IN) :: SA, k
double precision, intent(OUT) :: LL ! of SA
integer :: i, l, x, y, m,u, opPar(ns(SA,k)), GPY(2)
double precision :: PrL(nSnp), PrSA(3), PrXY(3, 3), PrGY(3,2)

opPar = 0
GPY = 0
do i=1, ns(SA,k)
  opPar(i) = Parent(SibID(i, SA, k), 3-k)
enddo
if (ANY(opPar > 0)) then
  LL = 444
  return
else
  do i=1, ns(SA,k)
    if (opPar(i) < 0) then
      if (ns(-opPar(i), 3-k) > COUNT(opPar == opPar(i))) then
        LL = 999  ! TODO?  
        return
      else
        do m=1,2
          if (GpID(m,-OpPar(i),3-k) /= GPY(m)) then
            if (GPY(m) == 0) then
              GPY(m) = GpID(m,-OpPar(i),3-k)
            else
              LL = 777
              return
            endif
          endif
        enddo
      endif
    endif
  enddo
endif

LL = 999D0
do l=1, nSnp
  call ParProb(l, -SA, k, -1, 0, PrSA)
  do m=1,2
    call ParProb(l, GPY(m), 3-k, 0, 0, PrGY(:,m))  
  enddo
  do x=1,3
    do y=1,3
      do u=1,3
        PrXY(x,y) = PrSA(x) * SUM(AKA2P(y, u, :) * PrGY(u,1) * PrGY(:,2))
      enddo
      do i=1, ns(SA,k)
        if (Genos(l, SibID(1,SA,k))/=-9) then
          PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l, SibID(i,SA,k)), x, y)
        endif
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine OppMerge

! #####################################################################                                                                      
subroutine FSmerge(SA,SB,k, LL)  
! calc LL if SA and SB merged via both pat & mat
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LL(5) ! 1:not, 2: via k, 3: via 3-k, 4:both, !!5: 3-k + par HS
integer :: l, x, y, i, u,v, G(2,2),z, m, Par(2), SX(2), MaybeOpp(2)
double precision :: PrL(nSnp,5), PrXY(3,3), PrUV(3,3), PrXV(3,3,3,3,5),&
  PrG(3,2,2), PrX(3,2), PrTmp(3), ALR, PrY(3,2)
logical :: DoParHS                                                                 

! TODO: currently assumes no gps of sibship 3-k, no close inbreeding
LL = 999D0
! check if all FS
SX = (/SA, SB/)
MaybeOpp = 0
do i=1,2
  call getFSpar(SX(i), k, .TRUE., Par(i))
  if (Par(i)>0) cycle
  if (Par(i)==0 .and. ANY(Parent(SibID(1:nS(SX(i),k),SX(i),k),3-k)>0))&
    cycle
  MaybeOpp(i) = 1
enddo
if (MaybeOpp(1)==1 .and. MaybeOpp(2)==1) then   
  if (Par(1)==Par(2) .and. Par(1)/=0) then 
    if (ALL(Parent(SibID(1:nS(SA,k),SA,k),3-k)==Par(1)) .and. &
      ALL(Parent(SibID(1:nS(SB,k),SB,k),3-k)==Par(2))) then
      MaybeOpp = 0
    endif
  else if (Par(1)<0 .and. Par(2)<0) then
    call CalcALRmerge(-Par(1), -Par(2), 3-k, ALR)
    if (ALR==777) MaybeOpp = 0
    if (nS(-Par(1),3-k) > ns(SA,k) .or. nS(-Par(2),3-k) > ns(SB,k)) then
      LL = 444
    endif
  endif
endif
if (ANY(MaybeOpp==0) .or. ALL(LL==444)) return

G = 0
do i=1,2
  if (GpID(i,SA,k)/=0) then
    if(GpID(i,SA,k)/=GpID(i,SB,k) .and. GpID(i,SB,k)/=0) then
      G(i,k) = 0  ! shouldn't happen
    else
      G(i,k) = GpID(i,SA,k)
    endif
  else
    G(i,k) = GpID(i,SB,k)
  endif
  if (Par(1)<0) then
    G(i,3-k) = GpID(i, -Par(1),3-k)
    if (GpID(i, -Par(1),3-k)/=0) then
      if (Par(2) < 0) then
        if (GpID(i, -Par(2),3-k) /= G(i,3-k) .and. &
          GpID(i, -Par(2),3-k)/=0) then
          G(i,3-k) = 0   ! shouldn't happen
        else if (G(i,3-k)==0 .and. GpID(i, -Par(2),3-k)/=0) then
          G(i,3-k) = GpID(i, -Par(2),3-k)
        endif
      endif
    endif
  endif
enddo

 if (ALL(GPID(:,SA,k)==0) .and. ALL(GPID(:,SB,k)==0)) then
  DoParHS = .TRUE.
else
  DoParHS = .FALSE.
endif    

PrL = 0D0
do l=1,nSnp 
  do m=1,2
    do i=1,2
      call ParProb(l, G(i,m), i, 0, 0, PrG(:,i,m))
    enddo
    do x=1,3
      do z=1,3
        PrTmp(z) = SUM(AKA2P(x,:,z) * PrG(:,1,m) * PrG(z,2,m))
      enddo
      PrX(x,m) = SUM(PrTmp)
    enddo
  enddo
  do i=1,2
    call ParProb(l, Par(i), 3-k, -1, 0, PrY(:,i))
  enddo
  do x=1,3  ! P1
    do y=1,3  ! P2
      PrXY(x,y) = 1  ! XPr(2,x,l, sA,k) * AHWE(y,l)
      do i=1,nS(SA,k)
        if (Genos(l, SibID(i,SA,k))/=-9) then
          PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l, SibID(i,SA,k)), x, y)
        endif
      enddo
    enddo
  enddo
  do u=1,3
    do v=1,3
      PrUV(u,v) = 1  ! XPr(2,u,l, sB,k) * AHWE(v,l)
      do i=1,nS(SB,k)
        if (Genos(l, SibID(i,SB,k))/=-9) then
          PrUV(u,v) = PrUV(u,v) * OKA2P(Genos(l, SibID(i,SB,k)), u, v)
        endif
      enddo
    enddo
  enddo

  PrXV = 0
  do x=1,3 
    do y=1,3
      do u=1,3
        do v=1,3
          PrXV(x,y,u,v,1) = PrXY(x,y) * XPr(2,x,l, sA,k) * PrY(y,1) * &
            PrUV(u,v) * XPr(2,u,l, sB,k) * PrY(v,2)
          PrXV(x,y,x,v,2) = PrXY(x,y) * PrX(x,k) * PrY(y,1) * &
            PrUV(x,v) * PrY(v,2)
        enddo
        PrXV(x,y,u,y,3) = PrXY(x,y) * XPr(2,x,l, sA,k) * PrX(y,3-k) * &
          PrUV(u,y) * XPr(2,u,l, sB,k)
        if (DoParHS) then
          do z=1,3
!           PrTmp(z) = SUM(AKA2P(x,z,:) *AKA2P(u,z,:) *AHWE(z,l)*AHWE(:,l))
            PrTmp(z) = AKAP(x,z,l) * AKAP(u,z,l) * AHWE(z,l)
          enddo
          PrXV(x,y,u,y,5) = PrXY(x,y) * PrX(y,3-k) * PrUV(u,y) * &
            SUM(PrTMP)
        endif
      enddo           
      PrXV(x,y,x,y,4) = PrXY(x,y) * PrX(x,k) * PrX(y,3-k) * PrUV(x,y)
    enddo
  enddo
  do x=1,5
    PrL(l,x) = LOG10(SUM(PrXV(:,:,:,:,x)))
  enddo
enddo
LL = SUM(PrL,DIM=1)
if (.not. DoParHS)  LL(5) = 777D0                                 

end subroutine FSmerge

! #####################################################################

subroutine MakeFS(A, B)
use Global
implicit none

integer, intent(IN) :: A,B
integer :: x, i, j, Ai, Bj

Ai = 0
Bj = 0
if (nFS(A)>0) then
  Ai = A
else
  Ai = FSID(maxSibSize+1, A)
endif
if (nFS(B)>0) then
  Bj = B
else
  Bj = FSID(maxSibSize+1, B)
endif

if (ANY(FSID(1:nFS(Ai),Ai)==B) .or. ANY(FSID(1:nFS(Bj),Bj)==A)) then
  return ! already are FS.
endif

i = MIN(Ai,Bj)
j = MAX(Ai,Bj)
do x=1, nFS(j)   
  FSID(nFS(i)+x, i) = FSID(x, j)
  FSID(maxSibSize+1, FSID(x,j)) = i
enddo
nFS(i) = nFS(i) + nFS(j)
FSID(maxSibSize+1,i) = i    ! 'primary' sib
FSID(:,j) = 0
FSID(1,j) = j
FSID(maxSibSize+1,j) = i
nFS(j) = 0

end subroutine MakeFS

! #####################################################################
subroutine DoAdd(A, SB, k)
use Global
implicit none

integer, intent(IN) :: A, SB, k
integer :: i, n, j

if (nS(SB,k) +1 >= maxSibSize) then
  call Erstop("reached MaxSibshipSize, please increase") 
endif

Parent(A, k) = -SB
if (.not. ANY(SibID(1:nS(SB,k),SB,k)==A)) then
  SibID(nS(SB,k)+1, SB, k) = A  ! add A to sibship
  nS(SB,k) = nS(SB,k) + 1
endif
do n=1, nS(SB,k)  ! check for FS
  i = SibID(n,SB,k)
  if (i==A) cycle 
  if (nFS(i)==0) cycle
  if (Parent(A, 3-k)/=0 .and. Parent(A, 3-k)==Parent(i, 3-k)) then
    call MakeFS(A, i)
  endif
enddo
 call calcCLL(SB,k)
 call CalcLind(A)

do n=1,nS(SB,k)   ! update LL of connected sibships
  i = SibID(n,SB,k)
  if (Parent(i,3-k) < 0) then
    call CalcCLL(-Parent(i,3-k), 3-k)
  endif                    
  call CalcLind(i)
enddo
 call calcCLL(SB,k)
 call CalcLind(A)

end subroutine DoAdd

! #####################################################################
recursive subroutine DoMerge(SA, SB, k, valid)
use Global
implicit none

integer, intent(IN) :: SA, SB, k
logical, intent(OUT) :: valid                             
integer :: i, j, n, m, x, AncA(2,mxA), AncB(2,mxA)
logical :: OK                                                    

if (SA == SB) return

valid = .TRUE.
if (SA/=0) then
  call GetAncest(-SA, k, AncA)
  call GetAncest(-SB, k, AncB)
  do i=1, nS(SA, k)
    if (ANY(AncB == SibID(i,SA, k))) then
      valid = .FALSE.
      exit
    endif
  enddo
  do j=1, nS(SB, k)
    if (ANY(AncA == SibID(j,SB, k))) then
      valid = .FALSE.
      exit
    endif
  enddo
endif
if (valid .eqv. .FALSE.) then
  call Erstop("Pedigree loop created")
endif

if (SA/=0) then
  if (nS(SA,k) + nS(SB,k) >= maxSibSize) then
    call Erstop("reached maxSibSize")
  endif    
  do n=1,nS(SA,k)    ! check for FS
    i = SibID(n,SA,k)
    if (nFS(i)==0) cycle
    do m=1, nS(SB,k)
      j = SibID(m,SB,k)
      if (nFS(j)==0) cycle
      if (Parent(i, 3-k)/=0 .and. Parent(i, 3-k)==Parent(j, 3-k)) then
        call MakeFS(i, j)
      endif
    enddo
  enddo
  do m=1,nS(SB,k)  ! add sibship SB to SA
    SibID(nS(SA,k)+m, SA, k) = SibID(m, SB, k)
    Parent(SibID(m, SB, k), k) = -SA
    do i=1,2
      if (GpID(i, SA, k)==0 .and. GpID(i, SB, k)/=0) then
        GpID(i, SA, k) = GpID(i, SB, k)  !checked for mismatches earlier
      endif  ! else keep GpID(i,SA,k)
    enddo
  enddo
  nS(SA,k) = nS(SA,k) + nS(SB,k)
  
  call calcCLL(SA,k)
  do n=1,nS(SA,k) 
    i = SibID(n,SA,k)
    if (Parent(i,3-k) < 0) then
      call CalcCLL(-Parent(i,3-k), 3-k)
    endif                    
    call CalcLind(i)
  enddo
  call calcCLL(SA,k)
  do n=1,nS(SA,k)                 
    call CalcLind(SibID(n,SA,k))
  enddo
endif
 
do x=SB, nC(k)-1  !remove cluster SB, shift all subsequent ones
  SibID(:, x, k) = SibID(:, x+1, k)
  nS(x, k) = nS(x+1, k)
  GpID(:, x,k) = GpID(:, x+1,k)
  do n=1, nS(x,k)
    Parent(SibID(n,x,k),k) = -x ! shift towards zero.
  enddo
  CLL(x,k) = CLL(x+1, k)
  XPr(:,:,:,x,k) = XPr(:,:,:,x+1,k)
  DumP(:,:,x,k) = DumP(:,:,x+1,k)
  DumBY(:,x,k) = DumBY(:,x+1,k)
enddo
SibID(:,nC(k),k) = 0
GpID(:,nC(k),k) = 0
nS(nC(k), k) = 0

do m=1,2  !fix GPs
  do n=1, nC(m)
    if (GpID(k, n, m) == -SB)  GpID(k, n, m) = -SA
    if (all(GpID(:,n,m)==0) .and. ns(n,m)==1) then
      call DoMerge(0, n, m, OK)
    endif                        
    do x=SB+1, nC(k)  
      if (GpID(k, n, m) == -x)  GpID(k, n, m) = -x+1 
    enddo
  enddo
enddo
nC(k) = nC(k) -1

end subroutine DoMerge

! #####################################################################

subroutine DoFSMerge(SA, SB, k)   ! merge via k .and. k-3
use Global
implicit none

integer, intent(IN) :: SA, SB, k
integer :: i, j, ParA, ParB
logical :: OK            

! assume all checks done beforehand
 call getFSpar(SA, k, .TRUE., ParA)
 call getFSpar(SB, k, .TRUE., ParB)
OK = .TRUE.
if (ParA < 0 .and. ParB < 0) then
    call DoMerge(-ParA, -ParB, 3-k, OK)
    if (.not. OK)  return                       
else if (ParA==0 .and. ParB==0) then
    nC(3-k) = nC(3-k)+1  ! new sibship       
    nS(nC(3-k), 3-k) = nS(SA,k) + nS(SB,k)
    j = 0
    do i=1, nS(SA, k)
        j = j+1
        SibID(j, nC(3-k), 3-k) = SibID(i, SA, k)                                       
        call NewPar(SibID(i, SA, k), 3-k, -nC(3-k))
    enddo
    do i=1, nS(SB, k)
        j = j+1
        SIbID(j, nC(3-k), 3-k) = SibID(i, SB, k)
        call NewPar(SibID(i, SB, k), 3-k, -nC(3-k))
    enddo 
    call CalcCLL(nC(3-k), 3-k)
    call CalcCLL(SA, k)
    call CalcCLL(SB, k)
    call CalcCLL(nC(3-k), 3-k)
else if (ParA < 0 .and. ParB == 0) then
    do i=1, nS(SB, k)
        call DoAdd(SibID(i, SB, k), -ParA, 3-k)
    enddo
else if (ParB < 0 .and. ParA == 0) then
    do i=1, nS(SA, k)
        call DoAdd(SibID(i, SA, k), -ParB, 3-k)
    enddo
! else not implemented yet
endif

 call DoMerge(SA, SB, k, OK)  ! takes care of MakeFS

end subroutine DoFSMerge

! #####################################################################

subroutine getOff(P, kP, dums, nOff, Off, sxOff)  ! list all offspring for parent P
use Global
implicit none

integer, intent(IN) :: P, kP
logical, intent(IN) :: dums  ! include dummy offspring
integer, intent(OUT) :: nOff, sxOff(maxSibSize), Off(maxSibSize)
integer :: i, k, m, s

nOff = 0
Off = 0
if (P==0) return                
do k=1,2
  if (P>0 .and. kP/=1 .and. kP/=2) then
    if (Sex(P)/=3 .and. Sex(P)/=k) cycle
  else if (k/=kP) then 
    cycle
  endif
  do i=1, nInd
    if (Parent(i,k) == P) then
      nOff = nOff + 1
      Off(nOff) = i
      sxOff(nOff) = Sex(i)
    endif
    if (nOff == maxSibSize) then
      call rexit("reached MaxSibshipSize, please increase")
    endif                            
  enddo
  if (dums) then
    do m=1,2
      do s=1,nC(m)
        if (GpID(k,s,m) == P) then
          nOff = nOff + 1
          Off(nOff) = -s
          sxOff(nOff) = m 
        endif
        if (nOff == maxSibSize) then
          call rexit("reached MaxSibshipSize, please increase")
        endif 
      enddo
    enddo
  endif
enddo    

end subroutine getOff

! #####################################################################
subroutine CalcU(A, kAIN, B, kBIN, LL)  ! A, SB, k, SA, LL
use Global
implicit none

integer, intent(IN) :: A, kAIN, B, kBIN
double precision, intent(OUT) :: LL
integer :: m, n, cat, par(2), Ai, Bj, SA, SB, kA, kB, i, tmpGP
logical :: swap, con, OpG, conP

LL = 999D0
con = .FALSE.
if (A>0) then
  call CalcLind(A)
else if (A<0) then
  kA = kAIN
  call CalcCLL(-A, kA)
endif
if (B>0) then
  call CalcLind(B)
else if (B<0) then
  kB = kBIN
  call CalcCLL(-B, kB)
endif
!==================================

if (A==0) then
  if (B==0) then
    LL = 0
  else if (B>0) then
    LL = Lind(B)
  else if (B<0) then
    LL = CLL(-B, kB)
  endif
  return
else if (B==0) then
  if (A>0) then
    LL = Lind(A)
  else if (A<0) then
    LL = CLL(-A,kA)
  endif
  return
else if (A>0 .and. B<0) then
  if (Parent(A,kB)==B) then
    LL = CLL(-B,kB)
    return
  else if (ANY(GpID(:,-B,kB) == A)) then
    LL = CLL(-B,kB) + Lind(A)  ! CLL already conditional on A
    return
  else if (ALL(Parent(A,:)>=0)) then  
    LL = Lind(A) + CLL(-B, kB)
    return
  else
    call Connected(A,1,B,kB, con)
    if (.not. con) then
      LL = Lind(A) + CLL(-B, kB)
      return
    endif
  endif
else if (B>0 .and. A<0) then
  if (Parent(B,kA)==A) then
    LL = CLL(-A, kA)
    return
  else if (ANY(GpID(:,-A,kA) == B)) then
    LL = CLL(-A,kA) + Lind(B)  
    return
  else if (ALL(Parent(B,:)>=0)) then 
    LL = CLL(-A, kA) + Lind(B) 
    return
  else
    call Connected(B,1,A,kA, con)
    if (.not. con) then
      LL = CLL(-A, kA) + Lind(B) 
      return
    endif
  endif
endif

!==================================
! determine relationship between focal individuals/clusters

Ai = 0
Bj = 0
SA = 0
SB = 0
cat = 0

if (A>0 .and. B>0) then  ! == pairs ==
  do m=1,2
    if (Parent(A,m)/=0 .and. Parent(A,m)==Parent(B,m)) then
       par(m) = Parent(A,m)
    else
      par(m) = 0  ! unknown or unequal
    endif
  enddo

  if (par(1)/=0 .and. par(2)/=0) then
    cat = 2  ! FS
  else if (par(1)/=0 .or. par(2)/=0) then
    cat = 3  ! HS
  else 
    do m=1,2
      if (parent(A,m) < 0) then
        do n=1,2
          if (GpID(n, -parent(A,m), m) == B) then
            cat = 0 !4  ! already conditioned on.
          else
            if (GpID(n, -parent(A,m), m)==Parent(B, n) .and. &
             Parent(B, n)<0) then
              cat = 5
            endif
          endif
        enddo
      else if (parent(B,m) < 0) then
        if (ANY(GpID(:, -parent(B,m), m) == A)) then
          cat = 0 !4
          swap = .TRUE.
        else
          do n=1,2
            if (GpID(n, -parent(B,m), m) == Parent(A, n) .and. &
             Parent(A, n)<0) then
              cat = 5
              swap = .TRUE.
            endif
          enddo
        endif
      endif
    enddo
  endif

  if (cat==0 .or. cat==5) then  ! TODO? cat=5
    LL = Lind(A) + Lind(B)
    return
  else if (cat==2 .and. par(1)<0 .and. par(2)<0) then
    Ai = A
    Bj = B
    SA = -par(1)
    kA = 1
    SB = -par(2)
    kB = 2
    cat = 0
  else
    call Upair(A, B, cat, LL)
    return
  endif

else if (A>0 .and. B<0) then
  SB = -B
  Ai = A
  if (ALL(Parent(A,:) < 0)) then
    SA = -Parent(A,3-kB)
    kA = 3-kB
    do m=1,2
      do i=1,ns(-Parent(A,m),m)
        if (Parent(SibID(i,-Parent(A,m),m), 3-m)/=Parent(A,3-m)) then
          call Connected(SibID(i,-Parent(A,m),m),m,B,kB, conP)
          if (conP) then
            SA = -parent(A,m)
            kA = m
            exit
          endif
        endif
      enddo
    enddo
  else if (Parent(A,3-kB) < 0) then
    SA = -Parent(A,3-kB)
    kA = 3-kB    
  else if (Parent(A,kB) < 0) then
    SA = -Parent(A,kB)
    kA = kB
  endif ! else: Lind + CLL (earlier)
else if (B>0 .and. A<0) then
  SA = -A
  Bj = B
  if (ALL(Parent(B,:) < 0)) then
    SB = -Parent(B,3-kA)
    kB = 3-kA
    do m=1,2
      do i=1,ns(-Parent(B,m),m)
        if (Parent(SibID(i,-Parent(B,m),m), 3-m)/=Parent(B,3-m)) then
          call Connected(SibID(i,-Parent(B,m),m),m,A,kA, conP)
          if (conP) then
            SB = -parent(B,m)
            kB = m
            exit
          endif
        endif
      enddo
    enddo
  else if (Parent(B,3-kA) < 0) then
    SB = -Parent(B,3-kA)
    kB = 3-kA
  else if (Parent(B,kA) < 0) then
    SB = -Parent(B,kA)
    kB = kA 
  endif
else if (A<0 .and. B<0) then
  SA = -A  
  SB = -B
endif

cat = 0
swap = .FALSE.
if (GpID(kB, SA, kA) == -SB) then
  cat = 1  ! PO
else if (GpID(kA, SB, kB) == -SA) then
  cat = 1
  swap = .TRUE.
else 
  do m=1,2
    if (GpID(m, SA, kA)==GpID(m, SB, kB) .and. GpID(m, SA, kA)/=0) then
     if (GpID(3-m,SA,kA)==GpID(3-m,SB,kB) .and. GpID(3-m,SA,kA)/=0) then  
        cat = 2  ! FS
      else
        cat = 3  ! HS
      endif
    else 
      if (GpID(m, SA, kA)<0) then
        if (GpID(kB, -GpID(m, SA, kA), m) == -SB) then
          cat = 4  ! GP
        endif
      endif
      if (GpID(m, SB, kB)<0) then
        if (GpID(kA, -GpID(m, SB, kB),m) == -SA) then
          cat = 4
          swap = .TRUE.
        endif
      endif
    endif
  enddo  ! FA between SA, SB not currently considered.
endif

OpG = .FALSE.
if (con .and. cat==0) then
  if (A<0 .and. B>0) then
    do i=1, ns(-A,kA)
      if (Parent(SibID(i,-A,kA), 3-kA) < 0) then
        if (ANY(GpID(:,-Parent(SibID(i,-A,kA), 3-kA), 3-kA)==B)) then
          SB = -Parent(SibID(i,-A,kA), 3-kA)
          kB = 3-kA
          Bj = 0
          OpG = .TRUE.
        endif
      endif
    enddo
else if (A>0 .and. B<0) then
  do i=1, ns(-B,kB)
    if (Parent(SibID(i,-B,kB), 3-kB) < 0) then
        if (ANY(GpID(:,-Parent(SibID(i,-B,kB), 3-kB), 3-kB)==A)) then
          SA = -Parent(SibID(i,-B,kB), 3-kB)
          kA = 3-kB
          Ai = 0
          OpG = .TRUE.
          swap = .TRUE.
        endif
      endif
    enddo
  endif
endif

if (cat==0) then ! swap if BY(A) < BY(B)
  if (A>0 .and. B<0) then
    if (BY(A)>=0 .and. BY(A) < MINVAL(BY(SibID(1:nS(SB,kB),SB,kB)),&
    MASK=BY(SibID(1:nS(SB,kB),SB,kB))>=0)) then
      swap = .TRUE.
    endif
  else if (A<0 .and. B>0) then
    if (BY(B)>=0 .and. BY(B) < MINVAL(BY(SibID(1:nS(SA,kA),SA,kA)),&
    MASK=BY(SibID(1:nS(SA,kA),SA,kA))>=0)) then
      swap = .TRUE.
    endif
  else if (A<0 .and. B<0) then
    if (MAXVAL(BY(SibID(1:nS(SB,kB),SB,kB))) < MINVAL(BY(SibID(1:nS(SA,kA),SA,kA)),&
    MASK=BY(SibID(1:nS(SA,kA),SA,kA))>=0)) then
      swap = .TRUE.
    endif
  endif
endif

if (con .and. A<0 .and. B<0 .and. kA==kB) then
  do i=1,ns(-B,kB)
    if (Parent(SibID(i,-B,kB), 3-kB) < 0) then
      if (GpID(3-kA, -Parent(SibID(i,-B,kB), 3-kB), 3-kB) < 0) then
        tmpGP = GpID(3-kA, -Parent(SibID(i,-B,kB), 3-kB), 3-kB)
        if (ANY(Parent(SibID(1:ns(-A,kA),-A,kA), 3-kA) == tmpGP) .and. &
         .not. ANY(Parent(SibID(1:ns(-B,kB),-B,kB), 3-kB) == tmpGP)) then
          swap = .TRUE.
        endif
      endif
    endif
  enddo
endif                       
if (.not. swap) then
  call UClust(-SA, -SB, kA, kB, cat, Ai, Bj, LL)
else
  call UClust(-SB, -SA, kB, kA, cat, Bj, Ai, LL)
endif

if (opG) then
  if (B>0) then
    LL = LL - CLL(SB,kB) + Lind(B)
  else if (A>0) then
    LL = LL - CLL(SA,kA) + Lind(A)
  endif
endif

end subroutine CalcU

! #####################################################################

subroutine Upair(A, B, cat, LL)
use Global
implicit none

integer, intent(IN) :: A, B, cat
double precision, intent(OUT) :: LL
integer :: m, l, n, x, y, par(2)
double precision :: PrL(nSnp), PrP(3,2), PrPA(3), PrPB(3), PrXY(3,3)

LL = 999D0
 call CalcLind(A)
 call CalcLind(B)

do m=1,2
  if (Parent(A,m)/=0 .and. Parent(A,m)==Parent(B,m)) then
     par(m) = Parent(A,m)
  else
    par(m) = 0  ! unknown or unequal
  endif
enddo
      
PrL = 0D0
do l=1, nSnp
  if (Genos(l,A)==-9 .and. Genos(l,B)==-9) then
    cycle
  else if (Genos(l,A)==-9) then
    PrL(l) = LOG10(SUM(LindX(:,l,B)))
    cycle
  else if (Genos(l,B)==-9) then
    PrL(l) = LOG10(SUM(LindX(:,l,A)))
    cycle
  endif
  
  if (cat==2) then
    do m=1,2
      call ParProb(l, Par(m), m, A, B, PrP(:,m))
    enddo
    do x=1,3
      do y=1,3
        PrXY(x,y) = OKA2P(Genos(l,A),x,y) * OKA2P(Genos(l,B),x,y) * &
          PrP(x,1) * PrP(y,2)
      enddo
    enddo
  else if (cat==3) then  ! HS
    do m=1,2
      if (Par(m)==0) cycle
      call ParProb(l, Par(m), m, A, B, PrP(:,m))
      call ParProb(l, Parent(A, 3-m), 3-m, A, 0, PrPA)
      call ParProb(l, Parent(B, 3-m), 3-m, B, 0, PrPB)
      do x=1,3  ! shared parent
        do y=1,3  ! parent A
          PrXY(x,y) = OKA2P(Genos(l,A),x,y) * PrP(x,m) * PrPA(y) * &
             SUM(OKA2P(Genos(l,B),x,:) * PrPB)
        enddo
      enddo
    enddo
  else if (cat==4) then
    do m=1,2
      if (Parent(A,m)<0) then
        do n=1,2
          if (GpID(n, -parent(A,m), m) == B) then
            call ParProb(l, parent(A,m), m, A, -4, PrP(:,m))  
            call ParProb(l, parent(A,3-m), 3-m, A, 0, PrPA)
            call ParProb(l, GpID(3-n, -parent(A,m), m), 3-n, 0, 0, PrPB)
            call ParProb(l, B, n, 0, 0, PrP(:,3-m))
            do x=1,3  ! in-between parent
              do y=1,3  ! other parent of A
                PrXY(x,y) = SUM(OKA2P(Genos(l,A),x,:) *PrPA) *PrP(x,m)*&
                   SUM(AKA2P(x, y,:) * PrP(y,3-m) * PrPB)
              enddo
            enddo
          endif
        enddo
      endif
    enddo     
  endif
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine Upair

! #####################################################################

subroutine UClust(A, B, kA, kB, cat, Ai, Bj, LL)
use Global
implicit none

integer, intent(IN) :: A, B, kA, kB, cat, Ai, Bj
double precision, intent(OUT) :: LL
integer, allocatable, dimension(:) :: AA, BB
integer :: nA, nB, l,x,y, v, i, j, z,m, f, e, u, DoneA(maxSibSize), &
  catA(maxSibSize), catB(maxSibSize), GA(2), GB(2), g, Ei, AB(2*maxSibSize)
double precision :: PrL(nSnp,2), PrGA(3,2), PrGB(3,2), PrGGP(3), &
  PrUZ(3,3, 3,3,3,3,2), PrE(3), PrH(3), PrW(3)
integer, allocatable, dimension(:) :: UseEE, MateABpar, TypeEE
double precision, allocatable, dimension(:,:) :: PrEE                                                                                                                    

 call CalcCLL(-A, kA)
 call CalcCLL(-B, kB)
LL = 999D0

nA = nS(-A, kA)
allocate(AA(nA))
AA = SibID(1:nA, -A, kA)
GA = GpID(:, -A, kA)

nB = nS(-B, kB)
allocate(BB(nB))
BB = SibID(1:nB, -B, kB)
GB = GpID(:, -B, kB)

!============================================

allocate(UseEE(nA+nB))
allocate(TypeEE(nA+nB))
allocate(MateABpar(nA+nB))
allocate(PrEE(3, nA+nB))
UseEE = 0

if (kA==kB) then
  AB(1:nB) = BB
  AB((nB+1):(nB+nA)) = AA
  call FindEE(AB(1:(nB+nA)), nB, nA, kB, UseEE, MateABpar) 
  BB = AB(1:nB)
  AA = AB((nB+1):(nB+nA))
  TypeEE = 3-kB
  do i=1, nB  ! safety net
    if (UseEE(i) > nB) then
      UseEE(i) = 0  ! else use before store
    endif
  enddo
else if (kA/=kB) then
  call FindEE(BB, nB, 0, kB, UseEE(1:nB), MateABpar(1:nB))  ! may reorder BB
  call FindEE(AA, nA, 0, kA, UseEE((nB+1):(nB+nA)), MateABpar((nB+1):(nB+nA)))
  do i=1, nA
    if (UseEE(nB+i)/=0) then
      UseEE(nB+i) = nB + UseEE(nB+i)
    endif
  enddo
  TypeEE(1:nB) = 3-kB
  TypeEE((nB+1):(nB+nA)) = 3-kA
endif

!============================================                                      
catA = 0
catB = 0
do i = 1, nA
  do j = 1, nB
    if (kA /= kB) then
      if (Parent(AA(i), kB) == B) then
        catA(i) = 1
        UseEE(nB+i) = 0           
      endif
      if (Parent(BB(j), kA) == A) then
        catB(j) = 1
        UseEE(j) = 0          
      endif
    else if (kA == kB) then
      if (Parent(AA(i), 3-kA) == Parent(BB(j), 3-kB) .and. &
        Parent(BB(j), 3-kB)<0) then  
        catA(i) = 7
        catB(j) = 7
      endif
    endif
    do m=1,2
      if (GA(m) < 0) then
        if (GA(m) == Parent(AA(i), m) .and. m/=kA) then
          catA(i) = 2
        endif
        if (GA(m) == Parent(BB(j), m) .and. m/=kB) then
          catB(j) = 2
        endif
      endif
      if (GB(m) < 0) then
        if (GB(m) == Parent(AA(i), m) .and. m/=kA) then
          catA(i) = 3
        endif
        if (GB(m) == Parent(BB(j), m) .and. m/=kB) then
          catB(j) = 3  
        endif
      endif
    enddo
    if (Parent(AA(i), 3-kA) < 0) then
      if (GpID(kA, -Parent(AA(i), 3-kA), 3-kA)==A) then
        catA(i) = 6
      else if (GpID(kB, -Parent(AA(i), 3-kA), 3-kA)==B) then
        catA(i) = 8                                               
      endif
    endif
    if (Parent(BB(j), 3-kB) < 0) then
      if (GpID(kA, -Parent(BB(j),3-kB), 3-kB)==A) then
        catB(j) = 6
      else if (GpID(kB, -Parent(BB(j),3-kB), 3-kB)==B) then
        catB(j) = 8                                          
      endif 
    endif
  enddo
enddo

!==================================
if (cat==0 .and. ALL(catA==0) .and. ALL(CatB==0) .and. ALL(UseEE==0)) then
  if (A<0 .and.B>0) then
    LL = CLL(-A,kA) + Lind(B)
  else if (A>0 .and. B<0) then
    LL = Lind(A) + CLL(-B,kB)
  else if (A<0 .and. B<0) then
    LL = CLL(-A,kA) + CLL(-B,kB)
  endif
  return
endif

!==================================

PrL = 0D0
do l=1, nSnp
  PrUZ = 0
  
  do m=1, 2
    if ((ANY(catA==2) .and. m/=kA) .or. (ANY(catB==2) .and. m/=kB)) then
      call ParProb(l, GA(m), m, -1, 0, PrGA(:, m))
    else
      call ParProb(l, GA(m), m, 0, 0, PrGA(:, m))
    endif
    if ((ANY(catA==3) .and. m/=kA) .or. (ANY(catB==3) .and. m/=kB)) then
      call ParProb(l, GB(m), m, -1, 0, PrGB(:, m))
    else
      call ParProb(l, GB(m), m, 0, 0, PrGB(:, m))
    endif
  enddo
  if (cat==4) then
    do m=1,2
      if (GA(m)<0) then
        if (GpID(kB, -GA(m), m) == B) then
          call ParProb(l, GpID(3-kB, -GA(m), m), 3-kB, 0, 0, PrGGP) 
        endif
      endif
    enddo
  endif
  
  ! == grandparents ==
  do x=1,3 
    do y=1,3
      do u=1,3  ! GP A, kB
        do z=1,3  ! GP A, 3-kB
          do v=1,3  ! GP B, kB
            if (cat == 1) then
             if (kA==kB .and. GA(3-kB)==GB(3-kB) .and. GA(3-kB)/=0) then 
                PrUZ(x,y,y,z,v,z,1) = AKA2P(x,y,z) * AKA2P(y,v,z) *&
                 PrGA(z,3-kB) * PrGB(v,kB)
              else
                PrUZ(x,y,y,z,v,:,1) = AKA2P(x,y,z) * AKA2P(y,v,:) * &
                  PrGA(z,3-kB) * PrGB(v,kB) * PrGB(:,3-kB)
              endif
            else if (cat==2) then
              PrUZ(x,y,u,z,u,z,1) = AKA2P(x,u,z) * AKA2P(y,u,z) *&
               PrGA(u,kB) * PrGA(z,3-kB)
            else if (cat==3) then
              do m=1,2
                if (GA(m)/=0 .and. GA(m) == GB(m)) then
                  if (m==kB) then
                    PrUZ(x,y,u,z,u,:,1) = AKA2P(x,u,z) * AKA2P(y,u,:) *&
                       PrGA(u,m) * PrGA(z,3-m) * PrGB(:,3-m)   
                  else
                    PrUZ(x,y,u,z,v,z,1) = AKA2P(x,u,z) * AKA2P(y,v,z) *&
                       PrGA(u,3-m) * PrGA(z,m) * PrGB(v,3-m)
                  endif
                endif
              enddo
            else if (cat==4) then 
              do m=1,2
                if (GA(m)<0) then
                  if (GpID(kB, -GA(m), m) == B) then
                    if (m==kB) then
                      PrUZ(x,y,u,z,v,:,1) =AKA2P(x,u,z) *&
                       SUM(AKA2P(u,y,:) *PrGGP) *PrGA(z,3-kB) *&
                       AKA2P(y,v,:) * PrGB(v,kB) * PrGB(:, 3-kB) 
                    else
                      PrUZ(x,y,u,z,v,:,1) = AKA2P(x,u,z) *&
                     SUM(AKA2P(z,y,:) *PrGGP) *PrGA(u,kB) *AKA2P(y,v,:)&
                     * PrGB(v,kB) * PrGB(:, 3-kB)
                    endif
                  endif
                endif
              enddo
            else
              PrUZ(x,y,u,z,v,:,1) = AKA2P(x,u,z) * AKA2P(y,v,:) * &
                PrGA(u,kB) * PrGA(z,3-kB) * PrGB(v,kB) * PrGB(:, 3-kB)
            endif
            PrUZ(x,y,u,z,v,:,2) = PrUZ(x,y,u,z,v,:,1)
          enddo
        enddo
      enddo
    enddo
  enddo
   
  ! == siblings ==   
  if (ALL(catA==0) .and. ALL(catB==0) .and. Ai==0 .and. Bj==0 .and. &
    ALL(UseEE==0)) then       
    do x=1,3 
      do y=1,3
        PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * XPr(1,x,l,-A,kA) *&
         XPr(1,y,l,-B,kB)
      enddo  ! TODO: needs special for cat<4 ?
    enddo
  
  else

  do x=1,3  ! SA
    doneA = 0
    do y=1,3  ! SB
      PrEE = 0            
      do j=1, nB
        if (nFS(BB(j))==0) cycle
        if (catB(j)==1 .or. catB(j)==2 .or. catB(j)==3) then
          PrE = 1
        else if (catB(j)==6) then  
          call ParProb(l, GpID(3-kA,-Parent(BB(j),3-kB),3-kB),3-kA, 0, 0, PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,x,:) * PrH) 
          enddo
        else if (catB(j)==8) then  
          call ParProb(l, GpID(3-kB,-Parent(BB(j),3-kB),3-kB),3-kB, 0, 0, PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,y,:) * PrH) 
          enddo
        else if (UseEE(j)/=0) then
          call ParProb(l, MateABpar(j), 3-TypeEE(j), 0,0,PrH)
          do e=1,3
            do u=1, 3
              PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,UseEE(j)) * PrH)
            enddo
            PrE(e) = SUM(PrW)
          enddo
          PrE = PrE/SUM(PrE)
        else
          call ParProb(l, Parent(BB(j),3-kB), 3-kB, -1, 0, PrE)
        endif
        
        if (Parent(BB(j),3-kB) < 0 .and. catB(j)/=1) then 
          do e=1,3
            do g=1, nS(-Parent(BB(j),3-kB), 3-kB)
              Ei = SibID(g, -Parent(BB(j),3-kB), 3-kB)
              if (nFS(Ei) == 0) cycle
              if (Parent(Ei, kB) == B) cycle  
              if (Parent(Ei, kA) == A) cycle 
              call ParProb(l, Parent(Ei, kB), kB, Ei, -1, PrH) 
              do i=1, nFS(Ei)
                if (Genos(l,FSID(i, Ei))==-9) cycle
                PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
              enddo
              PrE(e) = PrE(e) * SUM(PrH)
            enddo
          enddo
        endif
        
        if (Bj/=0 .or. catB(j)==7 .or. (catB(j)==1 .and. Ai/=0)) then
          do f=1, nFS(BB(j))
            if (Bj==0 .or. FSID(f, BB(j))==Bj) cycle
            if (Parent(BB(j),kA)==A .and. (Ai==0 .or. FSID(f, BB(j))==Ai)) cycle
            if (Genos(l,FSID(f, BB(j)))==-9) cycle
            PrE = PrE * OKA2P(Genos(l,FSID(f,BB(j))), y, :)
          enddo
        endif
        
        if (catB(j)==7 .and. Ai/=0) then 
          do i=1,nA
            if (Parent(AA(i), kB) == B) cycle
            do f=1, nFS(AA(i))
              if (FSID(f, AA(i))==Ai) cycle
              if (Genos(l,FSID(f, AA(i)))==-9) cycle
              PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
            enddo
          enddo
        endif
        
        if (catB(j)==1) then  ! Parent(BB(j), 3-kB)==PA
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * PrE(x)
        else if (CatB(j)==2) then
          do z=1,3
            PrUZ(x,y,:,z,:,:,1) = PrUZ(x,y,:,z,:,:,1) * PrE(z)   
          enddo
        else if (CatB(j)==3) then
          do z=1,3
            PrUZ(x,y,:,:,:,z,1) = PrUZ(x,y,:,:,:,z,1) * PrE(z)   
          enddo                       
        else 
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * SUM(PrE)
        endif
        
        do f=1, nFS(BB(j)) ! includes some AA if cat=1 
          if (Bj==0 .or. FSID(f, BB(j))==Bj .or. &
            (Parent(BB(j),kA)==A .and. (Ai==0 .or. FSID(f, BB(j))==Ai))) then
            if (Genos(l,FSID(f, BB(j)))==-9) cycle
            PrE = PrE * OKA2P(Genos(l,FSID(f, BB(j))), y, :)
          endif
        enddo

        if (catB(j)==7) then 
          do i=1,nA
            if (Parent(AA(i), 3-kB) /= Parent(BB(j), 3-kB)) cycle
            if (Parent(AA(i), kB) == B) cycle
            if (Ai==0 .and. nFS(AA(i))==0) cycle
            do f=1, nFS(AA(i))
              if (Ai/=0 .and. FSID(f, AA(i))/=Ai) cycle
              if (Genos(l,FSID(f, AA(i)))==-9) cycle
              PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
              DoneA(i) = 1
            enddo
          enddo
        endif
        
        if (catB(j)==1) then  ! Parent(BB(j), 3-kB)==PA
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * PrE(x)
        else if (CatB(j)==2) then
          do z=1,3
            PrUZ(x,y,:,z,:,:,2) = PrUZ(x,y,:,z,:,:,2) * PrE(z)   
          enddo
        else if (CatB(j)==3) then
          do z=1,3
            PrUZ(x,y,:,:,:,z,2) = PrUZ(x,y,:,:,:,z,2) * PrE(z)   
          enddo                       
        else 
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * SUM(PrE)
        endif
        PrEE(:,j) = PrE
      enddo  ! B_j
    
      do i=1, nA
        if (DoneA(i)==1) cycle
        if (nFS(AA(i))==0) cycle
        if (Parent(AA(i),kB)==B) cycle
        if (catA(i)>1 .and. catA(i)<4) then  ! catA==1 already done
          PrE = 1
        else if (catA(i)==6) then ! Parent(AA(i), 3-kA) <0
          call ParProb(l, GpID(3-kA,-Parent(AA(i),3-kA), 3-kA), 3-kA, 0, 0, PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,x,:) * PrH) 
          enddo
        else if (catA(i)==8) then ! Parent(AA(i), 3-kA) <0
          call ParProb(l, GpID(3-kB,-Parent(AA(i),3-kA),3-kA),3-kB,0,0,PrH)
          do e=1,3
            PrE(e) = SUM(AKA2P(e,y,:) * PrH) 
          enddo
        else if (UseEE(nB+i)/=0) then
          call ParProb(l, MateABpar(nB+i), 3-TypeEE(nB+i), 0,0,PrH)
          do e=1,3
            do u=1, 3
              PrW(u) = SUM(AKA2P(e,u,:) * PrEE(u,UseEE(nB+i)) * PrH)
            enddo
            PrE(e) = SUM(PrW)
          enddo
          PrE = PrE/SUM(PrE) 
        else    
          call ParProb(l, Parent(AA(i), 3-kA), 3-kA, -1, 0, PrE)   
        endif
        
        if (Parent(AA(i), 3-kA) < 0) then 
          do e=1,3
            do g=1, nS(-Parent(AA(i), 3-kA), 3-kA)
              Ei = SibID(g, -Parent(AA(i), 3-kA), 3-kA)
              if (nFS(Ei) == 0) cycle
              if (Parent(Ei, kA) == A) cycle 
              if (Parent(Ei, kB) == B) cycle
              call ParProb(l, Parent(Ei, kA), kA, Ei, -1, PrH) 
              do j=1, nFS(Ei)                  
                if (Genos(l,FSID(j, Ei))==-9) cycle
                PrH = PrH * OKA2P(Genos(l,FSID(j,Ei)), :, e)
              enddo
              PrE(e) = PrE(e) * SUM(PrH)
            enddo
          enddo
        endif
        
        if (Ai/=0) then
          do f=1, nFS(AA(i))
            if (FSID(f, AA(i))==Ai) cycle
            if (Genos(l,FSID(f, AA(i)))==-9) cycle
            PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
          enddo
        endif
        
        if (catA(i)==2) then
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,z,:,:,1) = PrUZ(x,y,:,z,:,:,1) * PrE(z)
            else
              PrUZ(x,y,z,:,:,:,1) = PrUZ(x,y,z,:,:,:,1) * PrE(z)         
            endif
          enddo
        else if (catA(i)==3) then  
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,:,:,z,1) = PrUZ(x,y,:,:,:,z,1) * PrE(z)
            else
              PrUZ(x,y,:,:,z,:,1) = PrUZ(x,y,:,:,z,:,1) * PrE(z)
            endif
          enddo
        else
          PrUZ(x,y,:,:,:,:,1) = PrUZ(x,y,:,:,:,:,1) * SUM(PrE)    
        endif
        
        do f=1, nFS(AA(i)) 
          if (Ai/=0 .and. FSID(f, AA(i))/=Ai) cycle 
          DoneA(i)=2                                           
          if (Genos(l,FSID(f, AA(i)))==-9) cycle
          PrE = PrE * OKA2P(Genos(l,FSID(f,AA(i))), x, :)
        enddo
        
        if (catA(i)==2) then
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,z,:,:,2) = PrUZ(x,y,:,z,:,:,2) * PrE(z)
            else
              PrUZ(x,y,z,:,:,:,2) = PrUZ(x,y,z,:,:,:,2) * PrE(z)         
            endif
          enddo
        else if (catA(i)==3) then  
          do z=1,3
            if (kA==kB) then
              PrUZ(x,y,:,:,:,z,2) = PrUZ(x,y,:,:,:,z,2) * PrE(z)
            else
              PrUZ(x,y,:,:,z,:,2) = PrUZ(x,y,:,:,z,:,2) * PrE(z)
            endif
          enddo
        else
          PrUZ(x,y,:,:,:,:,2) = PrUZ(x,y,:,:,:,:,2) * SUM(PrE)    
        endif
        PrEE(:,nB+i) = PrE     
      enddo  ! i
    enddo  ! x
  enddo  ! y
  endif
  do f=1,2
    PrL(l,f) = LOG10(SUM(PrUZ(:,:,:,:,:,:,f)))
  enddo
enddo
LL = SUM(PrL(:,2)) - SUM(PrL(:,1))

deallocate(AA)
deallocate(BB)
deallocate(UseEE)
deallocate(TypeEE)
deallocate(PrEE)
deallocate(MateABpar)

end subroutine UClust

! #####################################################################

subroutine FindEE(AB, nA, nB, k, UseEE, MatePar)  ! find PO pairs among mates
use Global
use qsort_c_module
implicit none

integer, intent(IN) :: nA, nB, k
integer, intent(INOUT) :: AB(nA+nB)
integer, intent(OUT) :: UseEE(nA+nB), MatePar(nA+nB)
integer :: i, j,x, nAB(2), Mate(MAX(nA,nB), 2), GGK(2), Order(2*maxSibSize),&
  ABM(MAX(nA,nB),2), UseM(MAX(nA,nB),2), MateI
logical :: reorder, OrderAgain
double precision :: EEtmp(2*maxSibSize)

UseEE = 0
MatePar = 0

ABM = 0
ABM(1:nA, 1) = AB(1:nA)
ABM(1:nB, 2) = AB((nA+1):(nA+nB))
nAB = (/nA, nB/)
Mate = 0
do x=1,2
  do i=1, nAB(x)
    if (nFS(ABM(i,x))==0 .and. nAB(x)>1)  cycle
    Mate(i,x) = Parent(ABM(i,x), 3-k)
  enddo
enddo
if ((nAB(1)==1 .and. nAB(2)<2) .or. COUNT(Mate < 0) < 2) return

GGK = 0
do x=1,2
  if (ABM(1,x)==0)  cycle
  if (Parent(ABM(1,x),k) < 0) then  ! else not called?
    GGK(x) = GpID(3-k, -Parent(ABM(1,x),k), k)   
  endif
enddo

! re-order AA and BB, so that PrE calculated before used
UseM = 0
reorder = .FALSE.
do x=1,2
  if (nAB(x)<=1) cycle
  do i=1, nAB(x)
    if (Mate(i,x) < 0) then 
      if (GpID(3-k, -Mate(i,x), 3-k) < 0 .and. &
       .not. ANY(GGK == GpID(3-k, -Mate(i,x), 3-k))) then
        do j=1, nAB(x)
          if (Mate(j,x) == GpID(3-k, -Mate(i,x), 3-k)) then
            UseM(i,x) = j
            if (j > i) reorder = .TRUE.
            exit
          endif
        enddo
      endif
    endif
  enddo
  
  if (reorder) then
    EEtmp(1:nAB(x)) = UseM(1:nAB(x),x)
    Order = (/ (i, i=1, nAB(x), 1) /)
    call QsortC(EEtmp(1:nAB(x)), Order(1:nAB(x)))
    ABM(1:nAB(x),x) = ABM(Order(1:nAB(x)), x)
    UseM(1:nAB(x),x) = UseM(Order(1:nAB(x)),x)
    OrderAgain = .FALSE.
    do i=1, nAB(x)
      if (UseM(i,x) /= 0) then
        do j=1, nAB(x)
          if (UseM(i,x) == Order(j)) then
            UseM(i,x) = j
            if (j>i)  OrderAgain = .TRUE.
            exit
          endif
        enddo
      endif
    enddo
    if (OrderAgain) then
      EEtmp(1:nAB(x)) = UseM(1:nAB(x),x)
      Order = (/ (i, i=1, nAB(x), 1) /)
      call QsortC(EEtmp(1:nAB(x)), Order(1:nAB(x)))
      ABM(1:nAB(x),x) = ABM(Order(1:nAB(x)), x)
    endif
  endif
enddo

AB = 0
AB(1:nA) = ABM(1:nA, 1)
AB((nA+1):(nA+nB)) = ABM(1:nB, 2)

do i=1, nA+nB
  if (nFS(AB(i))==0 .and. ((i<=nA .and. nA>1) .or. (i>nA .and. nB>1)))  cycle
  MateI = Parent(AB(i), 3-k)
  if (MateI < 0) then 
    if (GpID(3-k, -MateI, 3-k) < 0 .and. &
     .not. ANY(GGK == GpID(3-k, -MateI, 3-k))) then
      do j=1, i
      if (nFS(AB(j))==0 .and. ((j<=nA .and. nA>1) .or. (j>nA .and. nB>1)))  cycle
        if (Parent(AB(j), 3-k) == GpID(3-k, -MateI, 3-k)) then
          UseEE(i) = j
          MatePar(i) = GpID(k, -MateI, 3-k)
          exit
        endif
      enddo
    endif
  endif
enddo

end subroutine FindEE

! #####################################################################

subroutine AddSib(A, SB, k, LL)  
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x, Bj, m, f, AncB(2,mxA), OpPar, j,v,  Ei, i, z
double precision :: PrL(nSnp), PrX(3), PrXb(3,2), PrY(3), PrZ(3), PrE(3)
logical :: Inbr, AllFS                                                                      

LL = 999D0
if (Parent(A,k)==-SB) then
  LL = 888
else if (Parent(A,k)/=0) then
  LL = 777
endif
if (LL/=999) return

do f=1, nS(SB,k)
  Bj = SibID(f, SB, k)
  if (Parent(A, 3-k) /= 0) then
    if (Parent(Bj, 3-k) == Parent(A, 3-k)) then
      LL = 777  ! use addFS() instead
    endif
  endif
  if (AgeDiff(A,Bj)==999) cycle
  if (AgePriorM(ABS(AgeDiff(A, Bj))+1,k)==0.0) then   
    LL=777
  endif 
enddo
if (LL/=999) return

 call GetAncest(-SB, k, AncB)
if (ANY(AncB == A)) then  ! A>0
  LL = 777
else if (BY(A)/=-999) then
  do x=3, mxA
    do m=1,2
      if (AncB(m, x) > 0) then
        if (AgeDiff(A, AncB(m, x)) <=0) then  ! A older than anc
          LL = 777
        endif
      endif
    enddo
  enddo
endif
if (LL == 777) return

Inbr = .FALSE.
if (Parent(A,3-k) < 0) then
  if (Parent(A,3-k) == GpID(3-k, SB, k)) then
    Inbr = .TRUE.  ! inbreeding loop created
  else if (GpID(k,-Parent(A,3-k),3-k) == -SB) then
    Inbr = .TRUE.
  endif
endif
do f=1, nS(SB,k)
  Bj = SibID(f, SB, k)
  if (Parent(A,3-k) == Bj)  Inbr = .TRUE.
  if (Parent(Bj,3-k) == A)  Inbr = .TRUE.
  if (Parent(Bj,3-k) == GpID(3-k,SB,k) .and. Parent(Bj,3-k)/=0) Inbr = .TRUE.
enddo

AllFS = .FALSE.
 call getFSpar(SB, k, .TRUE., OpPar)
if (OpPar < 0 .or. ns(SB,k)==0) then
  AllFS = .TRUE.
endif

if (.not. Inbr .and. .not. AllFS) then
  PrL = 0D0
  do l=1,nSnp
    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrY)          
    if (nS(SB,k)==0) then
      PrX = XPr(2,:,l, SB,k)
    else if (.not. AllFS) then
      PrX = XPr(3,:,l, SB,k)
    endif
    do x=1,3
      if (Genos(l,A) /= -9) then
        PrX(x) = PrX(x) * SUM(OKA2P(Genos(l,A), x, :) * PrY)
      endif
    enddo
    PrL(l) = LOG10(SUM(PrX))   
  enddo
  LL = SUM(PrL)

else if (.not. Inbr) then
  PrL = 0D0
  do l=1,nSnp 
    call ParProb(l,-SB,k,-1,0,PrXb(:,1))  ! GPs
    PrXb(:,2) = PrXb(:,1)
    do x=1,3
      do j=1, nS(SB,k)
        Bj = SibID(j,SB,k)
        if (nFS(Bj)==0) cycle
        call ParProb(l, Parent(Bj,3-k), 3-k, -1, 0, PrZ)
        do z=1,3
          if (Parent(Bj,3-k)<0) then
            do v=1, nS(-Parent(Bj, 3-k), 3-k)
              Ei = SibID(v, -Parent(Bj, 3-k), 3-k)  
              if (NFS(Ei) == 0) cycle
              if (Parent(Ei, k) == -SB) cycle
              call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)        
              do i=1, nFS(Ei) 
                if (Genos(l,FSID(i, Ei))/=-9) then
                  PrE = PrE * OKA2P(Genos(l,FSID(i,Ei)), :, z)
                endif
              enddo
              PrZ(z) = PrZ(z) * SUM(PrE)  
            enddo  
          endif
        enddo
        PrXb(x,1) = PrXb(x,1) * SUM(PrZ)
          
        do v=1, nFS(Bj)
          if (Genos(l, FSID(v,Bj))==-9) cycle
          PrZ = PrZ * OKA2P(Genos(l,FSID(v,Bj)),x,:)
        enddo
        PrXb(x,2) = PrXb(x,2) * SUM(PrZ)
      enddo
    enddo

    call ParProb(l, Parent(A,3-k), 3-k, A, 0, PrY)
    do x=1,3
      if (Genos(l,A) /= -9) then
        PrXb(x,2) = PrXb(x,2) * SUM(OKA2P(Genos(l,A), x, :) * PrY)
      endif
    enddo
    PrL(l) = LOG10(SUM(PrXb(:,2))) - LOG10(SUM(PrXb(:,1)))  
  enddo
  LL = SUM(PrL)
else
  call DoAdd(A, SB, k)
  LL = CLL(SB,k)
  call RemoveSib(A, SB, k)
endif

end subroutine AddSib

! #####################################################################

subroutine AddSibInbr(A,SB,k,LL)
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL(3)
integer :: l, x, y, GA, GG, Par, i, u
double precision :: PrL(nSnp,3), PrXY(3,3), PrZ(3), PrPA(3), LLtmp(3), &
  ALR(3), LLU(4), PrXYU(3,3,3), PrLU(nSnp,3)
logical :: maybe(3)

! 1: Par(Parent(A,3-k),k)=SB, 2: Parent(A,3-k)=GpID(3-k,SB,k)
! 3: as 1, A FS of B's (PA == DB)

LL = 999D0
maybe = .TRUE.
! TODO: check all ancestors. 
if (Parent(A,3-k)>0) then
  if (Parent(Parent(A,3-k),k)/=0 .and. Parent(Parent(A,3-k),k)/=-SB) then
    maybe(1) = .FALSE.
  endif
  if (Parent(Parent(A,3-k),k)==-SB) then
    maybe(2) = .FALSE.
  endif
  GA = Parent(Parent(A,3-k), 3-k)
else if (Parent(A,3-k)<0) then
  if (GpID(k,-Parent(A,3-k),3-k)/=0 .and. GpID(k,-Parent(A,3-k),3-k)/=-SB) then
    maybe(1) = .FALSE.
  endif
  GA = GpID(3-k,-Parent(A,3-k),3-k)
else
  GA = 0
endif
GG = GpID(k,SB,k)
if (GpID(3-k,SB,k)/=0 .and. GpID(3-k,SB,k)/=Parent(A,3-k)) then
  maybe(2) = .FALSE.
else if (ANY(SibID(1:ns(SB,k),SB,k)==Parent(A,3-k))) then
  maybe(2) = .FALSE.
endif
if (.not. ANY(maybe(1:2))) then
  LL = 777
  return
endif

Par = 0
if (maybe(1)) then
  call getFSpar(SB, k, .TRUE., Par)
  if (Par==0) then
    maybe(3) = .FALSE.
  else if (Parent(A,3-k)/=Par .and. Parent(A,3-k)/=0) then
    maybe(3) = .FALSE.
  else if (Par>0) then
    if (Parent(Par, k)/=SB .and. Parent(Par, k)/=0) then
      maybe(3) = .FALSE.
    endif
  else if (Par<0) then
    if (GpID(k, -Par, 3-k)/=SB .and. GpID(k, -Par, 3-k)/=0) then
      maybe(3) = .FALSE.
    endif
  endif
  if (maybe(3) .and. GA==0) then
    if (Par>0)  GA = Parent(Par, 3-k)
    if (Par<0)  GA = GpID(3-k, -Par, 3-k)
  endif
endif

call CalcAgeLR(Parent(A,3-k),3-k, SB,k, 0,1, .TRUE., ALR(1))
call CalcAgeLR(SB,k, Parent(A,3-k),3-k,  0,1, .TRUE., ALR(2))
call CalcAgeLR(Par,3-k, SB,k, 0,1, .TRUE., ALR(3))
do x=1,3
  if (ALR(x) == 777 .or. ALR(x) < 3*TF)  maybe(x) = .FALSE.
enddo

if (.not. ANY(maybe)) then
  LL = 777
  return
endif

PrL = 0D0
PrLU = 0D0
do l=1,nSnp
  if (maybe(1)) then
    if (Parent(A,3-k)>0) then
      call ParProb(l, GA, 3-k, Parent(A,3-k), 0, PrZ) 
      call ParProb(l, Parent(A,3-k),3-k,0,0,PrPA)
    else
      call ParProb(l, GA, 3-k, 0, 0, PrZ) 
      call ParProb(l, Parent(A,3-k),3-k,A,-4,PrPA)
    endif
    if (Parent(A,3-k)==0)   PrPA = 1
    do x=1,3
      do y=1,3
        PrXY(x,y) = XPr(3,x,l, SB,k) * PrPA(y) * SUM(AKA2P(y,x,:) * PrZ)
        if (Genos(l,A)/=-9) then
          PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
        endif
        do u=1,3
          PrXYU(x,y,u) = XPr(3,x,l, SB,k) * PrPA(y) * SUM(AKA2P(y,u,:) * PrZ)
          if (Genos(l,A)/=-9) then
            PrXYU(x,y,u) = PrXYU(x,y,u) * OKA2P(Genos(l,A), u, y)
          endif
        enddo
      enddo
    enddo
    PrL(l,1) = LOG10(SUM(PrXY))   ! Parent(A,3-k) offspring of SB 
    PrLU(l,1) = LOG10(SUM(PrXYU)) 
  endif
  
  !===
  if(maybe(3)) then
    do i=1, ns(SB, k)
      if (nFS(SibID(i,SB,k))==0 .or. Parent(SibID(i,SB,k),3-k)/=Par)  cycle
      call ParProb(l, Par,3-k,SibID(i,SB,k),-5,PrPA)  ! exclude both GPs & Bi & FS of Bi
    enddo
    if (Par < 0) then
      if (Parent(A,3-k)==Par .and. Genos(l,A)/=-9) then
        PrPA = PrPA/OKAP(Genos(l,A),:,l)
        PrPA = PrPA/SUM(PrPA)
      endif
      call ParProb(l, GA, 3-k, 0, 0, PrZ) 
    else
      call ParProb(l, GA, 3-k, Par, 0, PrZ) 
    endif
    do x=1,3
      do y=1,3
        PrXY(x,y) = PrPA(y) * SUM(AKA2P(y,x,:) * PrZ) * XPr(2,x,l, SB,k)  !Xpr(2,) = GP's
        do i=1, ns(SB, k)
           if (Genos(l,SibID(i,SB,k))/=-9) then
            PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,SibID(i,SB,k)), x, y)
          endif
        enddo
        if (Genos(l,A)/=-9) then
          PrXYU(x,y,:) = PrXYU(x,y,:) * OKAP(Genos(l,A), :, y)
          PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
        endif
      enddo
    enddo
    PrL(l,3) = LOG10(SUM(PrXY)) 
    PrLU(l,3) = LOG10(SUM(PrXYU)) 
  endif
  
  !===
  if (maybe(2)) then
    call ParProb(l, GG, k, 0, 0, PrZ)
    call ParProb(l, Parent(A,3-k),3-k, A,0,PrPA)
    do x=1,3
      do y=1,3
        PrXY(x,y) = XPr(1,x,l, SB,k) * PrPA(y) * SUM(AKA2P(x,y,:) * PrZ)
        if (Genos(l,A)/=-9) then
          PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,A), x, y)
        endif
        do u=1,3
          PrXYU(x,y,u) = XPr(1,x,l, SB,k) * PrPA(y) * SUM(AKA2P(u,y,:) * PrZ)
          if (Genos(l,A)/=-9) then
            PrXYU(x,y,u) = PrXYU(x,y,u) * OKA2P(Genos(l,A), u, y)
          endif
        enddo
      enddo
    enddo
    PrL(l,2) = LOG10(SUM(PrXY))   ! SB offspring of Parent(A,3-k)
    PrLU(l,2) = LOG10(SUM(PrXYU)) 
  endif
enddo
LLtmp = SUM(PrL, dim=1)
LLU(1:3) = SUM(PrLU, dim=1)
call CalcU(A,k, -SB, k, LLU(4))   ! unrelated & A non-inbred
if (maybe(1) .and. Parent(A,3-k)>0) then
  LLtmp(1) = LLtmp(1) - Lind(Parent(A,3-k))
endif

LL = 777D0
do x=1,3
  if (.not. maybe(x))  cycle
  if (LLU(x) > LLU(4)) then
    LL(x) = LLtmp(x) - LLU(x) + LLU(4)
  else
    LL(x) = LLtmp(x)
  endif
enddo

end subroutine AddSibInbr

! #####################################################################
subroutine MergeSibs(SA, SB, k, LL)  
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: LL
integer :: l,x,y, r,v, Bj, Ai, i, G(2), m, Ei, e,f, j, &
  AncA(2,mxA), AncB(2,mxA), catG, catB(ns(SB,k)), catA(ns(SA,k)), &
  nAB(2), AB(2,maxSibsize), catGG(2), z, GGP(2)
double precision :: PrL(nSnp, 2), PrG(3,2), PrXY(3,3,3,2), PrE(3), PrH(3)

LL = 999D0
G = 0  
do m=1,2  
  if (GpID(m,SA,k) /= 0) then
    if (GpID(m,SB,k) /= 0 .and. GpID(m,SA,k)/=GpID(m,SB,k)) then
      LL = 777  ! incompatible grandparents
    else
      G(m) = GpID(m,SA,k)  ! including if GP(B) is dummy
    endif
  else if (GpID(m,SA,k) == 0) then
    G(m) = GpID(m,SB,k)
  endif
enddo
if (GpID(k, SA,k)==-SB .or. GpID(k, SB, k)==-SA) then
  LL = 777
endif
if (LL==777) return
 call GetAncest(-SA, k, AncA)
 call GetAncest(-SB, k, AncB)
if (ANY(AncA(k, 3:mxA) == -SB)) LL = 777
if (ANY(AncB(k, 3:mxA) == -SA)) LL = 777
if (LL==777) return 

nAB(1) = ns(SA,k)
nAB(2) = ns(SB, k)
AB(1, 1:ns(SA,k)) = SibID(1:ns(SA,k), SA,k)
AB(2, 1:ns(SB,k)) = SibID(1:ns(SB,k), SB, k)

catG = 0
catGG = 0
GGP = 0                                            
if (ANY(G/=0)) then
  do j=2,1,-1
    do i=1,nAB(j)
      if (nFS(AB(j,i))==0) cycle
      if (Parent(AB(j,i), 3-k)==0) cycle
      if (Parent(AB(j,i), 3-k) == G(3-k)) then
       if (catG==0) then
          catG = AB(j,i)
          exit  
        endif  
      endif
      do m=1,2
        if (G(m)==0) cycle
        if (Parent(AB(j,i), 3-k) > 0) then
          if (Parent(Parent(AB(j,i),3-k), m) == G(m)) then
            catGG(m) = AB(j,i)
            GGP(3-m) = Parent(Parent(AB(j,i),3-k), 3-m)
          endif
        else if (Parent(AB(j,i), 3-k) < 0) then
          if (GpID(m, -Parent(AB(j,i), 3-k),3-k) == G(m)) then
            catGG(m) = AB(j,i)
            GGP(3-m) = GpID(3-m, -Parent(AB(j,i), 3-k),3-k)
          endif
        endif
      enddo
    enddo
    if (catG /= 0) exit                           
  enddo
endif

catB = 0
catA = 0
do r = 1, nS(SB, k)
  Bj = SibID(r, SB, k) 
  if (nFS(Bj)==0 .or. Parent(Bj,3-k)==0) cycle
  do v=1,nS(SA,k)
    Ai = SibID(v, SA, k)   
    if (nFS(Ai)==0) cycle
    if (Parent(Ai,3-k) == Parent(Bj,3-k)) then
      catA(v) = r
      catB(r) = 1 
    endif
  enddo
enddo

PrL = 0D0
do l=1,nSnp
  do m=1,2
    if (m/=k .and. catG/=0) then
      call ParProb(l, G(m), m, -1, 0, PrG(:,m))
    else if (catGG(m)/=0) then  
      if (Parent(catGG(m),3-k) > 0) then
        call ParProb(l, G(m), m, Parent(catGG(m),3-k), 0, PrG(:,m))
      else
        call ParProb(l, G(m), m, 0, 0, PrG(:,m))
      endif
    else
      call ParProb(l, G(m), m, 0, 0, PrG(:,m))
    endif
  enddo
  do x=1,3
    do y=1,3
      do z=1,3
        PrXY(x,y,z,:) = AKA2P(x, y, z) * PrG(y,3-k) * PrG(z,k)
      enddo
    enddo
  enddo
  
  do z=1,3  !  =y if catGG(3-k)/=0
    do y=1,3
      if ((y>1 .or. z>1) .and. ALL(catGG==0)) cycle
  do x=1,3
    do j=1,2
      do r=1, nAB(j)
        if (j==2) then
          if (catB(r)==1) cycle ! done as FS of an A
        endif
        Ai = AB(j,r)
        if (NFS(Ai) == 0) cycle
        if (catG==Ai) then
          PrE = 1
        else if (ANY(catGG == Ai)) then
          if (ALL(catGG == Ai)) then
            PrE = AKA2P(:,z,y)
          else
            do m=1,2
              if (catGG(m)==Ai) then
                if (Parent(Ai,3-k)>0) then
                  call ParProb(l, GGP(3-m), 3-m, Parent(Ai,3-k),0,PrH)
                else
                  call ParProb(l, GGP(3-m), 3-m, 0,0,PrH)
                endif
                do e=1,3
                  if (m==k) then
                    PrE(e) = SUM(AKA2P(e, z, :) * PrH)
                  else
                    PrE(e) = SUM(AKA2P(e, y, :) * PrH)
                  endif
                enddo
              endif
            enddo
          endif
          if (Parent(Ai,3-k)>0) then
            if (Genos(l, Parent(Ai,3-k))/=-9) then
              PrE = PrE * OcA(Genos(l, Parent(Ai,3-k)), :)
            endif
          endif
!          PrE = PrE/SUM(PrE)                                     
        else
          call ParProb(l, Parent(Ai,3-k), 3-k, -1, 0, PrE)
        endif

        if (Parent(Ai,3-k) < 0) then 
          do e=1,3
            do f=1, nS(-Parent(Ai,3-k), 3-k)
              Ei = SibID(f, -Parent(Ai,3-k), 3-k)
              if (nFS(Ei) == 0) cycle
              if (Parent(Ei, k)==-SB .or. Parent(Ei,k)==-SA) cycle  
              if (catGG(3-k)>0) then
                if (Parent(catGG(3-k),3-k) == Ei)  cycle
              endif    
              call ParProb(l, Parent(Ei, k), k, Ei, -1, PrH) 
              do i=1, nFS(Ei)
                if (Genos(l,FSID(i, Ei))==-9) cycle
                PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
              enddo
              PrE(e) = PrE(e) * SUM(PrH)
            enddo
          enddo
        endif
        
        if (SUM(PrE)<3) then
          if (catG==Ai) then
            if (ANY(catGG/=0)) then
              PrXY(x,y,z,1) = PrXY(x,y,z,1) * PrE(y)
            else
              do e=1,3
                PrXY(x,e,:,1) = PrXY(x,e,:,1) * PrE(e)
              enddo
            endif
          else if (ANY(catGG/=0)) then
            PrXY(x,y,z,1) = PrXY(x,y,z,1) * SUM(PrE)
          else 
            PrXY(x,:,:,1) = PrXY(x,:,:,1) * SUM(PrE)
          endif                                 
        endif

        do f=1, nFS(Ai)
          if (Genos(l,FSID(f, Ai))==-9) cycle
          PrE = PrE * OKA2P(Genos(l,FSID(f, Ai)), x, :)
        enddo

        if (j==1) then
          if (catA(r)/=0) then
            Bj = AB(2, catA(r))
            do f=1, nFS(Bj) 
              if (Genos(l,FSID(f, Bj))==-9) cycle
              PrE = PrE * OKA2P(Genos(l,FSID(f, Bj)), x, :)
            enddo
          endif
        endif

        if (SUM(PrE)<3) then
          if (catG==Ai) then
            if (ANY(catGG/=0)) then
              PrXY(x,y,z,2) = PrXY(x,y,z,2) * PrE(y)
            else
              do e=1,3
                PrXY(x,e,:,2) = PrXY(x,e,:,2) * PrE(e)
              enddo
            endif
          else if (ANY(catGG/=0)) then
            PrXY(x,y,z,2) = PrXY(x,y,z,2) * SUM(PrE)                                                                                                    
          else 
            PrXY(x,:,:,2) = PrXY(x,:,:,2) * SUM(PrE)
          endif
        endif
      enddo  ! r
    enddo  ! j
  enddo  ! x
  enddo  ! y (catGG>0 only)
  enddo  ! z (catGG>0 only)
  do m=1,2
    PrL(l,m) = LOG10(SUM(PrXY(:,:,:,m)))! - LOG10(SUM(PrXY(:,:,:,1)))
  enddo
enddo
LL = SUM(PrL(:,2)) - SUM(PrL(:,1))
               
end subroutine MergeSibs

! #####################################################################
subroutine AddFS(A, SB, kB, SA, kA, LL, TopSib, dLL)  ! A/SA FS with any B?
use Global
implicit none

integer, intent(IN) :: A, SB, kB, SA, kA
integer, intent(OUT) :: TopSib    ! most likely FS of A within SB
double precision, intent(OUT) :: LL, dLL(maxSibSize)  ! dLL
integer :: l, x, y, Par(nS(SB,kB)), i, Bj, Ei, f, g,MaybeFS(nS(SB,kB)),&
  z, PA, AncA(2,mxA), AncB(2,mxA), AncPF(2,mxA), h, Inbr(nS(SB,kB)), InbrX      
double precision :: PrL(nSnp, nS(SB,kB),2), PrY(3,2), PrX(3,2), PrZ(3),&
  LLtmp(2), LLUX, PrW(3), ALR

PrL = 0D0
LL = 999D0
dLL = 999D0
TopSib = 0        

Par = 0  ! shared parent 3-kB  (cand. parent(kB) == SB)
MaybeFS = 1

if (nS(SB,kB)==0) then
  LL = 777
  return   ! nobody to be FS with
endif

 call GetAncest(-SB, kB , AncB)
PA = 0
if (A /= 0) then
  PA = Parent(A, 3-kB)
  call GetAncest(A, kA, AncA)
else if (SA /= 0) then   ! TODO: does it matter if kA=kB?
  PA = GpID(3-kB, SA, kA)
   call GetAncest(-SA, kA, AncA)
endif

if (A/=0) then
  if (Parent(A,kB)/=0 .and. Parent(A,kB)/=-SB) then
    LL = 777
  else if (ANY(AncB == A)) then 
    LL = 777
  else
    call CalcAgeLR(A, sex(A), -SB, kB, kB, 1, .TRUE., ALR)
    if (ALR == 777)  LL = 777
  endif
  if (Parent(A,3-kB)==GpID(3-kB,SB,kB) .and. Parent(A,3-kB)/=0) then
    LL = 777
  endif 
else if (SA/=0) then
if (GpID(kB, SA, kA)/=0 .and. GpID(kB, SA, kA)/=-SB) then
    LL = 777
  else if (ANY(AncB(kA, 2:mxA) == -SA)) then
    LL = 777
  else
    call CalcAgeLR(-SA, kA, -SB, kB, kB, 1, .TRUE., ALR)
    if (ALR == 777)  LL = 777
  endif
endif
if (LL /= 999D0) return

InbrX = 0
if (ANY(AncA(kB, 3:mxA) == -SB)) then
  if (A>0 .and. AncA(kB,5-kB)==-SB) then  ! P-O mating
    if (Parent(A,3-kB)<0) then
      InbrX = -1  !
    else
      InbrX = 0  ! or: find out which sibID (matters when low CR)
    endif
  else if (A>0) then
    if (Parent(A,3-kB)/=0 .and. &
      ANY(Parent(SibID(1:nS(SB,kB),SB,kB),3-kB)==Parent(A,3-kB))) then
        InbrX = -2  
    else
      LL = 444
    endif
  else 
    LL = 444   ! TODO: check       
  endif
endif
if (LL /= 999D0) return

if (SA/=0 .and. GpID(3-kA, SB, kB)/=0) then
  do i=1, ns(SA,kA)
    if (Parent(SibID(i,SA,kA), 3-kA) == GpID(3-kA, SB, kB)) then
      LL = 444
      return
    endif
  enddo
endif

Inbr = 0       
do f=1, nS(SB,kB)
  if (NFS(SibID(f, SB, kB))==0) then
    MaybeFS(f) = -1
    cycle
  endif   
  do i=1,nFS(SibID(f, SB, kB))
    Bj = FSID(i, SibID(f, SB, kB))
    if (A == Bj) then
      LL = 888
    else if (A >0) then
      if (Parent(A,3-kB) == Bj) then  
         MaybeFS(f) = 0     ! can't be FS with own parent    
      else if (Parent(Bj, 3-kB) == A) then
        MaybeFS(f) = 0
      else if (Parent(Bj, 3-kB)/=0) then
        call CalcAgeLR(A, Sex(A), Parent(Bj, 3-kB), 3-kB, 3-kB, 1, .TRUE., ALR)
        if (ALR==777)  MaybeFS(f) = 0
      endif
      call CalcAgeLR(A, Sex(A), Bj, Sex(Bj), kB, 2, .TRUE., ALR)
      if (ALR==777)  MaybeFS(f) = 0
      
    else if (SA/=0) then
      if (kA/=kB .and. Parent(Bj, 3-kB) == -SA) then
        MaybeFS(f) = 0  ! cannot be FS with own parent
        LL = 444   ! TODO: implement. 
        cycle
      else if (Parent(Bj, 3-kB)/=0) then
        call CalcAgeLR(-SA, kA, Parent(Bj, 3-kB), 3-kB, 3-kB, 1, .TRUE., ALR)
        if (ALR==777)  MaybeFS(f) = 0
      endif
      call CalcAgeLR(-SA, kA, Bj, Sex(Bj), kB, 2, .TRUE., ALR)
      if (ALR==777)  MaybeFS(f) = 0
    endif
    if (Bj == PA .or. (A/=0 .and. A == Parent(Bj, 3-kB))) then
      MaybeFS(f) = 0
      cycle
    endif
    if (PA>0) then
      if (Parent(PA,1)==Bj .or. Parent(PA,2)==Bj) then
        MaybeFS(f) = 0
        cycle
      endif
    endif
    
    Par(f) = Parent(Bj, 3-kB)
    if (PA/=0 .and. PA/=Par(f) .and. Par(f)/=0) then
      MaybeFS(f) = 0
    else if (Par(f)==0) then
      Par(f) = PA
    endif
  enddo
  if(SA/=0 .and. kA==kB .and. Par(f)/=0) then
  else if (Par(f)<0 .and. A>0) then
    if (GpID(kB, -Par(f), 3-kB) == -SB) then
      Inbr(f) = 1
    endif
  endif
enddo

if (LL /= 999D0) return
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
  LL = 777
  return
endif

do f=1, nS(SB,kB)
  if (nFS(SibID(f, SB, kB))==0) cycle
  if (MaybeFS(f)<1 .or. Par(f)==0 .or. Par(f)==PA) cycle
  call getAncest(Par(f), 3-kB, AncPF)
  if (Par(f)>0) then
    if (A/=0 .and. ANY(AncPF(:,2:mxA)==A)) MaybeFS(f) = 0
    if (SA/=0 .and. ANY(AncPF(kA,2:mxA)==-SA)) MaybeFS(f) = 0
  else if (Par(f) < 0) then
    if (A/=0 .and. ANY(AncPF(:,3:mxA)==A)) MaybeFS(f) = 0
    if (SA/=0 .and. ANY(AncPF(kA,3:mxA)==-SA)) MaybeFS(f) = 0
  endif
enddo
if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
  LL = 777
  return
endif     

if (A/=0 .and. nAgeClasses>1) then  
  do f=1, nS(SB,kB)
    if (MaybeFS(f)<1 .or. Par(f)/=0 .or. &
     Parent(SibID(f, SB, kB), 3-kB)/=0) cycle  
    if (AgeDiff(SibID(f, SB, kB), A) <= 0) cycle
    call CalcU(-SB, kB, A, kB, LLtmp(1))
    call NewPar(SibID(f, SB, kB), 3-kB, A)
    call CalcU(-SB, kB, A, kB, LLtmp(2))
    call NewPar(SibID(f, SB, kB), 3-kB, 0)
    if (LLtmp(1) - LLtmp(2) < TA) then
      MaybeFS(f) = 0
    endif
  enddo
endif

if (ALL(MaybeFS==0 .or. MaybeFS==-1)) then
  LL = 777
  return
endif

dLL = 999D0
if (InbrX > -2) then  ! not: PO or GP-GO mating                                               
do l=1,nSnp
  do f=1, nS(SB,kB)
    if (MaybeFS(f) < 1) cycle
    do x=1,3
      PrX(x,:) = XPr(2,x,l, SB, kB)
      do g=1,nS(SB,kB)
        Bj = SibID(g, SB, kB)
        if (NFS(Bj) == 0) cycle
        if (g==f) then
          if (Inbr(g)==1) then
            call ParProb(l, GpID(3-kB,-Par(g),3-kB), 3-kB, 0,0,PrW)
            do y=1,3
              PrY(y,1) = SUM(AKA2P(y,:,x) * PrW)
            enddo
          else
            call ParProb(l, Par(g), 3-kB, -1,0, PrY(:,1))
          endif
        else
          if (Inbr(g)==1 .and. Parent(Bj, 3-kB) == Par(g)) then
            call ParProb(l, GpID(3-kB,-Par(g),3-kB), 3-kB, 0,0,PrW)
            call ParProb(l, Par(g), 3-kB, Bj,-5, PrY(:,1))  ! no GPs & no Bj & no FS of Bj
            do y=1,3
              PrY(y,1) = PrY(y,1) * SUM(AKA2P(y,:,x) * PrW)
            enddo
          else
            call ParProb(l, Parent(Bj, 3-kB), 3-kB, Bj,-1, PrY(:,1))
          endif
        endif            
        PrY(:,2) = PrY(:,1)  ! 1: FS, 2: HS via 3-k, 3: U 
        do y=1,3 
          do i=1,nFS(Bj)
            if (Genos(l,FSID(i, Bj))==-9) cycle
            PrY(y,:) = PrY(y,:) * OKA2P(Genos(l,FSID(i,Bj)), x, y)
          enddo
          if (g==f) then
            if (Par(g) < 0) then 
              do h = 1, nS(-Par(g), 3-kB)
                Ei = SibID(h, -Par(g), 3-kB)                               
                if (Parent(Ei, kB) == -SB) cycle  
                if (NFS(Ei) == 0) cycle  
                call ParProb(l, Parent(Ei,kB), kB, Ei,-1, PrZ)
                do z=1,3                 
                  do i=1, nFS(Ei)
                    if (FSID(i,Ei) == A) cycle
                    if (Genos(l, FSID(i,Ei))/=-9) then
                      PrZ(z) = PrZ(z) * OKA2P(Genos(l,FSID(i,Ei)),y,z)
                    endif
                  enddo
                enddo
                PrY(y,:) = PrY(y,:) * SUM(PrZ)
              enddo
            endif

            if (A/=0) then
              if (Genos(l, A)/=-9) then
                PrY(y,1) = PrY(y,1) * OKA2P(Genos(l,A), x, y)
                if (PA/=0) then
                  PrY(y,2) = PrY(y,2) * OKAP(Genos(l,A), y, l)
                else
                  PrY(y,2) = PrY(y,2) * OHWE(Genos(l,A), l)
                endif
              endif
            else if (SA/=0) then
              PrY(y,1) = PrY(y,1) * SUM(XPr(1,:,l, SA,kA) *AKA2P(:,x,y))
              if (PA/=0) then
                PrY(y,2) =PrY(y,2) *SUM(XPr(1,:,l, SA,kA) *AKAP(:,y, l))
              else
                PrY(y,2) = PrY(y,2) * SUM(XPr(1,:,l, SA,kA) * AHWE(:,l))
              endif
            endif 
          endif
        enddo  ! y
        do i=1,2
          PrX(x,i) = PrX(x,i) * SUM(PrY(:,i))
        enddo
      enddo  ! g
    enddo  ! x
    PrL(l,f,:) = LOG10(SUM(PrX, DIM=1))
  enddo  ! f   
enddo

  dLL = 777
  do f = 1, nS(SB, kB)
    if (MaybeFS(f)<1) cycle
    Bj = SibID(f, SB, kB)
    if (NFS(Bj) == 0) then
      cycle
    else if (nFS(Bj) == 1) then
      dLL(f) = SUM(PrL(:, f,1)) - SUM(PrL(:,f,2))
    else
      do g=1, nS(SB,kB)
        if (Parent(SibID(g, SB, kB), 3-kB) == Parent(Bj, 3-kB)) then
          dLL(g) = SUM(PrL(:, f,1)) - SUM(PrL(:,f,2))
        endif
      enddo
    endif
  enddo
  
if (A/=0) then
  call CalcU(A,kA, -SB, kB, LLUX)
  LL = MAXVAL(dLL, MASK=dLL/=777) + LLUX   
  TopSib = MAXLOC(dLL, MASK=dLL/=777, DIM=1)                                          
  if(TopSib>0)  TopSib = SibID(TopSib, SB, kB)                                                                 
else if (SA/=0) then
  call CalcU(-SA, kA, -SB, kB, LLUX)
  do f = 1, nS(SB, kB)
    if (dLL(f)==777) cycle
    if (Par(f)==0 .and. nS(SA,kA)>1) then
      dLL(f) = dLL(f) + LLUX
    else  ! consider changes in SA (e.g. inbreeding loops) 
      GpID(3-kB, SA, kA) = Par(f)  
      call CalcCLL(SA,kA)
      call PairUA(-SA, -SB, kA, kB, dLL(f))  
    endif
  enddo
  GpID(3-kB, SA, kA) = PA
  call CalcCLL(SA,kA)
  LL = MaxLL(dLL) 
endif

 else if (InbrX==-2 .and. A>0 .and. SB/=0) then  ! A already HS via 3-k
  call DoAdd(A, SB, kB)
  call CalcU(-SB, kB, A, kB, LL)
  call RemoveSib(A, SB, kB)
else
  LL = 444
endif                                                                     
end subroutine AddFS

! #####################################################################

subroutine AddParent(A, SB, k, LL)  ! is A parent of sibship SB?
use Global
implicit none

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l,x,y,m, G(2), Inbr, n,r, Bi, AncA(2,mxA)
double precision :: PrL(nSnp), PrXY(3,3), PrG(3, 2), PrP(3)

LL = 999D0
Inbr = 0                                                          
AncA = 0        
if (ANY(AgeDiff(SibID(1:nS(SB,k),SB,k), A)<=0)) then 
  LL = 777
  return
else  ! if age A unknown, or all age B's unknown
  call GetAncest(A, k, AncA)
  if (ANY(AncA(k, :) == -SB)) then
    LL = 777
    return
  endif
endif  

G = 0     
do m=1,2
  if (Parent(A,m)/= 0) then  
    if (GpID(m,SB,k)/= 0 .and. GpID(m,SB,k) /= Parent(A,m)) then
      LL = 777
      return
    else
      G(m) = Parent(A,m)
    endif
  else if(GpID(m,SB,k)/=0) then
    G(m) = GpID(m,SB,k)
  endif
enddo

if (G(3-k)/=0) then
  do n=1, nS(SB, k)
    if (nFS(SibID(n,SB,k))==0) cycle
    if (Parent(SibID(n,SB,k), 3-k)==G(3-k)) then
      Inbr = n
    endif
  enddo 
endif

PrL = 0D0
do l=1,nSnp
  if (Inbr==0) then                      
    do m=1,2
      call ParProb(l, G(m), m, A,0, PrG(:,m))    
    enddo
    do x=1,3
      do y=1,3
        PrXY(x,y) = XPr(1,x,l, SB,k) * &
          SUM(AKA2P(x, y, :) *  PrG(y,3-k) * PrG(:, k))
      enddo
    enddo
  else
    call ParProb(l, G(k), k, A,0, PrG(:,k))
    call ParProb(l, G(3-k), 3-k, -1,0, PrG(:,3-k))
    do x=1,3
      do y=1,3
        PrXY(x,y) = SUM(AKA2P(x, y, :) *  PrG(y,3-k) * PrG(:, k))
        do n=1, nS(SB,k)
          Bi = SibID(n, SB, k)
          if (nFS(Bi)==0) cycle
          if (Inbr == n) then
            do r=1, nFS(Bi)
              if (Genos(l, FSID(r, Bi))/=-9) then
                PrXY(x,y) = PrXY(x,y) * OKA2P(Genos(l,FSID(r,Bi)),x,y)
              endif
            enddo
          else
            call ParProb(l, Parent(Bi,3-k), 3-k, Bi,-1, PrP)
            do r=1, nFS(Bi)
              if (Genos(l, FSID(r, Bi))/=-9) then
                PrP = PrP * OKA2P(Genos(l, FSID(r, Bi)), x, :)
              endif
            enddo
            PrXY(x,y) = PrXY(x,y) * SUM(PrP)
          endif
        enddo
      enddo
    enddo
  endif
  do x=1,3
    if (Genos(l,A)/=-9) then
      PrXY(x,:) = PrXY(x,:) * OcA(Genos(l,A), x)
    endif
  enddo
  PrL(l) = LOG10(SUM(PrXY))
enddo
LL = SUM(PrL)

end subroutine AddParent

! #####################################################################

subroutine AddGP(A, SB, k, LL)  ! add A as a grandparent to sibship SB
! TODO: check if sharing a parent 3-k
use Global
implicit none 

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x,y, m, i, cat, catG, curGP, z, g, Bi, r, OpPar, w, Bw
double precision :: PrL(nSnp), PrY(3), PrXY(3,3), LLtmp(2), PrG(3), &
  PrPA(3,2), PrXYZ(3,3,3,2), PrTmp(3), PrP(3), PrLU(nSnp), LLU(2), &
  PrW(3), PrXYW(3,3,3)
logical :: AllFS
  
LL = 999D0
if (Sex(A)/=3) then
  m = Sex(A)
!    if (GpID(m,SB,k) /= 0) then  ! allow replacement 
!        LL = 777
!    endif
else if (GpID(1,SB,k)==0) then
  m = 1
else if (GpID(2,SB,k)==0) then
  m = 2
else
  LL = 777
endif
if (LL==777) return

 cat = 0
 catG = 0      
if (Parent(A, 3-k)==GpID(3-k,SB,k) .and. Parent(A, 3-k) /= 0) then                                        
  cat = 1
else
  do i=1,nS(SB,k)         
    if (nFS(SibID(i,SB,k))==0) cycle
    if (Parent(SibID(i, SB, k), 3-k) == 0) cycle
    if (Parent(SibID(i, SB, k), 3-k) == Parent(A, 3-k)) then
      catG = i
    else if (Parent(SibID(i, SB, k), 3-k) == GpID(3-k,SB,k)) then
      cat = 1
      exit
    else if (Parent(SibID(i,SB,k),3-k) < 0) then
      if (GpID(k, -Parent(SibID(i,SB,k),3-k), 3-k) == -SB) then
        cat = 2
        exit
      endif
    endif
  enddo
endif                                                              

do i=1, nS(SB, k)
  if (AgeDiff(A, SibID(i,SB,k))==999) cycle
  if (AgeDiff(SibID(i,SB,k), A)<=0) then  ! A too young 
    LL = 777
  else if (k==1 .and. m==1 .and. &
   AgePriorM(AgeDiff(SibID(i,SB,k), A)+1, 3)==0.0) then
    LL = 777
  else if (k==2 .and. m==2 .and. &
   AgePriorM(AgeDiff(SibID(i,SB,k), A)+1, 4)==0.0) then
    LL = 777 
  else if (k/=m .and. AgePriorM(AgeDiff(SibID(i,SB,k),A)+1,5)==0.0) then
    LL = 777         
  endif
enddo
if (LL==777) return

AllFS = .FALSE.
Bw = 0
call getFSpar(SB, k, .TRUE., OpPar)
if (OpPar /= 0) then
  AllFS = .TRUE.
  do i=1, nS(SB,k)
    if (nFS(SibID(i,SB,k)) == 0) cycle
    if (Parent(SibID(i,SB,k), 3-k) == OpPar) then
      Bw = SibID(i,SB,k)
    endif
  enddo
endif  

PrL = 0D0
PrLU = 0D0
LLU = 999D0
LLtmp = 999D0             
if (cat==0 .and. catG==0) then  !  .and. .not. AllFS
  do l=1,nSnp
    call ParProb(l, A, 0, 0, 0, PrY)
    call ParProb(l, GpID(3-m, SB, k), 3-m, 0, 0, PrG)
    if (AllFS)  call ParProb(l, OpPar, 3-k, BW, -1, PrW)
    do x=1,3  ! SB
      do y=1,3  ! A
        if (.not. AllFS) then
          PrXY(x,y) = XPr(1,x,l, SB,k) * SUM(AKA2P(x, :, y) * PrG *PrY(y))
        else
          do w=1,3  !OpPar of SB
            PrXYW(x,y,w) = SUM(AKA2P(x, :, y) * PrG *PrY(y)) * PrW(w)
            do i=1, nS(SB,k)
              Bi = SibID(i, SB, k)
              if (Genos(l, Bi) == -9) cycle
              if (Parent(Bi, 3-k) == OpPar) then
                PrXYW(x,y,w) = PrXYW(x,y,w) * OKA2P(Genos(l,Bi),x,w)
              else if (Parent(Bi, 3-k) == 0) then
                PrXYW(x,y,w) = PrXYW(x,y,w) * OKAP(Genos(l,Bi),x,l)
              endif
            enddo
          enddo
        endif
      enddo
    enddo
    if (.not. AllFS) then
      PrL(l) = LOG10(SUM(PrXY))
    else
      PrL(l) = LOG10(SUM(PrXYW))
    endif
  enddo
  LL = SUM(PrL) + Lind(A)
  
else if (cat==0 .and. catG/=0) then 
  do l=1,nSnp
    call ParProb(l, GpID(3-m, SB, k), 3-m, 0, 0, PrG)
    call ParProb(l, Parent(A,3-k), 3-k, -1, 0, PrPA(:,3-k))
    call ParProb(l, Parent(A,k), k, A, 0, PrPA(:,k))
    do x=1,3
      do y=1,3
        do z=1,3
          do g=1,3
            PrTmp(g) = AKA2P(x,g,y) * PrG(g) * SUM(AKA2P(y,z,:) * &
              PrPA(z,3-k) * PrPA(:,k))
          enddo
          PrXYZ(x,y,z,:) = SUM(PrTmp)
          do i=1, nS(SB,k)
            Bi = SibID(i, SB, k)
            if (nFS(Bi)==0) cycle
            if (catG == i) then
              do r=1, nFS(Bi)
                if (Genos(l, FSID(r, Bi))/=-9) then
                  PrXYZ(x,y,z,:) = PrXYZ(x,y,z,:) * OKA2P(Genos(l,FSID(r,Bi)),x,z)
                endif
              enddo
            else
              call ParProb(l, Parent(Bi,3-k), 3-k, Bi,-1, PrP)
              do r=1, nFS(Bi)
                if (Genos(l, FSID(r, Bi))/=-9) then
                  PrP = PrP * OKA2P(Genos(l, FSID(r, Bi)), x, :)
                endif
              enddo
              PrXYZ(x,y,z,:) = PrXYZ(x,y,z,:) * SUM(PrP)
            endif
          enddo 
        enddo
        if (Genos(l,A)/=-9) then
        PrXYZ(x,y,:,1) = PrXYZ(x,y,:,1) * OcA(Genos(l,A), y)
      endif
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ(:,:,:,1)))
    PrLU(l) = LOG10(SUM(PrXYZ(:,:,:,2)))
  enddo
  LL = SUM(PrL)
  LLU(1) = SUM(PrLU) + Lind(A)
  call CalcU(A,k, -SB, k, LLU(2))
  LL = LL - MaxLL(LLU) + LLU(2)                         

else if (cat/=0) then  ! inbreeding loop present / will be created   .or. AllFS
  LLtmp = 999D0
  if (GpID(3-m, SB, k) < 0) then
    call CalcU(-SB, k, GpID(3-m, SB, k), 3-m, LLtmp(1))
  endif
  curGP = GPID(m, SB, k)
  GpID(m, SB, k) = A
  call CalcCLL(SB,k)
  if (curGP<0) then            
    call CalcCLL(-curGP,m)
  endif
  call CalcU(-SB, k, A, 3-k, LL)
  if (GpID(3-m, SB, k) < 0) then
    call CalcU(-SB, k, GpID(3-m, SB, k), 3-m, LLtmp(2))
    LL = LL + (LLtmp(2) - LLtmp(1))
  endif
  GPID(m,SB,k) = CurGP
  call CalcCLL(SB,k)
  if (curGP<0) then          
    call CalcCLL(-curGP,m)
  endif
  call CalcCLL(SB, k)
endif
    
end subroutine AddGP

! #####################################################################
subroutine AddGGP(A, SB, k, LL)
use Global
implicit none
! A a GGP of sibship SB? (only calculating over non-gp-assigned)

integer, intent(IN) :: A, SB, k
double precision, intent(OUT) :: LL
integer :: l, x, m, y,z, AncG(2,mxA), i, v, catG, n,GG
double precision :: PrL(nSnp), PrXYZ(3,3,3), PrZ(3),PrA(3),PrP(3),PrV(3)

LL = 999D0
AncG = 0
if (GpID(1, SB,k)/=0) then
  if (GpID(2, SB,k)/=0) then  ! should be assigned as parent-of-gp
    LL = 777   !(or 888)
  else
    m = 2
  endif
else
  m = 1  ! doesn't really matter (?); GpID(m, SB, k) == 0.
endif
if (LL==777) return

if (Sex(A)<3) then
  n = Sex(A)
else
  n = 1
endif

catG =0
GG = GpID(3-m, SB, k)
if (GG/=0) then
  if (GG==Parent(A,3-m)) then
    catG = 1
  endif
  call GetAncest(GG, 3-m, AncG)   
  if (ANY(AncG == A)) then
    if ((GG>0 .and. AncG(n,2)==A) .or. (GG<0 .and. &
     AncG(n,3-m+2)==A)) then  
      catG = 2  ! already GGP via 3-m; check if double ggp
    else
      LL = 444  ! possible; not yet implemented
    endif
  else
    do v=1,2
      if (Parent(A,v)==0) cycle
      if ((GG<0 .and. ANY(AncG(v, 2:4)==Parent(A,v))) .or. &
       (GG>0 .and. ANY(AncG(v, 1:2)==Parent(A,v)))) then
        LL = 444   ! TODO: stricter implementation?                                                    
      endif
    enddo
  endif
endif
if (GpID(3-k,SB,k) < 0) then
  do i=1,nS(SB,k)
    if (Parent(SibID(i, SB, k), 3-k) == GpID(3-k,SB,k)) then 
      LL = 444
      exit
    endif
  enddo
endif
if (Parent(A,3-k)<0) then
  do i=1,nS(SB,k)
    if (Parent(SibID(i, SB, k), 3-k) == Parent(A,3-k)) then 
      LL = 444
      exit
    endif
  enddo
endif  
if (LL==444) return ! TODO: implement.

if (catG/=0) then  ! age check
  if (GG>0) then
    if (AgeDiff(GG, A) >= 0 .and. catG==1 .and. AgeDiff(GG, A)/=999) &
     LL = 777  ! A older than GG
    if (AgeDiff(GG, A) <= 0 .and. catG==2)  LL = 777  ! GG older than A
  else if (GG<0) then
    do v=1, nS(-GG, 3-m)  ! TODO? AgePriorM; ancestors
      if (AgeDiff(SibID(v,-GG,3-m), A) >= 0 .and. catG==1 .and. &
        AgeDiff(SibID(v,-GG,3-m), A)/=999)  LL = 777
      if (AgeDiff(SibID(v,-GG,3-m), A) <= 0 .and. catG==2)  LL = 777
    enddo
  endif
endif
if (LL == 777) return

PrL = 0D0
do l=1,nSnp
  call ParProb(l, A, 0, 0, 0, PrA)
  if (catG==1) then
    call ParProb(l, GG, 3-m, A, 0, PrZ)
    call ParProb(l, Parent(A,m), m, -1, 0, PrP)
  else if (catG==2) then
    call ParProb(l, GG, 3-m, -4, 0, PrZ)
    if (GG > 0) then  ! TODO: use AncG
      call ParProb(l, Parent(GG,3-n), 3-n, GG, 0, PrP) 
    else if (GG < 0) then
      call ParProb(l, GpID(3-n, -GG,3-m), 3-n, 0, 0, PrP)
    else
      PrP = AHWE(:,l)
    endif
  else
    call ParProb(l, GG, 3-m, 0, 0, PrZ)
  endif
  do x=1,3  ! sibship parent
    do y=1,3
      do z=1,3
        PrXYZ(x,y,z) = XPr(1,x,l, SB,k) * AKA2P(x, y, z) * PrZ(z)  
        if (catG==1) then
          do v=1,3
            PrV(v) = AKAP(y, v, l) * PrA(v) * SUM(AKA2P(v,z,:) * PrP)
          enddo
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(PrV) 
        else if (catG==2) then
          do v=1,3
            PrV(v) = AKAP(y, v, l) * PrA(v) * SUM(AKA2P(z,v,:) * PrP)
          enddo
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(PrV) 
        else
          PrXYZ(x,y,z) = PrXYZ(x,y,z) * SUM(AKAP(y, :, l) * PrA)
        endif
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXYZ))           
enddo
if (catG==1) then
  LL = SUM(PrL)
else if (catG==2 .and. GG>0) then
  LL = SUM(PrL) + Lind(A) - Lind(GG)
else
  LL = SUM(PrL) + Lind(A)
endif

end subroutine AddGGP

! #####################################################################

subroutine addGAU(A, SB, k, m, LL)  ! great-full-avuncular
use Global
implicit none

integer, intent(IN) :: A, SB, k, m
double precision, intent(OUT) :: LL
double precision :: PrL(nSnp), PrGGG(3,2), PrXV(3,3,3,3), PrG(3)
integer :: l, x, y, z, v,g

LL = 999D0
if (GpID(m,SB,k)/=0) then
  LL = 444
else 
  do g=1,2
    if (ANY(SibID(1:ns(SB,k),SB,k)==Parent(A,g))) then
      LL = 777  
    else if (g/=k .and. Parent(A,g)/=0 .and. & 
     ANY(Parent(SibID(1:ns(SB,k),SB,k),g)==Parent(A,g))) then
      LL = 444
    endif
  enddo
endif
if (GpID(3-m,SB,k) /= 0) then
  if(ANY(Parent(SibID(1:ns(SB,k),SB,k), 3-k) == GpID(3-m,SB,k))) then
    LL = 444  ! inbreeding loops, approx. below invalid
  endif
endif    
if (LL /= 999) return

PrL = 0D0
do l=1,nSnp
  do g=1,2
    call ParProb(l, Parent(A,g), g, A, 0, PrGGG(:,g))
  enddo
  call ParProb(l, GpID(3-m,SB,k),3-m,0,0,PrG)
  do x=1,3  ! sibship parent
    do y=1,3
      do z=1,3
        do v=1,3
          PrXV(x,y,z,v) = XPr(1,x,l, SB,k) * SUM(AKA2P(x,y,:) * PrG) * &
           AKA2P(y,z,v) * PrGGG(z,1) * PrGGG(v,2)
          if (Genos(l,A)/=-9) then
            PrXV(x,y,z,v) = PrXV(x,y,z,v) * OKA2P(Genos(l,A), z, v)
          endif
        enddo
      enddo
    enddo
  enddo
  PrL(l) = LOG10(SUM(PrXV))           
enddo
LL = SUM(PrL)

end subroutine addGAU   
                                 
! #####################################################################
subroutine ParentHFS(A, SA, kA, SB, kB, hf, LL)  
! parents of SA and SB HS/FS?
use Global
implicit none

integer, intent(IN) :: A, SA, kA, SB, kB, hf
double precision, intent(OUT) :: LL
integer :: m, G(2), l, x, y, u,v, AncA(2,mxA), AncB(2,mxA), i, j,z, r,&
 Ei, GA, GB,e, DoneA(MaxSibSize), Ai, Bj, nA, AA(maxSibSize),&
 catA(maxSibSize), catB(nS(SB,kB)+1), catG, GGP, OpPar, PA
double precision :: PrG(3,2), PrL(nSnp, 2), PrXV(3,3,3,3,3,2), PrPA(3, 2),&
 LLm(2),PrGA(3), PrGB(3), PrE(3), PrH(3), PrGG(3), ALR
logical :: BallFS

PA = 0
nA = 0
LLm = 999D0       
if (A/=0) then
  PA = Parent(A,kA)
  call GetAncest(A, kA, AncA)
else
  PA = -SA
  call GetAncest(-SA, kA, AncA)
endif
do m=1,2
  if (m/=hf .and. hf/=3) cycle
  if (PA < 0) then
    if (GpID(m,-PA,kA)/=GpID(m,SB,kB) .and. GpID(m,-PA,kA)/=0 .and. GpID(m,SB,kB)/=0) then
      LLm(m) = 777
    endif
    nA = nS(-PA,kA)
    AA(1:nA) = SibID(1:nA, -PA, kA)
  else 
    nA = 1
    AA(1) = A
    if (PA > 0) then
      if (Parent(PA,m)/=GpID(m,SB,kB) .and. parent(PA,m)/=0 .and. GpID(m,SB,kB)/=0) then
        LLm(m) = 777
      endif
    endif
  endif
enddo
if (ALL(LLm == 777)) then
  LL = 777
  return
endif
call GetAncest(-SB, kB, AncB)
 
G = 0
GA = 0
GB = 0       
do m=1,2
  if (m/=hf .and. hf/=3) cycle
  if (AncA(m, kA+2)/=0) then
    if (AncA(m, kA+2) == -SB) then
      LLm(m) = 777
    else if (AncB(m, kB+2)/=0 .and. AncA(m, kA+2)/=AncB(m, kB+2)) then
      LLm(m) = 777
    else if (AncB(m, kB+2)==0) then
      G(m) = AncA(m, kA+2)
    else if (AncB(m, kB+2)==AncA(m, kA+2)) then
      G(m) = AncA(m, kA+2)
      LLm(m) = 888  ! already are sibs
    else
      LLm(m) = 777
    endif
  else
    if (AncB(m,kB+2)/=0 .and. AncB(m,kB+2) == AncA(m, 2)) then
      LLm(m) = 777
    else 
      G(m) = AncB(m, kB+2)
    endif
  endif
  if (hf==3) then  ! FS
    if (ANY(AncA(kB, 3:mxA) == -SB)) then
      LLm = 777
    else if (A>0) then
      if (ANY(AncB == A)) then
        LLm = 777
      endif
    else if (SA/=0) then
      if (ANY(AncB(kA,3:mxA) == -SA)) then
        LLm = 777
      endif
    endif
  endif
    if (A>0) then
    if (hf < 3) then
      call CalcAgeLR(A, Sex(A), -SB, kB, m, 6, .TRUE., ALR)     
    else
      call CalcAgeLR(A, Sex(A), -SB, kB, 0, 5, .TRUE., ALR)    
    endif
  else if (SA/=0) then
    if (hf < 3) then
      call CalcAgeLR(-SA, kA, -SB, kB, m, 3, .TRUE., ALR)
    else
      call CalcAgeLR(-SA, kA, -SB, kB, 0, 2, .TRUE., ALR)
    endif
  endif
  if (ALR == 777)  LLm = 777
enddo
if (ALL(LLm == 777)) then
  LL = 777
  return
endif

if (hf==3) then
  if (LLm(1)==777 .or. LLm(2)==777) then
    LL = 777
  else if (LLm(1)==888 .and. LLm(2)==888) then
    LL = 888  ! already are FS
  endif
else
  if (LLm(hf)==777) then 
    LL = 777
  else if (LLm(hf)==888) then
    LL = 888        
  else if (LLm(3-hf)==888) then
    LL = 777   ! already HS, would become FS
  else 
    GA = AncA(3-hf, kA+2)
    GB = AncB(3-hf, kB+2)
  endif
endif

if (ANY(AncA(kB, 5:mxA)==-SB)) then
  LL = 444  ! highly unlikely (but not strictly impossible: TODO)
else if (AncA(kA,2)/=0 .and. ANY(AncB(kA, 5:mxA) == AncA(kA,2))) then 
  LL = 777
endif
if (LL /=999) return  
   
 catA = 0  
 catB = 0
do i=1, nA
  if (kA/=kB) then
    if (Parent(AA(i), kB) == AncB(kB, 2) .and. AncB(kB, 2)<0) then
      catA(i) = 1
    endif
  else if (kA == kB .and. Parent(AA(i), 3-kA)/=0) then  
    do j=1, nS(SB, kB)
      if (Parent(AA(i), 3-kA) == Parent(SibID(j,SB,kB), 3-kB)) then
        catA(i) = 2
        catB(j) = 2
      endif
    enddo
  endif
  if (Parent(AA(i), 3-kA) /= 0) then
    if (G(3-kA) == Parent(AA(i), 3-kA)) then  ! incl. hf==3
      if (kA==kB) then
        catA(i) = 3  ! (u) 3-kA = 3-kB == hf 
      else if (kA/=kB) then
        catA(i) = 4  ! (z)
      endif 
    else if (kA==hf .and. GA == Parent(AA(i), 3-kA)) then
      catA(i) = 4  ! (z)
    else if (kA==hf .and. GB == Parent(AA(i), 3-kA)) then
      catA(i) = 5  ! (v)
    endif
  endif
enddo    

do i=1, nS(SB, kB)
  if (kA/=kB) then
    if (Parent(SibID(i,SB,kB), kA) ==AncA(kA,2) .and. AncA(kA,2)<0) then
      catB(i) = 1
    endif
  endif
  if (Parent(SibID(i,SB,kB), 3-kB) /= 0) then
    if (G(3-kB) == Parent(SibID(i,SB,kB), 3-kB)) then
      catB(i) = 3  ! (u)  (for hf<3 .and. hf==3)
    else if (kB==hf .and. GA == Parent(SibID(i,SB,kB), 3-kB)) then
      catB(i) = 4  ! (z) (GA of type 3-kB if hf==kB) 
    else if (kB==hf .and. GB == Parent(SibID(i,SB,kB), 3-kB)) then
      catB(i) = 5  ! (v)
    endif
  endif
enddo 

catG = 0
GGP = 0
if (hf<3) then
  if (G(hf) > 0) then
    if (Parent(G(hf), 3-hf) == GA .and. GA/=0) then
      catG = 1
      GGP = Parent(G(hf), hf)
    else if (Parent(G(hf), 3-hf) == GB .and. GB/=0) then
      catG = 2
      GGP = Parent(G(hf), hf)
    endif
  else if (G(hf) < 0) then
    if (GpID(3-hf, -G(hf), hf) == GA .and. GA/=0) then
      catG = 1
      GGP = GpID(hf, -G(hf), hf)
    else if (GpID(3-hf, -G(hf), hf) == GB .and. GB/=0) then
      catG = 2
      GGP = GpID(hf, -G(hf), hf)
    endif
  endif
endif

 BallFS = .FALSE.
 call getFSpar(SB, kB, .TRUE., OpPar)
if (OpPar /= 0) then
  BallFS = .TRUE.
endif
                  
PrL = 0D0
do l=1,nSnp
  do m=1,2
    if (m/=hf .and. hf/=3) cycle
    if (ANY(CatA==3) .or. ANY(CatB==3)) then
      call ParProb(l, G(m), m, -1,0, PrG(:,m)) 
    else if (catG/=0) then
      call ParProb(l, G(m), m, -4, 0, PrG(:,m))
      if (G(m) > 0) then
        call ParProb(l, GGP, hf, G(m), 0, PrGG) 
      else
        call ParProb(l, GGP, hf, 0, 0, PrGG)
      endif
    else
      call ParProb(l, G(m), m, 0,0, PrG(:,m)) 
    endif
  enddo
  if (hf < 3) then
    if (ANY(CatA==4) .or. ANY(CatB==4)) then
      call ParProb(l, GA, 3-hf, -1,0, PrGA)
    else
      call ParProb(l, GA, 3-hf, 0,0, PrGA)
    endif
    if (ANY(CatA==5) .or. ANY(CatB==5)) then
      call ParProb(l, GB, 3-hf, -1,0, PrGB)
    else
      call ParProb(l, GB, 3-hf, 0,0, PrGB)
    endif
  endif
  if (A>0) then
    if (Genos(l,A)==-9) then
      PrL(l,2) = LOG10(SUM(XPr(3,:,l, SB,kB)))
      cycle
    else 
      call ParProb(l, Parent(A,kA), kA, A,-4, PrPA(:,kA))
      call ParProb(l, Parent(A,3-kA), 3-kA, A,0, PrPA(:,3-kA))
    endif
  endif
  
  PrXV = 0
  do x=1,3  ! SA/PA
    do y=1,3  ! SB
      do u=1,3  ! G_hf / G_3-kB (hf==3)
        do z=1,3  ! G_A (hf/=3) / G_kB (hf==3)
          do v=1,3 ! G_B (hf/=3)
            if (hf==3) then  ! 0 for z/=v
              PrXV(x,y,u,z,z,:) = AKA2P(x,u,z) * AKA2P(y,u,z) *&
               PrG(u,3-kB) * PrG(z, kB)
            else
              if (GA < 0 .and. GA == -SB) then
                PrXV(x,y,u,y,v,:) = AKA2P(x,u,y) * AKA2P(y,u,v) *&
                 PrG(u,hf) * PrGB(v)
              else if (GB < 0 .and. GB == -SA) then
                PrXV(x,y,u,z,x,:) = AKA2P(x,u,z) * AKA2P(y,u,x) *&
                 PrG(u,hf) * PrGA(z)
              else if (catG == 1) then
                PrXV(x,y,u,z,v,:) =AKA2P(x,u,z)*AKA2P(y,u,v)*PrG(u,hf)*&
                  SUM(AKA2P(u,z,:) * PrGG) * PrGA(z) * PrGB(v)
              else if (catG == 2) then
                PrXV(x,y,u,z,v,:) =AKA2P(x,u,z)*AKA2P(y,u,v)*PrG(u,hf)*&
                  SUM(AKA2P(u,v,:) * PrGG) * PrGA(z) * PrGB(v)
              else
                PrXV(x,y,u,z,v,:) = AKA2P(x,u,z)*AKA2P(y,u,v)*PrG(u,hf)&
                * PrGA(z) * PrGB(v)
              endif
            endif
            if (A /=0) then
              if (Parent(A, kA)/=0) then
                PrXV(x,y,u,z,v,:) = PrXV(x,y,u,z,v,:) * PrPA(x, kA)
              endif
            endif
          enddo
        enddo
      enddo
      
      DoneA = 0            
      if (ALL(catA==0) .and. ALL(catB==0) .and. nS(SB,kB)>2 .and. .not. BallFS) then
        if (SA/=0) then
          PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * XPr(1,x,l, SA,kA) *&
           XPr(1,y,l, SB,kB)
        else if (A>0) then
         PrXV(x,y,:,:,:,2) =PrXV(x,y,:,:,:,2)*SUM(OKA2P(Genos(l,A),x,:)&
          * PrPA(:,3-kA)) * XPr(1,y,l, SB,kB)
        endif            

      else
        do r=1, nS(SB,kB)
          Bj = SibID(r, SB, kB) 
          if (NFS(Bj) == 0) cycle 
          if (catB(r)==0 .or. catB(r)==2) then
            call ParProb(l, Parent(Bj,3-kB), 3-kB, -1, 0, PrE)
          else
            PrE = 1
          endif                                           

          if (Parent(Bj,3-kB) <0 .and. CatB(r)/=1) then
            do e=1,3
              do v=1, nS(-Parent(Bj,3-kB), 3-kB)
                Ei = SibID(v, -Parent(Bj,3-kB), 3-kB)
                if (nFS(Ei) == 0) cycle
                if (Parent(Ei, kB) == -SB) cycle
                if (A>0 .and. Ei==A) cycle
                if (Parent(Ei, kA) == Parent(AA(1),kA) .and. &
                 Parent(AA(1),kA)/=0) cycle
                if (Ei==Parent(AA(1),kA)) cycle
                call ParProb(l, Parent(Ei, kB), kB, Ei, -1, PrH)
                do i=1, nFS(Ei)
                  if (Genos(l,FSID(i, Ei))==-9) cycle
                  PrH = PrH * OKA2P(Genos(l,FSID(i,Ei)), :, e)
                enddo
                PrE(e) = PrE(e) * SUM(PrH)
              enddo
            enddo
          endif

          if (catB(r)==0 .or. catB(r)==2) then 
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * SUM(PrE)
          else if (catB(r)==1) then  ! Parent(Bj,3-kB) = PA, kA/=kB
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * PrE(x)
          else if (catB(r)==3) then  ! hf==kB, Parent(Bj,3-kB) = GA
            do u=1,3
              PrXV(x,y,u,:,:,1) = PrXV(x,y,u,:,:,1) * PrE(u)
            enddo
          else if (catB(r)==4) then  ! Parent(Bj,3-kB) = G
            do e=1,3
              PrXV(x,y,:,e,:,1) = PrXV(x,y,:,e,:,1) * PrE(e)
            enddo
          else if (catB(r)==5) then  ! hf==kB, Parent(Bj,3-kB) = GB
            do v=1,3
              PrXV(x,y,:,:,v,1) = PrXV(x,y,:,:,v,1) * PrE(v)
            enddo
          endif
     
          do j=1, nFS(Bj) 
            if (Genos(l,FSID(j, Bj))==-9) cycle
            PrE =  PrE * OKA2P(Genos(l,FSID(j,Bj)), y, :)
          enddo

          if (catB(r)==2) then  ! kA==kB, share parent 3-kB
            do v = 1, nA
              Ai = AA(v)
              if (SA/=0 .and. nFS(Ai) == 0) cycle
              if (Parent(Ai, 3-kA)/=Parent(Bj,3-kB)) cycle
              do i=1, MAX(nFS(Ai), 1) 
                if (A/=0 .and. FSID(i, Ai)/=A) cycle
                if (Genos(l,FSID(i, Ai))==-9) cycle
                PrE =  PrE * OKA2P(Genos(l,FSID(i,Ai)), x, :)
              enddo
              doneA(v) = 1
            enddo
          endif
          
          if (catB(r)==0 .or. catB(r)==2) then 
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * SUM(PrE)
          else if (catB(r)==1) then  ! Parent(Bj,3-kB) = PA, kA/=kB
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * PrE(x)
          else if (catB(r)==3) then  ! hf==kB, Parent(Bj,3-kB) = GA
            do u=1,3
              PrXV(x,y,u,:,:,2) = PrXV(x,y,u,:,:,2) * PrE(u)
            enddo
          else if (catB(r)==4) then  ! Parent(Bj,3-kB) = G
            do e=1,3
              PrXV(x,y,:,e,:,2) = PrXV(x,y,:,e,:,2) * PrE(e)
            enddo
          else if (catB(r)==5) then  ! hf==kB, Parent(Bj,3-kB) = GB
            do v=1,3
              PrXV(x,y,:,:,v,2) = PrXV(x,y,:,:,v,2) * PrE(v)
            enddo
          endif
        enddo
        
        do r = 1, nA
          if (doneA(r)==1) cycle
          if (SA/=0 .and. NFS(AA(r)) == 0) cycle
          if (kA/=kB .and. Parent(AA(r),3-kA)==-SB) cycle  ! done
          if (catA(r)==0) then
            call ParProb(l, Parent(AA(r),3-kA), 3-kA, -1, 0, PrE)
          else
            PrE = 1
          endif

          if (Parent(AA(r), 3-kA) < 0 .and. (SA/=0 .or. &
           ANY(FSID(1:nFS(AA(r)), AA(r))==A))) then
            do e=1,3
              do i=1, nS(-Parent(AA(r), 3-kA), 3-kA)
                Ei = SibID(i, -Parent(AA(r), 3-kA), 3-kA)
                if (nFS(Ei) == 0) cycle
                if (Parent(Ei, kB) == -SB) cycle
                if (A>0 .and. Ei==A) cycle
                if (SA/=0 .and. Parent(Ei, kA) == -SA) cycle
                call ParProb(l, Parent(Ei, kA), kA, Ei, -1, PrH)  
                do j=1, nFS(Ei)
                  if (A/=0 .and. FSID(i, Ei)==A) cycle
                  if (Genos(l,FSID(j, Ei))/=-9) then
                    PrH = PrH * OKA2P(Genos(l,FSID(j,Ei)), :, e)
                  endif
                enddo
                PrE(e) = PrE(e) * SUM(PrH)
              enddo
            enddo
          endif
          
          if (catA(r)<3) then
            PrXV(x,y,:,:,:,1) = PrXV(x,y,:,:,:,1) * SUM(PrE)
           else if (catA(r)==3) then 
            do u=1,3
              PrXV(x,y,u,:,:,1) = PrXV(x,y,u,:,:,1) * PrE(u)
            enddo
          else if (catA(r)==4) then 
            do z=1,3
              PrXV(x,y,:,z,:,1) = PrXV(x,y,:,z,:,1) * PrE(z)
            enddo
          else if (catA(r)==5) then 
            do v=1,3
              PrXV(x,y,:,:,v,1) = PrXV(x,y,:,:,v,1) * PrE(v)
            enddo
          endif
          
          do i=1, MAX(1, nFS(AA(r)))
            if (Genos(l,FSID(i, AA(r)))==-9) cycle
            if (SA/=0 .or. FSID(i, AA(r))==A) then 
              PrE =  PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)
              doneA(r) = 2
            else
              PrE =  PrE * OKA2P(Genos(l,FSID(i,AA(r))), x, :)
            endif
          enddo
          
          if (catA(r)<3) then
            PrXV(x,y,:,:,:,2) = PrXV(x,y,:,:,:,2) * SUM(PrE)
           else if (catA(r)==3) then 
            do u=1,3
              PrXV(x,y,u,:,:,2) = PrXV(x,y,u,:,:,2) * PrE(u)
            enddo
          else if (catA(r)==4) then 
            do z=1,3
              PrXV(x,y,:,z,:,2) = PrXV(x,y,:,z,:,2) * PrE(z)
            enddo
          else if (catA(r)==5) then 
            do v=1,3
              PrXV(x,y,:,:,v,2) = PrXV(x,y,:,:,v,2) * PrE(v)
            enddo
          endif
        enddo
      endif
    enddo
  enddo
  PrL(l,1) = LOG10(SUM(PrXV(:,:,:,:,:,1)))
  PrL(l,2) = LOG10(SUM(PrXV(:,:,:,:,:,2)))                                        
enddo

if (ALL(catA==0) .and. ALL(catB==0) .and. ns(SB,kB)>2) then
 LL = SUM(PrL(:,2))
else
  LL = SUM(PrL(:,2)) - SUM(PrL(:,1))
endif

end subroutine ParentHFS

! #####################################################################
subroutine DummyGP(SA, SB, kA, kB, LL)  
! SB GP of SA? via observed or unobserved
use Global
implicit none

integer, intent(IN) :: SA, SB, kA, kB
double precision, intent(OUT) :: LL
integer :: i, m, l, x, y, z, G(2), ggp(2), v, w, Bi, r, AncB(2,mxA),&
  catB, catA, Ai, ix, catG, AncBi(2, mxA)
double precision :: LLGX(2,2), PrL(nSnp), PrZ(3), PrXYZ(3,3,3,3), PrG(3),&
  PrPG(3), PrW(3), LLtmp(maxSibSize, 2), dx(maxSibSize)
logical :: MaybeGP(maxSibSize)                                                       

LL = 999D0
do i=1, nS(SB,kB)
  if (kA /= kB) then
    if (Parent(SibID(i,SB,kB), kA) == -SA) then
      LL = 444
    endif
  else if (kA == kB) then
    do r= 1, nS(SA, kA)
      if (Parent(SibID(i,SB,kB), 3-kB)==Parent(SibID(r,SA,kA), 3-kA) &
       .and. Parent(SibID(i,SB,kB), 3-kB)/=0) then
        LL = 444
      endif
    enddo
  endif
enddo
if (LL == 444)  return  ! TODO

 call GetAncest(-SB, kB, AncB)
if (ANY(AncB(kA,3:mxA) == -SA)) then
  LL = 777 
  return
else if (ANY(AncB(3-kA, 3:mxA) < 0)) then
  do r= 1, nS(SA, kA)
    if (ANY(AncB(3-kA, 3:mxA)/=0 .and. AncB(3-kA, 3:mxA) == &
     Parent(SibID(r,SA,kA), 3-kA))) then
      LL = 444 ! not implemented
      return
    endif
  enddo
endif
G = GpID(:, SA, kA)

MaybeGP = .TRUE.   ! Bi potential GP of SA?
do r=1, ns(SB,kB)
  Bi = SibID(r, sB, kB)
  if (Parent(Bi, 3-kB)==0) cycle
  call getAncest(Bi, kB, AncBi)
  if (ANY(AncBi(kA,:) == -SA)) then
    MaybeGP(r) = .FALSE.
  endif
enddo
   
catA = 0
catB = 0
do r = 1, nS(sB,kB)   ! check if overlap
  Bi = SibID(r, sB, kB)
  if (NFS(Bi) == 0) cycle
  if (Parent(Bi, 3-kB) == G(3-kB) .and. G(3-kB)/=0) then
    catB = Bi
  endif
enddo
do r = 1, nS(sA,kA)   ! check if inbreeding loop
  Ai = SibID(r,SA,kA)
  if (NFS(Ai)==0) cycle
  if (Parent(Ai, 3-kA) == G(3-kA) .and. G(3-kA)/=0) then 
    catA = Ai
  endif
enddo           

catG = 0
do m=1,2
  if (GpID(m, SA, kA) == GpID(m, SB, kB) .and. GpID(m, SA, kA)/=0) then
    catG = m
  endif
enddo
if (catG/=0 .and. catB/=0) then
  LL = 444
  return
endif

LLGX = 999D0
LLtmp = 999D0
GGP = 0
do m=1,2
  if (m==kB .and. GpID(kB, SA, kA) == -SB) then
    LLGX(m,:) = 777
    cycle
  else if (G(m) > 0) then
    if (Parent(G(m), kB) /=0) then
      LLGX(m,:) = 777
    else 
      call AddSib(G(m), SB, kB, LLtmp(1,m))
      call AddFS(G(m), SB, kB,0,m, LLtmp(2,m), ix, dx)
      if (MaxLL(LLtmp(:,m)) < 0) then
        LLGX(m,1) = MaxLL(LLtmp(:,m)) + CLL(SA, kA)
        if (Parent(G(m), kB) /= -SB) then
          LLGX(m,1) = LLGX(m,1) - Lind(G(m))
        endif
      else
        LLGX(m,1) = MaxLL(LLtmp(:,m))
      endif
    endif
    cycle
  else if (G(m) == 0) then
    do r=1, nS(sB,kB)
      Bi = SibID(r, sB, kB) 
      if (Sex(Bi)/=m .and. Sex(Bi)/=3) cycle
      if (.not. MaybeGP(r)) cycle                                      
      call AddGP(Bi, SA, kA, LLtmp(r,m))
      if (LLtmp(r,m) < 0) then
        LLtmp(r,m) = LLtmp(r,m) - Lind(Bi) + CLL(SB,kB)
      endif
    enddo
    LLGX(m,1) = MaxLL(LLtmp(:,m))
  else if (G(m) < 0) then 
    if (GpID(kB, -G(m), m) /=0) then
      LL = 777
    else
      GGP(m) = GpID(3-kB, -G(m), m)
    endif
  endif
  
  PrL = 0D0
  do l=1,nSnp
    if (catB /= 0 .and. m==kB) then
      call ParProb(l, G(3-m), 3-m, catB, -1, PrZ)
    else if (catA/=0) then   ! TODO: catB/=0 .and. catA/=0
      call ParProb(l, G(3-m), 3-m, catA, -1, PrZ)
    else
      call ParProb(l, G(3-m), 3-m, 0, 0, PrZ)
    endif
    if (G(m) /= 0) then
      call ParProb(l, G(m), m, -4, 0, PrG)  ! G(m)'s offspring only,if<0
      PrG = PrG/SUM(PrG)
    else
      PrG = 1
    endif
    call ParProb(l, GGP(m), 3-kB, 0, 0, PrPG)

    PrXYZ = 0
    do x=1,3  ! SA
      do y=1,3  ! parent of SA, offspr of SB
        do v=1,3   ! SB 
          do z=1,3  ! other parent of SA
            PrXYZ(x,y,z,v) = AKA2P(x, y, z) * PrG(y) * PrZ(z) * &
              SUM(AKA2P(y, v, :) * PrPG)
          enddo
          if (catA==0) then
            PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * XPr(1,x,l, SA,kA)   
          else 
            do r=1, nS(sA,kA)
              Ai = SibID(r, sA, kA)  
              if (NFS(Ai) == 0) cycle
              if (Bi == G(m)) cycle
              if (catA==Ai) then
                PrW = 1
              else
                call ParProb(l, Parent(Ai, 3-kA), 3-kA, Ai, -1, PrW)
              endif
              do w=1,3
                do i=1, nFS(Ai)
                  if (Genos(l,FSID(i, Ai))==-9) cycle
                  PrW(w) = PrW(w) * OKA2P(Genos(l,FSID(i,Ai)), x, w)
                enddo
              enddo
              if (catA==Ai) then
                if (m==kA) then
                  PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * PrW
                else if (m/=kA) then
                  PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * PrW(y)
                endif
              else
                PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * SUM(PrW)
              endif
            enddo
          endif
          
          if (catB==0 .or. m/=kB) then
            if (catG==0) then
              PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * XPr(3,v,l, SB,kB) 
            else  ! SA inbred
              call ParProb(l, GpID(3-catG,SB,kB), 3-catG, 0, 0, PrW)
              do z=1,3
                PrXYZ(x,y,z,v) = PrXYZ(x,y,z,v) * SUM(AKA2P(v,z,:) *PrW) *& 
                 XPr(1,v,l, SB,kB)
              enddo
            endif  
          else if (catB/=0) then
            PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * (XPr(2,v,l, SB,kB)/SUM(XPr(2,:,l, SB,kB)))   ! GPs
            do r=1, nS(sB,kB)
              Bi = SibID(r, sB, kB)  
              if (NFS(Bi) == 0) cycle
              if (catB==Bi) then
                PrW = 1
              else
                call ParProb(l, Parent(Bi, 3-kB), 3-kB, Bi, -1, PrW)
              endif
              do i=1, nFS(Bi)
                if (Genos(l,FSID(i, Bi))==-9) cycle
                PrW = PrW * OKA2P(Genos(l,FSID(i,Bi)), v, :)
              enddo
              if (catB==Bi) then
                PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * PrW
              else
                PrXYZ(x,y,:,v) = PrXYZ(x,y,:,v) * SUM(PrW)
              endif
            enddo
          endif  
        enddo
      enddo
    enddo
    PrL(l) = LOG10(SUM(PrXYZ))         
  enddo
  LLGX(m,2) = SUM(PrL)
enddo
LL = MaxLL((/LLGX(:,1), LLGX(:,2)/))

end subroutine DummyGP

! #####################################################################
!       Age priors  
! #####################################################################

subroutine EstBY(s, k, A, DumRel, LPBY)     
! age prior for dummy parent, based on offspring & GP BY. 
! updates DumBY(:, s, k) as side effect
use Global
implicit none

integer, intent(IN) :: s, k, A
logical, intent(IN) :: DumRel                             
double precision, intent(OUT) :: LPBY(nAgeClasses+ maxAgePO)
integer :: i, y, m, x,nOff, Offspr(maxSibSize), sxOff(maxSibSize), GG(2)
double precision :: zero, BYS(maxSibSize, nAgeClasses+ maxAgePO), &
 BYGP(2, nAgeClasses+ maxAgePO), tmpX(nAgeClasses+ maxAgePO)
logical :: ParHasBY, OffHasBY

zero = 0.0D0
if (nAgeClasses==1 .and. A==0) then
  LPBY(1) = zero
  LPBY(2) =  LOG10(zero)
  DumBY(:,s,k) = LPBY
  return
else if (A/=0) then
  LPBY = LOG10(1.D0/(nAgeClasses+ maxAgePO))  ! CHECK
  ParHasBY = .FALSE.
  OffHasBY = .FALSE.
  do m=1,2
    if (Parent(A,m) > 0) then  
      if (BY(Parent(A,m)) > 0)   ParHasBY = .TRUE.
    else if (Parent(A,m) < 0 .and. DumRel) then
      ParHasBY = .TRUE.
    endif
  enddo
  call getOff(A,sex(A),DumRel, nOff, Offspr, sxOff)
  if (noff>0) then
    do i=1, nOff
      if (Offspr(i)>0) then
        if (BY(Offspr(i)) >0) then
          OffHasBY = .TRUE.
          exit
        endif
      else if (Offspr(i)<0 .and. DumRel) then
        OffHasBY = .TRUE.
        exit
      endif
    enddo
  endif
  if (.not. ParHasBY .and. .not. OffHasBY)   return
else if (s/=0) then
  call getOff(-s,k, DumRel, nOff, Offspr, sxOff)  
endif

BYS = 0.0D0  ! number of offspring born in year y                                         
GG = 0
if (s/=0) then
  GG = GpID(:,s,k)
else if (A/=0) then
  GG = Parent(A,:)
endif

do i=1, nOff
  if (Offspr(i)>0) then
    if (BY(Offspr(i))>=0) then
      BYS(i, BY(Offspr(i)) +maxAgePO) = 1D0
    endif
  else if (Offspr(i)<0 .and. DumRel) then
    BYS(i, :) = 10**DumBY(:, -Offspr(i), sxOff(i))
  endif
enddo

BYGP = 0.0D0
do m=1,2
  if (GG(m)>0) then
    if (BY(GG(m))>=0) then
      BYGP(m, BY(GG(m)) +maxAgePO) = 1D0
    endif
  else if (GG(m)<0 .and. DumRel) then
    BYGP(m,:) = 10**DumBY(:, -GG(m), m)
  endif  
enddo

LPBY = 0D0  
do y=1,nAgeClasses +maxAgePO
  if (ANY(BYS(:, 1:y) >= 1D0)) then
    LPBY(y) = LOG10(zero) ! some offspring born in/prior to year y
  else if (ANY(BYGP(:, y:(nAgeClasses +maxAgePO)) >= 1D0)) then
    LPBY(y) = LOG10(zero)  ! some grandparents born in/after year y
  else
    do i=1, nOff
      if (ANY(BYS(i,:)>0)) then
        tmpX = 0
        do x=y+1, nAgeClasses +maxAgePO
          if (x-y >= nAgeClasses) cycle
          if (BYS(i,x)>0) then
            tmpX(x) = BYS(i,x) * AgePriorM(x-y+1, 6+k)
          endif
        enddo
        LPBY(y) = LPBY(y) + LOG10(SUM(tmpX))
      endif
    enddo
    do m=1,2
      if (ANY(BYGP(m,:) > 0)) then
        tmpX = 0
        do x=1, y-1
          if (y-x >= nAgeClasses) cycle
          if (BYGP(m,x)>0) then
            tmpX(x) = BYGP(m,x) * AgePriorM(y-x+1, 6+m)
          endif
        enddo
        LPBY(y) = LPBY(y) + LOG10(SUM(tmpX))
      endif
    enddo
  endif
enddo

LPBY = 10**LPBY
LPBY = LPBY / SUM(LPBY)  ! scale to sum to unity 
LPBY = LOG10(LPBY)
if (s/=0) then
  DumBY(:, s, k) = LPBY
endif

end subroutine EstBY

! #####################################################################

subroutine CalcALRmerge(SA, SB, k, ALR)  ! change in ALR when merging
use Global
implicit none

integer, intent(IN) :: SA, SB, k
double precision, intent(OUT) :: ALR
integer :: i, j

ALR = 0.0
do i = 1, nS(SA,k)
  if (BY(SibID(i,SA,k))<0) cycle
  do j=1, nS(SB, k)
    if (AgeDiff( SibID(i,SA,k), SibID(j,SB,k))/=999) then
      ALR = ALR + LOG10(AgePriorM(ABS(AgeDiff( SibID(i,SA,k), &
       SibID(j,SB,k) ))+1, k)) 
    endif
  enddo
enddo

if (ALR < -HUGE(0.0D0)) then
  ALR = 777
endif

end subroutine CalcALRmerge

! #####################################################################

subroutine CalcAgeLR(A, kA, B, kB, m, focal, DumRel, ALR) ! m: mat/pat relatives
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB, m, focal
logical, intent(IN) :: DumRel                             
double precision, intent(OUT) :: ALR
integer :: x, y, z, gcat, i, AB(2), kAB(2) 
double precision :: BYLR(nAgeClasses+ maxAgePO, 2), zero, &
  ALRtmp(nAgeClasses+ maxAgePO, nAgeClasses+ maxAgePO)

gcat = 0
if (focal==4) then
  if (m==1 .and. kB==1) then
    gcat = 3
  else if (m==2 .and. kB==2) then
    gcat = 4
  else
    gcat = 5
  endif
endif

if (A==0 .or. B==0) then
  ALR = 0
  return  
else if (A>0 .and. B>0) then
  if (AgeDiff(A,B)/=999) then
    if (focal == 1 .or. focal==4) then
      if (AgeDiff(A,B) <=0) then
        ALR = 777
        return
      endif
    endif
    if (focal == 1) then
      ALR = LOG10(AgePriorM(AgeDiff(A,B)+1, 6+kB))
    else if (focal == 2) then
      ALR = LOG10(AgePriorM(ABS(AgeDiff(A,B))+1, 9))
    else if (focal == 3) then
      ALR = LOG10(AgePriorM(ABS(AgeDiff(A,B))+1, m))
    else if (focal == 4) then
      ALR = LOG10(AgePriorM(AgeDiff(A,B)+1, gcat))
    else if (focal == 5 .or. focal==6) then  !FA / HA
      ALR = LOG10(AgePriorM(ABS(AgeDiff(A,B))+1, 6))
    endif
    if (ALR < -HUGE(0.0D0) .or. ALR/=ALR)   ALR = 777
    return
  endif
endif
  
! allocate(DumBY(nAgeClasses +maxAgePO, nInd/2, 2)) 
zero = 0.0
BYLR = LOG10(zero)  ! likelihood ratio to be born in year X 
ALRtmp = LOG10(zero)
ALR = zero                    
AB = (/ A, B /)
kAB = (/ kA, kB /)
do i=1,2
  if (AB(i) < 0) then
    call EstBY(-AB(i), kAB(i), 0, DumRel, BYLR(:, i))
  else if (AB(i)>0) then
    if (BY(AB(i)) >= 0) then
      BYLR(BY(AB(i)) + maxAgePO, i) = LOG10(1.0)
    else if (BY(AB(i)) < 0) then  ! unknown BY  - use offspring & parent BY
      call EstBY(0, 0, AB(i), DumRel, BYLR(:, i))
      if (ALL(BYLR(:, i) <= LOG10(1.D0/(nAgeClasses+ maxAgePO)))) then
        ALR = 0  ! no age info available
        return
      endif
    endif
  endif
enddo

if (focal==1 .or. focal==4) then  ! quick check
  do y=2, nAgeClasses + maxAgePO  ! B 
    if (BYLR(y,2) < -HUGE(0.0D0)) cycle
    ! at oldest possible BY of B:
    if (ALL(BYLR((y-1):(nAgeClasses + maxAgePO), 1) < -HUGE(0.0D0))) then
      ALR = 777
      return
    else
      exit
    endif
  enddo
endif

do y=1,nAgeClasses + maxAgePO  ! B
  if (BYLR(y,2) < -HUGE(0.0D0)) cycle
  do x=1, nAgeClasses + maxAgePO  ! A 
    if (BYLR(x,1) < -HUGE(0.0D0)) cycle
    ALRtmp(x,y) = BYLR(x,1) + BYLR(y,2) 
    if (focal==-1) then  ! A==B
      if (x /= y)   ALRtmp(x,y) = LOG10(zero)
    else if (focal==1 .or. focal==4) then
      if (x > y .and. (x-y) < nAgeClasses) then
        if (focal == 1) then  ! B parent of A
          ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(x - y +1, 6+kB))
        else if (focal == 4) then  ! B GP of A
          ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(x - y +1, gcat))
        endif
      else if (x <= y) then
        ALRtmp(x,y) = LOG10(zero)
      endif
    endif
    z = ABS(x - y) + 1
    if (z<1 .or. z > nAgeClasses) cycle
    if (focal == 2) then  ! FS
      ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(z, 9))
    else if (focal == 3) then  !HS
      ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(z, m))
    else if (focal == 5 .or. focal==6) then  !FA / HA
      if (ANY(AgePriorM(2:nAgeClasses, (/1,2,9/)) > 0)) then
        ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(z, 6))
      else  ! discrete generations
        if (x > y .and. (x-y) < nAgeClasses) then
          ALRtmp(x,y) = ALRtmp(x,y) + LOG10(AgePriorM(z, 6))
        else if (x <= y) then
          ALRtmp(x,y) = LOG10(zero)
        endif
      endif    
    endif  
  enddo
enddo
ALR = LOG10(SUM(10**ALRtmp))  ! sum across age differences
if (ALR < -HUGE(0.0D0) .or. ALR/=ALR)   ALR = 777

end subroutine CalcAgeLR

! #####################################################################

! #####################################################################

subroutine BestRel(LLIN, focal, X, dLL)
use Global
implicit none
! return which relationship is most likely, by threshold TA
! assuming order PO,FS,HS,GG,FAU,HAU,U in LL vector

double precision, intent(IN) :: LLIN(7)
integer, intent(IN) :: focal
integer, intent(OUT) :: X
double precision, intent(OUT) :: dLL   ! diff best vs next best
double precision :: LL(7)
integer :: i,j, maybe(6), Y(7)

X=0
dLL = 0
LL = LLIN

if (ALL(LL(1:6) > 0)) then
  X = 8
  return
endif  
if (focal/=2 .and. LL(2) < 0 .and. LL(2)>=LL(3)) then  
  LL(3) = 333
else if (focal==3 .and. LL(3)>LL(2) .and. LL(3)<0) then
  LL(2) = 333 ! want sib vs non-sib
endif

if ((LL(7) - MAXVAL(LL(1:6), MASK=LL(1:6)<0)) > TA) then  
  X = 7  ! unrelated
else
  maybe = 1
  do i=1,6
    if (LL(i)>0) then
      maybe(i) = 0
    else
      do j=1,7
        if (i==j) cycle
        if (LL(j)>0) cycle
        if ((LL(i) - LL(j)) < TA) then
          maybe(i) = 0   ! i has no longer highest LL
        endif
      enddo
    endif
  enddo
  if (SUM(maybe)==0) then
    X = 8  ! unclear
  else if (SUM(maybe)==1) then
    X = MAXLOC(LL(1:6), MASK=maybe==1, DIM=1)
  endif       
endif

Y = (/(i, i=1,7, 1)/)
if (X<8 .and. X>0) then
  dLL = LL(X) - MAXVAL(LL, MASK=(LL<0 .and. Y/=X))
endif  
      
end subroutine BestRel

! #####################################################################
subroutine BestRel2(LLIN, X, dLL)
use Global
implicit none
! as BestRel, but no threshold, and consider all 1st & 2nd degree rel

double precision, intent(IN) :: LLIN(7)
integer, intent(OUT) :: X
double precision, intent(OUT) :: dLL(2)   ! diff best vs next best
double precision :: LL(7)
integer :: i,j, maybe(6), Y(7)

 X = 0
dLL = 0.0D0
LL = LLIN

if (MAXVAL(LL(1:6), MASK=LL(1:6)<0) - LL(7) < TA .or. &
  ALL(LL(1:6) > 0)) then  
  X = 7  ! unrelated 
else
  maybe = 1
  do i=1,6
    if (LL(i)>0) then
      maybe(i) = 0
    else
      do j=1,7
        if (i==j) cycle
        if (LL(j)>0) cycle
        if ((LL(i) - LL(j)) < 0.01) then
          maybe(i) = 0   ! i has no longer highest LL
        endif
      enddo
    endif
  enddo
  if (SUM(maybe)==0) then
    if (ABS(MaxLL(LL(3:5))-MaxLL(LL))<0.01 .and. COUNT(LL(3:5)<0)>1)then
      X = 9  ! any 2nd degree relative
    else
      X = 8  ! unclear
    endif
  else if (SUM(maybe)==1) then
    X = MAXLOC(LL(1:6), MASK=maybe==1, DIM=1)
  else if (SUM(maybe)>1) then
   if (ABS(MaxLL(LL(3:5))-MaxLL(LL))<0.01 .and. COUNT(LL(3:5)<0)>1) then
      X = 9  ! any 2nd degree relative
    else
      X = 8
    endif
  endif       
endif

Y = (/(i, i=1,7, 1)/)
if (X<8) then
  dLL(1) = LL(X) - MAXVAL(LL, MASK=(LL<0 .and. Y/=X))
  dLL(2) = LL(X) - MaxLL(LL(6:7))
else if (X==9) then
  dLL(1) = MaxLL(LL(3:5)) - MaxLL(LL((/1,2,6,7/)))
  dLL(2) = MaxLL(LL(3:5)) - MaxLL(LL(6:7))
endif  
      
end subroutine BestRel2

! #####################################################################

subroutine UpdateAllProbs(n) 
use Global
implicit none

integer, intent(IN) :: n
integer :: i, k, s, x

do x=1,n
  do k=1,2
    do s=1,nC(k)
      call CalcCLL(s, k)
    enddo
  enddo
  do i=1,nInd 
    call CalcLind(i)
  enddo
enddo

end subroutine UpdateAllProbs

! #####################################################################

subroutine CalcLind(i)
use Global
implicit none

integer, intent(IN) :: i
integer :: l, x, y, k, z, FSi(MaxSibSize), nFSi, j, &
  Inbr(2), Anc(2,mxA), PIK
double precision :: PrL(nSnp), Px(3,2), PrXYZ(3,3,3), PrX(3), PrG(3)
logical :: doFS

if (COUNT(Genos(:,i)/=-9) < nSnp/2.0) then
  return
endif

nFSi = 1
doFS = .FALSE.
if (Parent(i,1)<0 .or. Parent(i,2)<0) then
  doFS = .TRUE.
  nFSi = nFS(FSID(maxSibSize+1,i))
  FSi = FSID(1:maxSibSize, FSID(maxSibSize+1,i))
endif

! PO- and GP-mating
Inbr = 0
Anc = 0
if (Parent(i,1)/=0 .and. Parent(i,2)/=0) then
  call getAncest(i,1,Anc)
  do k=1,2
    if (Anc(3-k,5-k) == Parent(i,3-k)) then
      Inbr(k) = 1
    endif
  enddo
endif
do k=1,2
  if (Inbr(k)/=0 .and. Parent(i,3-k)>0) then
    PIK = Parent(i,3-k)
  else
    PIK = 0
  endif
enddo

PrL = 0D0
do l=1,nSnp
  do k=1,2
    call ParProb(l, Parent(i,k), k, i,-1, Px(:,k))
    if (Inbr(k)==1) then
      call ParProb(l, Anc(k,k+2), k, PIK,0, PrG)
    endif
  enddo
  do y=1,3
    do z=1,3
      if (ANY(Inbr==1)) then
        do k=1,2
          if (Inbr(k)==-1) then
            PrXYZ(:,y,z) =AKA2P(:,y,z) * SUM(AKA2P(y,z,:)*PrG) *&
             Px(y,k) * Px(z,3-k)
          endif
        enddo
      else
        PrXYZ(:,y,z) = AKA2P(:, y, z) * Px(y,1) * Px(z,2)
      endif 
      if (doFS) then
        do j=1, nFSi
          if (FSi(j) == i) cycle
          if (Genos(l, FSi(j))/=-9) then
            PrXYZ(:,y,z) = PrXYZ(:,y,z) *OKA2P(Genos(l,FSi(j)),y, z)
          endif
        enddo
      endif
    enddo
  enddo
  
  do x=1,3
    PrX(x) = SUM(PrXYZ(x,:,:))/SUM(PrXYZ)
  enddo
  PrX = PrX / SUM(PrX)
  if (Genos(l,i)/=-9) then
    PrX = OcA(Genos(l,i), :) * PrX     
  endif   
  PrL(l) = LOG10(SUM(PrX))
  LindX(:,l, i) = PrX
  LindG(:, l, i) = PrX / SUM(PrX)  ! used in parprob
enddo
Lind(i) = SUM(PrL)

 if (Lind(i)< -99999 .or. Lind(i)> 0 .or. Lind(i)/=Lind(i)) then
  call Erstop("Invalid individual LL")
endif                                                              
end subroutine CalcLind

! #####################################################################
subroutine CalcCLL(s, k) 
use Global
implicit none
! returns XPr: likelihood;  DumP: probability, scaled  (no age prior.),
! split into 1: sibs only 2: gp effect only, 3: all

integer, intent(IN) :: s, k ! S: sibship number, k: mat(1),pat(2),unk(3)
integer :: l, x, i, Ei, r, y, z, g, Ri, v, cat, catG, e, &
  IsInbr(ns(s,k)), HasInbr(ns(s,k), ns(s,k)), AncR(2,mxA), &
  UseEE(ns(s,k)), Sibs(ns(s,k)), MatePar(ns(s,k)), FSX
double precision :: PrL(nSnp), PrY(3,2), PrYp(3,ns(s,k)), PrGG(3,2),&
 PrZ(3),PrXZ(3,3,2), PrE(3), PrEE(3, ns(s,k))

!if (ns(s,k) == 0) then
!  CLL(s,k) = 0D0   ! during CalcParentLLR
!  return
!endif
if (ALL(GpID(:,s,k)==0) .and. ALL(SibID(1:ns(s,k),s,k)==0)) then
  call Erstop("Empty sibship!")
else if (s > nC(k)) then
  call Erstop("s out of bounds!")
endif

Sibs = SibID(1:ns(s,k), s, k)
 call FindEE(Sibs, ns(s,k), 0, k, UseEE, MatePar)

!================= 
cat = 0
catG = 0
IsInbr = 0
HasInbr = 0
AncR = 0
do r=1,nS(s,k)
  Ri = Sibs(r)  !SibID(r, s, k)
  if (Parent(Ri, 3-k)==0) cycle
  if (Parent(Ri, 3-k)==GpID(3-k,s,k) .and. nFS(Ri)/=0) then  
    cat = Ri
    UseEE(r) = 0
  endif
  do v=1, nS(s,k)
    if (r==v) cycle
    if (nFS(Sibs(v))==0) cycle
    do i=1, nFS(Sibs(v))
      if (Parent(Ri, 3-k) == FSID(i, Sibs(v))) then
        IsInbr(r) = FSID(i, Sibs(v))
        HasInbr(v,i) = r !-1
      endif
    enddo
  enddo
  if (IsInbr(r)/=0) cycle
  call GetAncest(Ri,k,AncR)
  if (AncR(k, 5-k) == -s) then 
    IsInbr(r) = Parent(Ri, 3-k)   ! via dummy
  endif
enddo  

FSX = 0  
if (ALL(GpID(:,s,k)<0) .and. cat==0) then  ! check if sibship par inbred
  if (GPID(1,s,k) == GPID(1, -GPID(2,s,k),2)) then
    catG = 2
  else if (GPID(2,s,k) == GPID(2, -GPID(1,s,k),1)) then
    catG = 1
  else 
    do i=1, ns(-GpID(1,s,k), 1)
      if (Parent(SibID(i, -GpID(1,s,k), 1), 2) == GpID(2,s,k)) then  ! FS of dummy par
        catG = 3
        if (nFS(SibID(i, -GpID(1,s,k), 1))/=0) then
          FSX = SibID(i, -GpID(1,s,k), 1)
        endif
      endif
    enddo 
  endif
endif

PrL = 0D0       
do l=1,nSnp
  do g=1,2   !grandparents
    if (g/=k .and. cat>0) then
      call ParProb(l, GpID(g,s,k), g, -1,0, PrGG(:,g))
    else if (catG==g) then
      PrGG(:,g) = XPr(3,:,l, -GpID(g,s,k),g)
    else if (catG==3) then
      call ParProb(l, GpID(g,s,k), g, FSX,-1, PrGG(:,g))
    else
      call ParProb(l, GpID(g,s,k), g, 0,0, PrGG(:,g))
    endif
  enddo
  
  do x=1,3  ! genotype dummy parent
    do z=1,3
      if (catG==0) then
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k) * PrGG(:,k))  ! GPs
      else if (catG==k) then
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(:,k))
      else if (catG==3-k) then
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k))
      else if (catG==3) then
        PrY(:,1) = PrGG(:,k)
        do y=1,3
          do i=1, nFS(FSX)
            if (Genos(l,FSID(i,FSX))==-9) cycle
            PrY(y,1) = PrY(y,1) * OKA2P(Genos(l,FSID(i,FSX)), y, z)
          enddo
        enddo
        PrXZ(x,z,:) = SUM(AKA2P(x, z, :) * PrGG(z,3-k) * PrY(:,1))
      endif
    enddo
  enddo
  if (catG>0) then
    do r=1,2
      PrXZ(:,:,r) = PrXZ(:,:,r)/SUM(PrXZ(:,:,r)) 
    enddo
  endif
  
  do r=1, nS(s,k)
    Ri = Sibs(r)
    if (NFS(Ri) == 0) cycle
    if (cat==Ri) then
      PrYp(:,r) = 1
    else if (IsInbr(r) == 0) then
      call ParProb(l, Parent(Ri, 3-k), 3-k, -1,0, PrYp(:,r)) 
    endif
  enddo
  
  do z=1,3
    if (z>1 .and. cat==0) cycle
  do x=1,3
    XPr(2,x,l, s,k) = SUM(PrXZ(x,:,2))  ! GP 
    do r=1, nS(s,k)
      Ri = Sibs(r)  ! array with IDs
      if (NFS(Ri) == 0) cycle  ! moved to its FS
      if (IsInbr(r) > 0) then
        cycle
      else if (IsInbr(r) < 0) then
        call ParProb(l, GpID(3-k,-Parent(Ri, 3-k),3-k), 3-k, 0,0, PrZ) 
        do y=1,3
          PrYp(y,r) = SUM(AKA2P(y,x,:) * PrZ)
        enddo
      else if (UseEE(r) /= 0) then
        call ParProb(l, GpID(k,-Parent(Ri, 3-k),3-k), k, 0,0, PrZ)
        do y=1,3
          do e=1,3
            PrE(e) = SUM(AKA2P(y,e,:) * PrEE(e,UseEE(r)) * PrZ)
          enddo
          PrYp(y,r) = SUM(PrE)
        enddo
      endif

      do y=1,3   ! parent 3-k 
        if (cat==Ri .and. y/=z) cycle      
        PrY(y,:) = PrYp(y,r)  ! dim2: 1:all, 2:non-sibs only
        do i=1, nFS(Ri)  ! default: nFS = 1
          if (HasInbr(r,i)==0) then
            if (Genos(l,FSID(i, Ri))/=-9) then
              PrY(y,1) = PrY(y,1) * OKA2P(Genos(l,FSID(i,Ri)), x, y)
            endif
          else
            do e=1,3
              PrE(e) = AKA2P(e, x, y)
              if (Genos(l,FSID(i, Ri))/=-9) then
                PrE(e) = PrE(e) * OcA(Genos(l,FSID(i,Ri)), e)
              endif
              do v=1, nS(s,k)
                if (IsInbr(v)==FSID(i,Ri)) then  
                  if (Genos(l,Sibs(v))==-9) cycle
                  PrE(e) = PrE(e) * OKA2P(Genos(l,Sibs(v)), x, e)
                endif
              enddo
            enddo
            PrY(y,1) = PrY(y,1) * SUM(PrE)
          endif
        enddo                
        
        if (Parent(Ri, 3-k) < 0) then
          do v=1, nS(-Parent(Ri, 3-k), 3-k)
            Ei = SibID(v, -Parent(Ri, 3-k), 3-k)  
            if (NFS(Ei) == 0) cycle
            if (Parent(Ei, k) == -s) cycle
            call ParProb(l, Parent(Ei, k), k, Ei,-1, PrE)        
            do i=1, nFS(Ei) 
              if (Genos(l,FSID(i, Ei))/=-9) then
                PrE = PrE * OKA2P(Genos(l,FSID(i,Ei)), :, y)
              endif
            enddo
            PrY(y,:) = PrY(y,:) * SUM(PrE)  
          enddo
        endif
      enddo 
      
      do v=1,2     ! all; non-sibs only        
        if (cat==Ri) then
          PrXZ(x,z,v) = PrXZ(x,z,v) * PrY(z,v)   
        else if (cat/=0) then
          PrXZ(x,z,v) = PrXZ(x,z,v) * SUM(PrY(:,v))
        else
          PrXZ(x,:,v) = PrXZ(x,:,v) * SUM(PrY(:,v))
        endif
      enddo
      PrEE(:,r) = PrY(:,2)
    enddo ! r 
  enddo ! x
  enddo ! z (cat/=0 only)
  do x=1,3  ! account for GP, dumm offspr & connected sibships
    XPr(3,x,l, s,k) = SUM(PrXZ(x,:,1))/ SUM(PrXZ(:,:,2))
  enddo
  do x=1,3
    DumP(x,l, s,k) = XPr(3,x,l, s,k)/ SUM(XPr(3,:,l, s,k))
    XPr(1,x,l, s,k) = XPr(3,x,l, s,k) / XPr(2,x,l, s,k) 
  enddo 
  PrL(l) = LOG10(SUM(XPr(3,:,l, s,k)))
enddo
 CLL(s,k) = SUM(PrL) 
 
end subroutine CalcCLL

! #####################################################################

subroutine ParProb(l, i, k, A, B, prob)  ! TODO: B<0 : no GP
use Global
implicit none

integer, intent(IN) :: l, i, k, A,B
double precision, intent(OUT) :: prob(3)
integer :: x,j, AB(2), A1
double precision :: PrP(3, 2), PrY(3)

prob = AHWE(:, l)                                   
 if (i == 0) then  ! no parent
   if (A==-4 .or. B==-4 .or. B==-5) then
    prob = 1
  else
    prob = AHWE(:, l)
  endif
else if (i > 0) then  ! real parent
  if ((A==-4 .or. B==-4 .or. B==-5) .and. Lind(i)/=0) then
    if (Genos(l,i)/=-9) then
      prob = OcA(Genos(l,i),:)
    else
      prob = 1
    endif
  else
    prob = LindG(:, l, i)  ! =AHWE if Lind(i) not yet calculated
  endif
else if (i < 0) then  ! dummy parent
   if (ns(-i,k)==0) then  ! during CalcParentLLR
    if (ANY(GpID(:,-i,k)/=0)) then
      prob = XPr(2,:,l, -i, k) / sum(XPr(2,:,l, -i, k))
    else
      prob = AHWE(:,l)
    endif                              
  else if (A==0) then   ! probability
    prob = DumP(:,l, -i,k)    
  else if (A == -1) then  ! grandparent contribution only
    prob = XPr(2,:,l, -i, k) 
  else if (A==-4) then  ! offspring contribution only
    prob = XPr(1,:,l, -i, k)        
  else if (A>0) then   ! exclude indiv A from calc & standardise
    if ((Genos(l,A)==-9 .and. (nFS(A)<=1 .or. B>=0)) .or. &
     Parent(A,k)/=i) then
      prob = DumP(:,l, -i,k)
    else
      AB = (/ A, B /)
      do j=1,2
        if (j==2 .and. B<=0) cycle
        if (Parent(AB(j), 3-k)==0) then  
          PrP(:,j) = AHWE(:,l)
        else if (Parent(AB(j), 3-k)>0) then
          if (Genos(l,Parent(AB(j), 3-k)) /= -9) then
            PrP(:,j) = OcA(Genos(l, Parent(AB(j),3-k)), :)
          else
            PrP(:,j) = LindG(:, l, Parent(AB(j), 3-k))
          endif
        else if (Parent(AB(j), 3-k)<0) then  
          PrP(:,j) = DumP(:,l, -Parent(AB(j),3-k), 3-k) 
        endif
      enddo

      A1 = 0
      if (B < 0) then
        A1 = FSID(maxSibSize+1, A)
      endif               

      do x = 1, 3
        if (B>=0 .or. (B==-1 .and. nFS(A)<=1)) then
          prob(x) = XPr(3,x,l, -i, k)
          if (Genos(l,A)/=-9) then
            prob(x) =  prob(x)/ SUM(OKA2P(Genos(l,A),x,:) *PrP(:,1))    
          endif  
        else if (B==-1) then  ! exclude all FS of A
          PrY = PrP(:,1)
          do j=1, nFS(A1)
            if (Genos(l,FSID(j, A1))==-9) cycle
            PrY = PrY * OKA2P(Genos(l,FSID(j, A1)), x, :)
          enddo
          prob(x) = XPr(3,x,l, -i, k) / SUM(PrY)
        else if (B==-4) then ! exclude both GPs & A
          if (ns(-i,k)==1) then
            prob = 1
          else
            prob(x) = XPr(1,x,l,-i,k) !*AHWE(x,l)  !?
            if (Genos(l,A)/=-9) then
              prob(x) =  prob(x)/SUM(OKA2P(Genos(l,A),x,:)* PrP(:,1))
            endif
          endif     
        else if (B==-5) then ! exclude both GPs & A & FS of A
          if (ns(-i,k)==1) then
            prob = 1
          else
            PrY = PrP(:,1)
            do j=1, nFS(A1)
              if (Genos(l,FSID(j, A1))==-9) cycle
              PrY = PrY * OKA2P(Genos(l,FSID(j, A1)), x, :)
            enddo
            prob(x) = XPr(1,x,l,-i,k) / SUM(PrY)   ! *AHWE(x,l)
          endif
        endif
        if (B>0) then
          if (Genos(l,B)==-9) cycle
          prob(x) = prob(x) / SUM(OKA2P(Genos(l,B), x, :) * PrP(:,2)) 
        endif
      enddo
    endif
  endif
endif
if (.not. ALL(prob==1)) then
  prob = prob/SUM(prob)
endif    

 if (ANY(prob< 0)  .or. ANY(prob/=prob)) then  ! .or. ANY(prob>1)
  call Erstop("Invalid ParProb!")
endif 

end subroutine ParProb

! #####################################################################

subroutine Connected(A, kA, B, kB, Con)
use Global
implicit none

integer, intent(IN) :: A, kA, B, kB
logical, intent(OUT) :: Con
integer :: i, j, m, nA, nB, AA(maxSibsize), BB(maxSibsize), n

 Con = .FALSE.
if (A==0 .or. B==0) then
  Con = .FALSE.
  return
endif

if (A>0) then
  nA = 1
  AA(1) = A
else
  nA = nS(-A,kA)
  AA(1:nA) = SibID(1:nA, -A, kA)
endif

if (B>0) then
  nB = 1
  BB(1) = B
else
  nB = nS(-B,kB)
  BB(1:nB) = SibID(1:nB, -B, kB)
endif

do j=1, nB
  do i=1, nA
    do m=1,2  
      if (Parent(AA(i), m) < 0) then
        if (Parent(AA(i),m) == Parent(BB(j),m)) then
          Con = .TRUE.
          return
        else if(ANY(GpID(:,-Parent(AA(i), m),m) == BB(j))) then
          Con = .FALSE.
!          return
        else if (A<0 .and. m==kA) then
          if(ANY(GpID(:,-Parent(AA(i), m),m) < 0)) then
            do n=1,2
              if (GpID(n,-Parent(AA(i), m),m) == Parent(BB(j),n) .and. &
                Parent(BB(j),n)<0) then
                Con = .TRUE.
                return
              endif 
            enddo
          endif
        endif
      endif
      if (Parent(BB(j),m)<0) then
        if (ANY(GpID(:,-Parent(BB(j),m),m) == AA(i))) then
          Con = .FALSE. 
!          return
        else if (B<0 .and. m==kB) then
          if(ANY(GpID(:,-Parent(BB(j),m),m) < 0)) then
            do n=1,2
              if (GpID(n,-Parent(BB(j), m),m) == Parent(AA(i),n) .and. &
                Parent(AA(i),n)<0) then
                Con = .TRUE.
                return
              endif 
            enddo
          endif
        endif
      endif
    enddo
  enddo
enddo

end subroutine Connected

! #####################################################################

subroutine GetAncest(A, kIN, Anc)
use Global
implicit none

integer, intent(IN) :: A, kIN
integer, intent(OUT) :: Anc(2, mxA)  ! 32 = 5 generations
integer :: m, j, i, k

if (A==0) return
k = kIN
if (kIN<1 .or. kIN>2) then
  if (A>0)  k = 1  
endif

Anc = 0
if (A > 0) then  
  if (Sex(A)/=3) then
    Anc(Sex(A),1) = A
  else
    Anc(1, 1) = A
  endif
  Anc(:, 2) = Parent(A, :)
else if (A < 0) then
  Anc(k, 2) = A
endif
if (ALL(Anc(:, 2)==0)) return
do j = 2, mxA/2  
  do m = 1, 2
    i = 2 * (j-1) + m
    if (Anc(m, j) > 0) then
      Anc(:, i) = Parent(Anc(m, j), :)
    else if (Anc(m, j) < 0) then
      Anc(:, i) = GpID(:, -Anc(m, j), m)  
    endif
  enddo
  if (j==2 .and. ALL(Anc(:, 3:4) == 0))  return
  if (j==4 .and. ALL(Anc(:, 5:8) == 0))  return
  if (j==8 .and. ALL(Anc(:, 9:16) == 0))  return
enddo

if ((A>0 .and. ANY(Anc(:, 2:mxA)==A)) .or. (A<0 .and. &
 ANY(Anc(k,3:mxA)==A))) then
!  call intpr ( "Female ancestors: ", -1, Anc(1,1:8), 8)
!  call intpr ( "Male ancestors: ", -1, Anc(2,1:8), 8)
  call Erstop("An individual is its own ancestor!")
endif

end subroutine GetAncest

! #####################################################################

subroutine CalcParentLLR
! Calc parental LLR (vs next most likely relationship)
use Global
implicit none

integer :: i, k, s, CurPar(2), m, nonG(6), CurGP(2), g, curNFS
double precision :: LLA(7), LLtmp(2,2,2), LLAA(7)
logical :: FSM

do i=1, nInd
  if (Error/=0) return
  if (MODULO(i,200)==0) call rchkusr()
  if (quiet<1 .and. ANY(nC>0) .and. nInd>1500) then
    if (MODULO(i,500)==0)  call intpr(" ", 1, i, 1)
  endif
  if (Parent(i,1)==0 .and. Parent(i,2)==0) cycle
  CurPar = Parent(i,:)
  curNFS = nFS(i)
  do k=1,2  ! remove i from sibgroup
    call NewPar(i, k, 0)
  enddo
  
  LLtmp = 999D0
  do m=1,2  ! m=1: no opp. sex parent;  m=2: with opp. sex parent
    if (m==2 .and. (CurPar(1)==0 .or. CurPar(2)==0)) cycle
    do k=1,2  ! mother, father
      if (m==1 .and. CurPar(k) == 0) cycle
      if (m==2) then  ! temp. assign parent 3-k
        call NewPar(i, 3-k, curPar(3-k))
      endif
      
      if (CurPar(k) > 0) then
        call CheckPair(i, CurPar(k), k, 1, LLA, LLAA)
        LLtmp(1,k,m) = LLA(1)
        LLtmp(2,k,m) = MaxLL(LLA(2:7))
      else if (CurPar(k) < 0) then
        call CheckAdd(i, -CurPar(k), k, .FALSE., LLA, 7)
        if (m==1 .and. Complx>0) LLA(2) = 333   ! FS does not count here
        LLtmp(1,k,m) =  MaxLL(LLA(2:3))
        LLtmp(2,k,m) =  MaxLL((/LLA(1), LLA(4:7)/))
      endif
      
      if (m==2) then  
        call NewPar(i, 3-k, 0)
      endif
    enddo
  enddo
  do k=1,2  ! max with - max w/o 
    if (LLtmp(1,k,1)>0) then
      LR_parent(i,k) = LLtmp(1,k,1)  ! something wrong
    else
      LR_parent(i,k) = LLtmp(1,k,1)-LLtmp(2,k,1)
    endif
  enddo
  if (CurPar(1)/=0 .and. CurPar(2)/=0) then
    if (Complx>0) then
      LR_parent(i,3) = MIN(LLtmp(1,1,2) -MaxLL((/LLtmp(2,1,2), &
        LLtmp(:,1,1)/)), LLtmp(1,2,2) -MaxLL((/LLtmp(2,2,2), &
        LLtmp(:,2,1)/)))
    else
      LR_parent(i,3) = LR_parent(i,1)
    endif
  endif
  
  do k=1,2
    call NewPar(i, k, CurPar(k))    ! restore
  enddo
enddo

!parents of dummies (Sibship GPs)
nonG = (/1,2,3,5,6,7/)
do k = 1,2
  do s=1, nC(k)
    CurGP = GpID(:, s, k)
    GpID(:, s, k) = 0
    call CalcCLL(s,k)        
    LLtmp = 999D0
    do m=1,2
      if (m==2 .and. (CurGP(1)==0 .or. CurGP(2)==0)) cycle
      do g=1,2
        if (m==1 .and. CurGP(g) == 0) cycle
        if (m==2) then  ! temp. assign GP 3-g
          GpID(3-g, s, k) = CurGP(3-g)
          call CalcCLL(s,k)
        endif
        
        if (curGP(g) > 0) then
          call checkAdd(CurGP(g),s,k, .FALSE., LLA, 7)  !B=GP+CurGP(m)_7
        else if (curGP(g) < 0) then
          call checkMerge(s, -CurGP(g), k, g, .FALSE., LLA, 4, FSM)  
          if (m==1) then
            call PairUA(-s, CurGP(g), k, g, LLA(4))  
          endif
        endif
        LLtmp(1,g,m) = LLA(4)
        LLtmp(2,g,m) = MaxLL(LLA(nonG))             
        if (m==1) then
          LR_GP(g,s,k) = LLtmp(1,g,m) - LLtmp(2,g,m)
        else if (m==2) then  ! reset to 0
          GpID(3-g, s, k) = 0
          call CalcCLL(s,k)
        endif 
      enddo
    enddo
    if (CurGP(1)/=0 .and. CurGP(2)/=0) then
      LR_GP(3,s,k) =MINVAL(LLtmp(1,:,2) -MAX(LLtmp(1,:,1),LLtmp(2,:,2)))       
    endif   
    GpID(:,s,k) = CurGP 
    call CalcCLL(s, k)
  enddo
enddo

end subroutine CalcParentLLR

! #####################################################################

subroutine NewPar(A, k, P)
use Global
implicit none

integer, intent(IN) :: A, k, P
integer :: nOffP, OffP(maxSibSize), sxOffP(maxSibSize), i

if (Parent(A,k) == P .and. P==0) return

if (Parent(A, k) < 0) then
  call RemoveSib(A, -Parent(A, k), k)
else if (Parent(A,k) > 0) then
  call RemoveFS(A)
  Parent(A,k) = 0
endif
   
if (P < 0) then
  call DoAdd(A, -P, k)
else if (P > 0) then
  if (Parent(A, 3-k) /= 0) then
    nOffP = 0 
    if (ANY(Parent(:, k) == P)) then
      call getOff(P, k, .FALSE., nOffP, OffP, sxOffP)
    endif
    Parent(A, k) = P
    if (nOffP > 0) then
      do i=1, nOffP
        if (Parent(A,3-k) == Parent(OffP(i), 3-k)) then
          call MakeFS(A, OffP(i))
        endif
      enddo
    endif
    if (Parent(A,3-k) < 0) then
      call CalcCLL(-Parent(A,3-k), 3-k)
    endif
  else
    Parent(A, k) = P  
  endif
endif
 call CalcLind(A)
if (Parent(A,3-k) < 0) then
  call CalcCLL(-Parent(A,3-k), 3-k)
  call CalcLind(A)
endif   
    
end subroutine NewPar

! #####################################################################

subroutine RemoveFS(A)
use Global
implicit none

integer, intent(IN) :: A
integer :: op, np, i, j

op = A   ! needs initialising to avoid compiler warning, not used.
np = op
if (nFS(A) == 1) then
  return
else if (nFS(A) > 1) then
  op = A
  np = MINVAL(FSID(1:nFS(A), A), MASK=(FSID(1:nFS(A),A)/=A)) 
else if (nFS(A) == 0) then
  op = FSID(maxSibSize+1, A)  ! 'primary' sib
  np = op
endif

i = 2  ! 1st one stays op
do j=1, nFS(op)
  if (FSID(j,op)==A) cycle
  if (FSID(j,op)==np) cycle
  FSID(i, np) = FSID(j, op)
  if (op /= np) then
    FSID(maxSibSize+1, FSID(j, op)) = np
  endif
  i = i+1
enddo

nFS(np) = nFS(op)-1
FSID(maxSibSize+1, np) = np
nFS(A) = 1
FSID(:,A) = 0
FSID(1,A) = A
FSID(maxSibSize+1, A) = A

end subroutine RemoveFS

! #####################################################################

subroutine RemoveSib(A, s, k)  ! removes individual A from sibship s
use Global
implicit none

integer, intent(IN) :: A, s, k
integer :: u, j, h, inS

 call RemoveFS(A)

inS = 0
do u=1,nS(s,k)
  if (SibID(u,s,k)==A) then
    inS = 1
    if (u<nS(s,k)) then
      do h=u, nS(s, k)-1  ! drop HS
        SibID(h, s, k) = SibID(h+1, s, k)
      enddo
    endif
    SibID(nS(s,k), s, k) = 0
  endif
enddo
if (inS==1)  nS(s,k) = nS(s,k)-1

Parent(A, k) = 0
if (ns(s,k)>0) then
   call calcCLL(s, k)
   call CalcLind(A) 
  do u=1,nS(s,k)   ! update LL of connected sibships
    j = SibID(u,s,k)
    if (Parent(j,3-k) < 0) then
      call CalcCLL(-Parent(j,3-k), 3-k)
    endif                    
    call CalcLind(j)
  enddo
   call calcCLL(s, k)
endif
 call CalcLind(A)
 
end subroutine RemoveSib

! #####################################################################

! @@@@   INPUT & PRECALC PROB.   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! #####################################################################

subroutine PrepData
use Global
implicit none

integer :: i,j,k,l

!=================
! allele frequencies
allocate(AF(nSNP))
do l=1,nSnp
	if (ANY(Genos(l,:)/=-9)) then
		AF(l)=float(SUM(Genos(l,:)-1, MASK=Genos(l,:)/=-9))/& !(2*nInd)
			(COUNT(Genos(l,:)/=-9)*2)
	else
		AF(l) = 1
	endif
enddo

!=================
BY1 = 1
if (ANY(BY >= 0)) then
  BY1 = MINVAL(BY, MASK=BY>=0)
  WHERE (BY /= -999) BY = BY - BY1 +1
endif

allocate(AgeDiff(nInd,nInd))
AgeDiff=999
do i=1, nInd
  do j=1, nInd
    if (BY(i)/=-999 .and. BY(j)/=-999) then
      AgeDiff(i,j) = BY(i) - BY(j)   ! if >0, then j older than i
    endif   
  enddo
enddo

  !=================   
maxAgePO = 1  ! maximum PO age difference (needed for dummy parents)
if (nAgeClasses>1) then
  do k = 2, nAgeClasses
    if (AgePriorM(k, 7)>0 .or. AgePriorM(k, 8)>0) then
      maxAgePO = k - 1  ! first k is agediff of 0
    endif
  enddo
endif

!=================
! allocate arrays
allocate(Lind(nInd))
Lind = 0
allocate(LindX(3, nSnp, nInd))
LindX = 0
allocate(FSID(MaxSibSize+1, nInd))
FSID(1, :) = (/ (i, i=1, nInd) /)
FSID(MaxSibSize+1, :) = (/ (i, i=1, nInd) /)    ! 'primary' sib                                                                                                                                    
allocate(NFS(nInd))
NFS = 1
allocate(Parent(nInd,2))
Parent = 0
allocate(DumP(3,nSnp, nInd/2,2))
DumP = 0
allocate(XPr(3,3,nSNP, nInd/2,2))
XPr = 0  
allocate(GpID(2, nInd/2,2))
GpID = 0 
 
end subroutine PrepData

! #####################################################################

subroutine DumBYrange(EstDumBY)
use Global
implicit none

integer, intent(OUT) :: EstDumBY(3,nInd/2,2)
integer :: s, k, y, CI(2), mx
double precision :: DBYP(nAgeClasses +maxAgePO), cumProp

EstDumBY = -9
do k=1,2
  do s=1,nC(k)
    call EstBY(s, k, 0, .TRUE., DBYP)
    DBYP = 10**DBYP                   
    mx = MAXLOC(DBYP, DIM=1)
    cumProp = DBYP(mx)
    CI = (/ mx, mx /)
    do y=1, nAgeClasses +maxAgePO 
      if (cumProp > 0.95) then
        EstDumBY(1,s,k) = mx  
        EstDumBY(2,s,k) = CI(1)
        EstDumBY(3,s,k) = CI(2)
        exit
      endif        
      if (CI(1) > 1 .and. CI(2) < nAgeClasses +maxAgePO) then
        if (DBYP(CI(1)-1) > DBYP(CI(2)+1)) then
          CI(1) = CI(1)-1
          cumProp = cumProp + DBYP(CI(1))
        else 
          CI(2) = CI(2)+1
          cumProp = cumProp + DBYP(CI(2))
        endif
      else if (CI(1) > 1) then
        CI(1) = CI(1)-1
        cumProp = cumProp + DBYP(CI(1))
      else if (CI(2) < nAgeClasses +maxAgePO) then
        CI(2) = CI(2)+1
        cumProp = cumProp + DBYP(CI(2))
      endif
    enddo
  enddo
enddo
EstDumBY = EstDumBY +BY1 -1 -maxAgePO

end subroutine DumBYrange

! #####################################################################

subroutine PrecalcProbs
use Global
implicit none

integer :: h,i,j,k,l,m
double precision :: OjA(3,3,nSnp), Tmp1(3), Tmp2(3,3)

!###################
allocate(AHWE(3,nSnp))
allocate(OHWE(3,nSnp))

! Prob. observed conditional on actual
!OcA(1, 1:3) = (/ 1-Er, Er/2, 0.0D0 /)   ! obs=0
!OcA(2, 1:3) = (/ Er, 1-Er, Er /)    ! obs=1
!OcA(3, 1:3) = (/ 0.0D0, Er/2, 1-Er /)   ! obs=2

!OcA(1, 1:3) = (/ 1-Er-Er**2, Er, Er**2 /)   ! obs=0
!OcA(2, 1:3) = (/ Er, 1-Er, Er /)    ! obs=1
!OcA(3, 1:3) = (/ Er**2, Er, 1-Er-Er**2 /)   ! obs=2

OcA(:, 1) = (/ 1-Er-(Er/2)**2, Er, (Er/2)**2 /)  ! act=0
OcA(:, 2) = (/ Er/2, 1-Er, Er/2 /)               ! act=1
OcA(:, 3) = (/ (Er/2)**2, Er, 1-Er-(Er/2)**2 /)  ! act=2

! probabilities actual genotypes under HWE
do l=1,nSnp
  AHWE(1,l)=(1 - AF(l))**2 
  AHWE(2,l)=2*AF(l)*(1-AF(l)) 
  AHWE(3,l)=AF(l)**2 
enddo

! joined probabilities actual & observed under HWE
do l=1,nSnp
  do i=1,3    ! obs
    do j=1,3    ! act
      OjA(i, j, l) = OcA(i,j) * AHWE(j, l)
    enddo
  enddo
enddo

! marginal prob. observed genotypes
do l=1,nSnp
  do i=1,3
    OHWE(i, l) = SUM(OjA(i, 1:3, l))
  enddo
enddo

! ########################
! inheritance conditional on 1 parent
allocate(AKAP(3,3,nSnp))
allocate(OKAP(3,3,nSnp))
allocate(OKOP(3,3,nSnp))

do l=1,nSnp
  AKAP(1, 1:3, l) = (/ 1-AF(l), (1-AF(l))/2, 0.0D0 /)
  AKAP(2, 1:3, l) = (/ AF(l), 0.5D0, 1-AF(l) /)
  AKAP(3, 1:3, l) = (/ 0D0, AF(l)/2, AF(l) /)
enddo

do l=1,nSnp
  do i=1,3  ! obs offspring
    do j=1,3    ! act parent
      Tmp1=0
      do k=1,3    ! act offspring
        Tmp1(k) = OcA(i,k) * AKAP(k,j,l)
      enddo
      OKAP(i,j,l) = SUM(Tmp1)
    enddo
  enddo
enddo

do l=1,nSnp
  do i=1,3  ! obs offspring
    do j=1,3    ! obs parent
      Tmp2=0
      do k=1,3    ! act offspring
        do m=1,3    ! act parent
          Tmp2(k,m) = OcA(i,k) * OcA(j,m) * AKAP(k,m,l)
        enddo
      enddo
      OKOP(i,j,l) = SUM(Tmp2) 
    enddo
  enddo
enddo

! #########################
! inheritance conditional on both parents

AKA2P(1,1,:) = (/ 1.0, 0.5, 0.0 /)
AKA2P(1,2,:) = (/ 0.5, 0.25, 0.0 /)
AKA2P(1,3,:) = (/ 0.0, 0.0, 0.0 /)

AKA2P(2,1,:) = (/ 0.0, 0.5, 1.0 /)
AKA2P(2,2,:) = (/ 0.5, 0.5, 0.5 /)
AKA2P(2,3,:) = (/ 1.0, 0.5, 0.0 /)

AKA2P(3,1,:) = (/ 0.0, 0.0, 0.0 /)
AKA2P(3,2,:) = (/ 0.0, 0.25, 0.5 /)
AKA2P(3,3,:) = (/ 0.0, 0.5, 1.0 /)

do i=1,3  ! obs offspring
  do j=1,3    ! act parent 1
    do h=1,3    !act parent 2
      Tmp1=0
      do k=1,3    ! act offspring
        Tmp1(k) = OcA(i,k) * AKA2P(k,j,h) 
      enddo
      OKA2P(i,j,h) = SUM(Tmp1)
    enddo
  enddo
enddo

!=================
allocate(PHS(3,3,nSnp))
allocate(PFS(3,3,nSnp))
do l=1,nSnp
  do i=1,3  ! obs offspring 1
    do j=1,3    ! obs offspring 2
      Tmp1=0
      Tmp2=0
      do m=1,3    !act shared parent 
        Tmp1(m) = OKAP(i,m,l) * OKAP(j,m,l) * AHWE(m,l)
        do h=1,3
          Tmp2(m,h) = OKA2P(i,m,h) * OKA2P(j,m,h) * AHWE(m,l) *AHWE(h,l)
        enddo
      enddo
      PHS(i,j,l) = SUM(Tmp1) / (AHWE(i,l) * AHWE(j,l))
      PFS(i,j,l) = SUM(Tmp2) / (AHWE(i,l) * AHWE(j,l))
    enddo
  enddo
enddo

allocate(LindG(3, nSnp, nInd))  ! used when missing genotype
do l=1,nSnp
  do i=1,3
    LindG(i,l,:) = AHWE(i,l)  
  enddo
enddo

deallocate(AF)

end subroutine PrecalcProbs

! #####################################################################

subroutine deallocall
use Global
implicit none

! allocated in PrepData
if (allocated(Sex)) deallocate(Sex)
if (allocated(BY)) deallocate(BY)
if (allocated(PairType)) deallocate(PairType)
if (allocated(nFS)) deallocate(nFS)

if (allocated(Genos)) deallocate(Genos)
if (allocated(AgeDiff)) deallocate(AgeDiff)
if (allocated(Parent)) deallocate(Parent)
if (allocated(OppHomM)) deallocate(OppHomM)
if (allocated(nS)) deallocate(nS)
if (allocated(PairID)) deallocate(PairID)
if (allocated(FSID)) deallocate(FSID)

if (allocated(SibID)) deallocate(SibID)
if (allocated(GpID)) deallocate(GpID)

if (allocated(Lind)) deallocate(Lind)
if (allocated(PairDLLR)) deallocate(PairDLLR)
if (allocated(AF)) deallocate(AF)

if (allocated(AHWE)) deallocate(AHWE)
if (allocated(OHWE)) deallocate(OHWE) 
if (allocated(LLR_O)) deallocate(LLR_O)
if (allocated(LindX)) deallocate(LindX)
if (allocated(LR_parent)) deallocate(LR_parent)
if (allocated(AgePriorM)) deallocate(AgePriorM)
if (allocated(CLL)) deallocate(CLL)

if (allocated(AKAP)) deallocate(AKAP)
if (allocated(OKAP)) deallocate(OKAP)
if (allocated(OKOP)) deallocate(OKOP)
if (allocated(LR_GP)) deallocate(LR_GP)
if (allocated(LindG)) deallocate(LindG)
if (allocated(PHS)) deallocate(PHS)
if (allocated(PFS)) deallocate(PFS)
if (allocated(DumBY)) deallocate(DumBY)
if (allocated(DumP)) deallocate(DumP)
if (allocated(XPr)) deallocate(XPr)

!if (allocated(BYRank)) deallocate(BYRank)
!if (allocated(SortBY)) deallocate(SortBY)

end subroutine deallocall

! #####################################################################

! -9   NA
! 999  NA
! 888  Already assigned
! 777  impossible
! 444  not yet implemented (typically involves inbreeding)
! 222  as likely to go via opposite parent