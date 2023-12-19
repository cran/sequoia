! check if the pedigree is compatible with the genetic data
! pedigree: animal_id, dam_id, sire_id, mat_grandsire_id
! genotype: 1 row per individual, 1st column is ID, then 1 column per SNP; 0/1/2 copies of minor allele, -9=missing
! other parameters: runtime options (a.o. genotyping error)

! It's an extremely minimalistic version of sequoia, which should hopefully be useful for finding pedigree errors in both the nucleus herd as in the multiplication herds, as well as 

! compile: gfortran -std=f95 -fall-intrinsics -O3 pedigree_checker.f90 -o PedChecker

! ##############################################################################
! ##  Global variables  ##
! ##############################################################################

module ester_global_vars
  implicit none

  integer :: nInd, nSnp, ID_len
  integer, parameter :: maxnOff = 50
  integer,allocatable,dimension(:,:) :: Genos, Parent
  integer,allocatable,dimension(:) :: DupID
  double precision :: OcA(0:2,-1:2), OKA2P(-1:2,0:2,0:2), AKA2P(0:2,0:2,0:2)
  double precision, allocatable, dimension(:,:) :: AHWE
  double precision, allocatable, dimension(:,:,:,:) :: PrG

end module ester_global_vars


! ==============================================================================
! ##   input/output  ##
! ==============================================================================

subroutine ester(ng, nl, genov, parentsv, dups, errin, totll, cntobsact)
  use ester_global_vars
  implicit none
  
  integer, intent(IN) :: ng, nl 
  integer, intent(IN) :: genov(ng*nl), parentsv(2*ng), dups(ng)
  double precision, intent(IN) :: errin(3)
  double precision, intent(INOUT) :: totLL, cntObsAct(3*3*nl)
  integer:: i,j,l
  double precision :: totLL_tmp, PropsOA(3,3,nl)
  
  nInd = ng   ! else input variable cannot be global
  nSnp = nl  
  
  ! ! fold input vectors to matrices
  allocate(Genos(nSnp, nInd))
  Genos = -1
  j = 0
  do l=1,nSnp
    do i=1, nInd
      j = j+1
      if (Genov(j)/=-9) then
        Genos(l,i) = Genov(j) 
      endif
    enddo
  enddo
  
  allocate(Parent(2, nInd))   ! dam - sire
  Parent(1,:) = parentsv(1:nInd)
  parent(2,:) = parentsv((nInd+1) : (2*nInd))
  
  allocate(dupID(nInd))
  dupID = dups

  ! ! prep
  allocate(AHWE(3,nSnp))   
  call ester_CalcProbs(Errin)
    
  ! Estimate actual genotypes
  allocate(PrG(0:2,3,nSnp,0:nInd))  ! 0/1/2 - anc/desc/own - snp - ind
  PrG = 0D0
  call EstAllG()
  
  ! Calculate likelihood
  TotLL_tmp = 999D0
  call ester_CalcLL(TotLL_tmp, PropsOA)   ! , PropObsAct_array  causes crash??
  TotLL = TotLL_tmp
!  call dblepr('totLL:', 10, (/TotLL/), 1)
  
  
  ! fold array to vector to pass back to R
   cntobsact = RESHAPE(PropsOA, (/3*3*nSnp/) )  ! , ORDER=(/3,2,1/)
  
  deallocate(Genos)
  deallocate(Parent) 
  deallocate(DupID)
  deallocate(AHWE)
  deallocate(PrG)

end subroutine ester


! ==============================================================================
! ##   Prep  ##
! ==============================================================================

! Calc obs|act & Mendelian inheritance probs
subroutine ester_CalcProbs(Er)
  use ester_global_vars
  implicit none

  double precision, intent(IN) :: Er(3)
  integer :: h,i,j,k,l
  double precision :: AF(nSnp), Tmp(0:2)
  logical :: MSK(nInd), dupExclude(nInd)

  ! Estimated allele frequency, without correction for genotyping errors
  dupExclude = .FALSE.
  do i=1,nInd
    if (dupID(i)/=0 .and. dupID(i) < i)  dupExclude(i) = .TRUE.
  enddo
  
  do l=1,nSnp
    MSK = Genos(l,:)/=-1 .and. .not. dupExclude
    if (ANY(MSK)) then      
      AF(l)=float(SUM(Genos(l,:), MASK=MSK))/(COUNT(MSK)*2)
    else
      AF(l) = 1D0
    endif
  enddo


  ! Prob. observed (d2) conditional on actual (d1) 
  ! 'ErrFlavour' = 2.0
  OcA(:,-1) = 1.0D0      ! missing 
  OcA(0, 0:2) = (/ 1-Er(1)-Er(2), Er(2), Er(1) /)   
  OcA(1, 0:2) = (/ Er(3), 1-2*Er(3), Er(3) /)      
  OcA(2, 0:2) = (/ Er(1), Er(2), 1-Er(1)-Er(2) /)  

  
  ! probabilities actual genotypes under HWE
  do l=1,nSnp
    AHWE(1,l)=(1 - AF(l))**2 
    AHWE(2,l)=2*AF(l)*(1-AF(l)) 
    AHWE(3,l)=AF(l)**2 
  enddo
  

  ! inheritance conditional on both parents
  AKA2P(0,0,:) = dble((/ 1.0, 0.5, 0.0 /))
  AKA2P(0,1,:) = dble((/ 0.5, 0.25, 0.0 /))
  AKA2P(0,2,:) = dble((/ 0.0, 0.0, 0.0 /))

  AKA2P(1,0,:) = dble((/ 0.0, 0.5, 1.0 /))
  AKA2P(1,1,:) = dble((/ 0.5, 0.5, 0.5 /))
  AKA2P(1,2,:) = dble((/ 1.0, 0.5, 0.0 /))

  AKA2P(2,0,:) = dble((/ 0.0, 0.0, 0.0 /))
  AKA2P(2,1,:) = dble((/ 0.0, 0.25, 0.5 /))
  AKA2P(2,2,:) = dble((/ 0.0, 0.5, 1.0 /))


  do i=-1,2  ! obs offspring
    do j=0,2    ! act parent 1
      do h=0,2    !act parent 2
        Tmp=0D0
        do k=0,2    ! act offspring
          Tmp(k) = OcA(k,i) * AKA2P(k,j,h) 
        enddo
        OKA2P(i,j,h) = SUM(Tmp)
      enddo
    enddo
  enddo

end subroutine ester_CalcProbs


! ==============================================================================
! ##    Main program     ##
! ##   Calc likelihood   ##
! ==============================================================================

subroutine ester_CalcLL(TotLL, PropObsAct)  
  use ester_global_vars
  implicit none

  double precision, intent(OUT) :: TotLL, PropObsAct(0:2,0:2,nSnp)
  integer :: i, l
  double precision:: PrG_par(0:2), PrG_off(0:2), PrG_dup(0:2), PrG_act(0:2), &
     CntObsAct(0:2,0:2,nSnp)  ! actual - observed - snp
  double precision, allocatable :: LL(:,:)
  
  allocate(LL(nSnp, nInd))
  LL = 0D0  
  CntObsAct = 0D0
  
  do i=1, nInd
!    if (DupID(i)/=0 .and. DupID(i) < i)  cycle  ! already done as duplicate-of?
! TODO: properly do duplicates + missing geno 
    
    do l=1,nSnp
      if (Genos(l,i)==-1)  cycle  ! no observed genotype --> no (direct) contribution to LL     
    
      PrG_par = PrG(:,1,l,i)
      PrG_off = PrG(:,2,l,i)
      
      if (DupID(i)/=0) then
        PrG_dup = OcA(:,Genos(l,DupID(i)))
      else
        PrG_dup = 1D0 !/3
      endif
      
      PrG_act = PrG_par * PrG_off * PrG_dup
      PrG_act = PrG_act / SUM(PrG_act)
      
      CntObsAct(:, Genos(l,i),l) = CntObsAct(:, Genos(l,i),l) + PrG_act
      LL(l,i) = LOG10(SUM(OcA(:,Genos(l,i)) * PrG_act))
  
    enddo
  enddo
  
  TotLL = SUM(LL)
  
  do l=1, nSnp
    PropObsAct(:,:,l) = CntObsAct(:,:,l) / SUM(CntObsAct(:,:,l))
  enddo
  
  deallocate(LL)

end subroutine ester_CalcLL


! ==============================================================================
! ##  Calculate genotype probabilities  ##
! ==============================================================================


! ==============================================================================
! estimate actual genotype for all individuals, from founders to latest generation & vv

subroutine EstAllG
  use ester_global_vars 
  implicit none

  integer :: l, i, g, Gen(nInd), maxG
  
  call ester_getGenerations(Gen)
  maxG = MAXVAL(Gen, DIM=1)
  
  PrG = 1D0
  do l=1, nSnp
    PrG(:,1,l,0) = AHWE(:,l)
  enddo

  ! top -> bottom
  do g=0, maxG
    do i=1, nInd
      if (Gen(i) /= g)  cycle
      call EstG(i)
    enddo
  enddo

end subroutine EstAllG


! ==============================================================================
! probability of own actual genotype, conditional on 
! 1: parents & ancestors
! 2: offspring & descendants
! 3: own observed genotype

subroutine EstG(i)
  use ester_global_vars 
  implicit none
  
  integer, intent(IN) :: i  
  integer :: l, k, h, nOff, IDoff(maxnOff), IDmate(maxnOff), sexi, x
  double precision :: PrG_par(0:2), PrG_off(0:2), PrG_mate(0:2), PrTmp(0:2), PrG_self(0:2)
  
  call ester_getOff(i, nOff, IDoff, IDmate, sexi)
  
  do l=1, nSnp  

    call ester_ParProb(l, i, 1, PrG_par)   
    
    PrG_off = 1D0  ! TODO: separate subroutine
    if (nOff > 0) then
      do h=1,nOff
        call ester_ParProb(l, IDmate(h), 13, PrG_mate)
        ! may have multiple offspring with same mate
        do k=1, nOff
          if (IDmate(h) == 0 .and. h/=k)  cycle  ! unknown parent, assume all different
          if (IDmate(k) /= IDmate(h))  cycle     ! different parent
          if (IDmate(k) == IDmate(h) .and. IDmate(h) /= 0 .and. h > k)  cycle   ! already done, with a full sib
        
          do x=0,2
            PrTmp(x) = SUM(OKA2P(Genos(l,IDoff(k)), x, :) * PrG_mate)
          enddo        
          PrG_off = PrG_off * PrTmp
          ! TODO: descendants
        enddo
      enddo
    endif
    
    PrG_self = OcA(:,Genos(l,i))
    if (DupID(i)/=0) then
      PrG_self = PrG_self * OcA(:,Genos(l, DupID(i)))
    endif
    
    PrG(:,1,l,i) = PrG_par   
    PrG(:,2,l,i) = PrG_off 
    PrG(:,3,l,i) = PrG_self 
    
    do x=1,3
      PrG(:,x,l,i) = PrG(:,x,l,i) / sum(PrG(:,x,l,i))
    enddo
    
  enddo

end subroutine EstG


! ==============================================================================
! genotype probabilities conditional on parent genotype or allele frequency + HWE

subroutine ester_ParProb(l, i, typ, prob)
  use ester_global_vars
  implicit none
  
  integer, intent(IN) :: l, i, typ
  double precision, intent(OUT) :: prob(0:2)
  integer :: k, x,y
  double precision :: PrPar(0:2,2), PrTmp(0:2)
  
  if (i == 0) then  ! no parent
    prob = AHWE(:,l)  
  
  else if (i > 0) then
    if (typ == 3) then  ! only own genotype
      prob = OcA(:,Genos(l,i)) * AHWE(:,l)  
    else if (typ == 13 .or. typ == 1) then  ! own genotype + parents  / only parents
      do k=1,2
        PrPar(:,k) = PrG(:,1,l, Parent(k,i)) * PrG(:,3,l, Parent(k,i))
      enddo
      do x=0,2
        do y=0,2
          PrTmp(y) = SUM( AKA2P(x,y,:) * PrPar(y,1) * PrPar(:,2) )
        enddo
        prob(x) = SUM(PrTmp)
      enddo
      if (typ == 13) then
        prob = OcA(:,Genos(l,i)) * prob
      endif
    else
      call rexit("  ERROR! *** parprob typ not yet implemented *** ")
    endif
    
    prob = prob / sum(prob)
  endif

end subroutine ester_ParProb


! ==============================================================================
! ##  Helper subroutines  ##
! ==============================================================================


! ==============================================================================
! get offspring ids for an individual 

subroutine ester_getOff(i, nOff, IDoff, IDmate, sexi)
  use ester_global_vars
  implicit none

  integer, intent(IN) :: i
  integer, intent(OUT) :: nOff, IDoff(maxnOff), IDmate(maxnOff), sexi
  integer :: j,k
  
  nOff = 0
  IDoff = 0
  IDmate = 0   ! unique mates
  sexi = 3
  
  if (i==0)  return
  
  do j = 1, nInd
    do k=1,2
      if (Parent(k,j) == i) then
      
        if (sexi == 3) then
          sexi = k 
        else if (sexi /= k) then
          call rexit("  ERROR! *** Individual is both dam and sire! *** ")
        endif
       
        nOff = nOff +1
        IDoff(nOff) = j
        IDmate(nOff) = Parent(3-sexi, j)
        
        if (nOff == maxnOff) then
          call rexit("  ERROR! *** reached maximum number of offspring, please increase 'maxnOff' *** ")
        endif
      endif  
    enddo
  enddo

end subroutine ester_getOff


! ==============================================================================
! get generation number for each individual

subroutine ester_getGenerations(Gen)
use ester_global_vars
implicit none

integer, intent(OUT) :: Gen(nInd)
integer :: i, GenPar(2, nInd), g, Prevgen(nInd), m

Gen = -9
GenPar = 0   ! generation number of dam, sire
Prevgen = 0   ! ids of individuals in generations up to g-1

do i=1,nInd
  if (ALL(Parent(:,i)==0)) then
    Gen(i) = 0  ! founder
    Prevgen(i) = i
  endif
enddo

do g = 0, 1000
  do i=1, nInd
    if (Gen(i) >= 0)  cycle
    do m=1,2
      if (Parent(m,i)==0 .or. GenPar(m,i) > 0)  cycle
      if (ANY(Prevgen == Parent(m,i))) then
        GenPar(m,i) = g
      endif
    enddo
    
    if (ALL(GenPar(:,i) <= g)) then  ! including Parent(m,i)==0 --> GenPar(m,i)==0
      Gen(i) = g+1
      PrevGen(i) = i
    endif    
  enddo
  if (.not. ANY(Gen < 0))  exit   ! all done
enddo

end subroutine ester_getGenerations

!===============================================================================
! end.