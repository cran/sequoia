subroutine getrel(nind, pedrf, nrel, relv)
implicit none

integer, intent(IN) :: nind, nrel
integer, intent(INOUT) :: relv(nind*nind*nrel)  ! 0/1 matrix
integer, intent(IN) :: pedrf(nInd*2)
logical :: doGP 
integer :: i, j, x, y, r, GPr(2,2)
integer, allocatable :: ped(:,:), rel(:,:,:)



! relationships:
! 1 = self
! 2 = dam
! 3 = sire
! 4 = offspring

! 5 = full sib
! 6 = maternal half sib
! 7 = paternal half sib
! 8 = 'cross' half sib (hermaphrodite)

! 9  = MGM
! 10 = MGF
! 11 = PGM
! 12 = PGF
! 13 = GO

! 14 = FA  + dHA
! 15 = FN
! 16 = HA
! 17 = HN
! 18 = dFC1
! 19 = FC1

GPr(1,1) = 9
GPr(1,2) = 10
GPr(2,1) = 11
GPr(2,2) = 12

! fold pedigree
allocate(ped(nInd,2))
ped(:,1) = pedrf(1:nInd)
ped(:,2) = pedrf((nInd+1) : (2*nInd))

if (nrel == 8) then
  doGP = .FALSE.
else
  doGP = .TRUE.
endif

allocate(rel(nInd, nInd, nrel))
rel = 0
do i = 1, nInd
  rel(i,i,1) = 1   ! self
  if (ped(i,1) /= 0) then
    rel(i, ped(i,1), 2) = 1  ! dam
    rel(ped(i,1), i, 4) = 1  ! offspring
  endif
  
  if (ped(i,2) /= 0) then
    rel(i, ped(i,2), 3) = 1  ! sire
    rel(ped(i,2), i, 4) = 1  ! offspring
  endif
enddo

! sibs  
do i = 1, nInd
  if (all(ped(i,:) == 0))  cycle
  
  do j = i+1, nInd
    if (all(ped(j,:) == 0))  cycle
    
    do x=1,2
      if (ped(i,x) == ped(j,x) .and. ped(i,x)/=0) then
        rel(i,j,5+x) = 1  ! mat/pat HS
      endif
    enddo
    if (rel(i,j,6) == 1 .and. rel(i,j,7) == 1) then    ! full sibs
      rel(i,j,5) = 1
      rel(i,j, 6:7) = 0
    endif
    
    ! hermaphrodites
    do x=1,2
      if (ped(i,x) == ped(j,3-x) .and. ped(i,x)/=0) then
        if (rel(i,j,8) == 1) then
          rel(i,j,5) = 1  ! cross FS  
          rel(i,j,8) = 0
        else
          rel(i,j,8) = 1  ! cross HS
        endif
      endif
    enddo
    
    do r=5,8
      rel(j,i,r) = rel(i,j,r)
    enddo

  enddo
enddo

! grandparents
if (doGP) then
  do i = 1, nInd
    do x = 1,2  
      if (ped(i,x) == 0)  cycle
      do y = 1,2
        j = ped (ped(i,x), y)  ! ID number of GP
        if (j == 0)  cycle
        rel(i,j, GPr(x,y) ) = 1  
        rel(j, i, 13) = 1  ! GO 
      enddo
    enddo
  enddo
endif

! avuncular   
if (doGP) then
  do i = 1, nInd
    do x = 1,2
      if (ped(i,x) == 0)  cycle
      
      do j = 1, nInd
        if (rel(ped(i,x), j, 5) == 1) then
          rel(i,j,14) = 1  ! FA
          rel(j,i,15) = 1  ! FN
        else if (any(rel(ped(i,x), j, 6:8) == 1)) then
          if (x==2 .and. rel(i,j,16) == 1) then
            rel(i,j,14) = 1  ! double HA
            rel(j,i,15) = 1  
          endif
          rel(i,j,16) = 1  ! HA
          rel(j,i,17) = 1  ! HN
        endif
        
        do y = 1,2
          if (ped(j,y) == 0)  cycle
          if (rel(ped(i,x), ped(j,y), 5) == 1) then
            if (rel(i,j,19) == 1)  rel(i,j,18) = 1  ! double full cousins
            rel(i,j,19) = 1  ! full cousins
          endif
        enddo
      enddo
    enddo
  enddo
endif

! unfold rel
relV = 0
do r = 1,nrel
  do j = 1,nInd
    do i = 1,nInd  
      x = ((r-1)*nInd + (j-1))*nInd + i
      relV(x) = rel(i,j,r)
    enddo
  enddo
enddo

deallocate(ped)
deallocate(rel)


end subroutine getrel
