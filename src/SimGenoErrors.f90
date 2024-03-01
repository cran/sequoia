! Array storage in Fortran is column-major:
! the first index varies fastest

subroutine mkerrors(nind, nsnp, genofr, eprobfr, randomv)
implicit none

integer, intent(IN) :: nind, nsnp
integer, intent(INOUT) :: genofr(nind*nsnp)
double precision, intent(IN) :: eprobfr(nsnp*3*3), randomv(nsnp*nInd)
integer :: l, i, x, j, h
integer, allocatable :: Genos(:,:)
double precision :: EProb(3,3,nSnp), p(3)
double precision, allocatable :: r(:,:)

allocate(Genos(nSnp,nInd))
allocate(r(nSnp, nInd))
Genos = -9
j = 0
r = 0D0
do l=1,nSnp
  do i=1, nInd
    j = j+1
    if (GenoFR(j)>=0) then
      Genos(l,i) = GenoFR(j) +1
    endif
    r(l,i) = randomv(j)
  enddo
enddo

j=0
do i=1,3
  do h=1,3
    do l=1,nSnp
      j = j+1
      EProb(h,i,l) = EProbFR(j)
    enddo
  enddo
enddo

j=0
do l=1, nSnp
  do i=1,nInd
    if (Genos(l,i)==-9)  cycle
    p = Eprob(Genos(l,i), :, l) / sum(Eprob(Genos(l,i), :, l))
    if (r(l,i) < p(1)) then
      x = 0
    else if (r(l,i) < (p(1) + p(2))) then
      x = 1
    else
      x = 2
    endif
    
    j = j+1
    GenoFR(j) = x
  enddo
enddo

deallocate(Genos)
deallocate(r)

end subroutine MkErrors
