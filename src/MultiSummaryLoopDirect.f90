!
!	MultiSummaryLoopDirect.f90	
!
!	Created by Xia Shen on 25/02/16.
!

subroutine MultiSummaryLoopDirect(k, m, nn, tmat, invR, pil, fstat)

  implicit none
  integer :: i, j, k, m
  double precision, dimension(m, m) :: invR
  double precision, dimension(1, m) :: t
  double precision, dimension(k) :: nn, pil, fstat
  double precision, dimension(k, m) :: tmat
  double precision, dimension(1, 1) :: t2
  do j = 1, k
     do i = 1, m
        t(1,i) = tmat(j,i)
     end do
     t2 =  matmul(matmul(t, invR), transpose(t))
     pil(j) = t2(1,1)
     fstat(j) = pil(j)*(nn(j) - m - 1)/m/(nn(j) - 1)
  end do

end subroutine MultiSummaryLoopDirect

