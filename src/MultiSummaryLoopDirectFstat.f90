!
!	MultiSummaryLoopDirectFstat.f90	
!
!	Created by Xia Shen on 22/03/18.
!

subroutine MultiSummaryLoopDirectFstat(k, m, nn, betamat, invsemat, invse, invR, invV, pil, fstat)

  implicit none
  integer :: i, j, k, m
  double precision, dimension(m, m) :: invR, invse, invV
  double precision, dimension(1, m) :: beta
  double precision, dimension(k) :: nn, pil, fstat
  double precision, dimension(k, m) :: betamat, invsemat
  double precision, dimension(1, 1) :: t2
  do j = 1, k
     do i = 1, m
        beta(1,i) = betamat(j,i)
        invse(i,i) = invsemat(j,i)
     end do
     invV = matmul(matmul(invse, invR), invse)
     t2 = matmul(matmul(beta, invV), transpose(beta))
     pil(j) = t2(1,1)
     fstat(j) = pil(j)/m/(nn(j) - pil(j))*(nn(j) - m - 1)
  end do

end subroutine MultiSummaryLoopDirectFstat

