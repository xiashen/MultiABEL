!
!	MultiSummaryLoop.f90	
!
!	Created by Xia Shen on 06/04/16.

subroutine MultiSummaryLoop(k, m, nn, f, betamat, R, invR, D, sdY, invsdY, sY, b, s, pil, coef, ss)

  implicit none
  integer :: ii, i, j, k, m
  double precision, dimension(m, m) :: sdY, invsdY, R, invR, H, HinvE0, invE0, D, A, invA
  double precision, dimension(m - 1, m - 1) :: zz, sdY1, R1
  double precision, dimension(m, 1) :: zg, beta, ab, alphabeta, sY, tmp2, bb
  double precision, dimension(m - 1, 1) :: beta1, alpha, sY1, cor1
  double precision, dimension(1, m - 1) :: tbeta1  
  double precision, dimension(k) :: nn, f, pil
  double precision, dimension(k, m) :: betamat, b, s, coef, ss
  double precision, dimension(m) :: lambda
  double precision, dimension(1, 1) :: abDalphabeta0, Vb, bDbeta0
  double precision :: tmp, vg, abDalphabeta, bDbeta, gg
  do j = 1, k
     vg = 2*f(j)*(1 - f(j))
     do i = 1, m
        beta(i,1) = betamat(j,i)
     end do
     do i = 1, m
        if (i == 1) then
           R1 = R(2:m,2:m)
           sdY1 = sdY(2:m,2:m)
           beta1(1:(m - 1),1) = beta(2:m,1)
           sY1(1:(m - 1),1) = sY(2:m,1)
           cor1(1:(m - 1),1) = R(2:m,1)
        else if (i == m) then
           R1 = R(1:(m - 1),1:(m - 1))
           sdY1 = sdY(1:(m - 1),1:(m - 1))
           beta1(1:(m - 1),1) = beta(1:(m - 1),1) 
           sY1(1:(m - 1),1) = sY(1:(m - 1),1)
           cor1(1:(m - 1),1) = R(1:(m - 1),m) 
        else
           R1(1:(i - 1),1:(i - 1)) = R(1:(i - 1),1:(i - 1))
           R1(i:(m - 1),i:(m - 1)) = R((i + 1):m,(i + 1):m)
           R1(1:(i - 1),i:(m - 1)) = R(1:(i - 1),(i + 1):m)
           R1(i:(m - 1),1:(i - 1)) = R((i + 1):m,1:(i - 1))
           sdY1(1:(i - 1),1:(i - 1)) = sdY(1:(i - 1),1:(i - 1))
           sdY1(i:(m - 1),i:(m - 1)) = sdY((i + 1):m,(i + 1):m)
           sdY1(1:(i - 1),i:(m - 1)) = sdY(1:(i - 1),(i + 1):m)
           sdY1(i:(m - 1),1:(i - 1)) = sdY((i + 1):m,1:(i - 1))
           beta1(1:(i - 1),1) = beta(1:(i - 1),1)
           beta1(i:(m - 1),1) = beta((i + 1):m,1) 
           sY1(1:(i - 1),1) = sY(1:(i - 1),1)
           sY1(i:(m - 1),1) = sY((i + 1):m,1) 
           cor1(1:(i - 1),1) = R(1:(i - 1),i)
           cor1(i:(m - 1),1) = R((i + 1):m,i) 
        end if
        zg = beta*vg
        zz = matmul(matmul(sdY1, R1), sdY1)
        tmp = sdY(i,i)
        do ii = 1, m - 1
           alpha(ii,1) = cor1(ii,1)/sY1(ii,1)*tmp
        end do

        alphabeta(1:(m - 1),1) = alpha(1:(m - 1),1)
        alphabeta(m,1) = beta(i,1)
        A(1:(m - 1),1:(m - 1)) = zz(1:(m - 1), 1:(m - 1))
        A(1:(m - 1),m) = beta1(1:(m - 1), 1)*vg
        tbeta1 =  transpose(beta1)
        A(m,1:(m - 1)) = tbeta1(1,1:(m - 1))*vg
        A(m,m) = vg
        do ii = 1, m
           D(ii,ii) = A(ii,ii)
        end do
        invA = inverse(A, m)
        ab = matmul(matmul(invA, D), alphabeta)
        abDalphabeta0 = matmul(matmul(transpose(ab), D), alphabeta)
        abDalphabeta = abDalphabeta0(1,1)
        Vb = (sdY(i,i)**2 - abDalphabeta)/(nn(j) - 3)*invA(m, m)  
        b(j,i) = ab(m,1)
        s(j,i) = sqrt(Vb(1,1))
     end do

     tmp2 = matmul(invR, beta)*vg
     do i = 1, m
        coef(j,i) = tmp2(i,1)
        bb(i,1) = tmp2(i,1)
        D(i,i) = sY(i,1)**2
     end do
     gg = vg*nn(j)
     bDbeta0 = matmul(matmul(transpose(bb), D), beta)
     bDbeta = bDbeta0(1,1)*vg
     invE0 = matmul(invsdY, matmul(invR, invsdY))
     do i = 1, m
        ss(j,i) = (gg - bDbeta)/(nn(j) - m)*invE0(i,i)/nn(j)
     end do
 
     H = matmul(beta, transpose(beta))*nn(j)*vg
     HinvE0 = matmul(matmul(H, invsdY), matmul(invR, invsdY))/nn(j)
     do i = 1, m
        lambda(i) = HinvE0(i,i)
     end do
     pil(j) = sum(lambda)
  end do

contains

function inverse(a, n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer :: n
double precision :: a(n,n)
double precision :: c(n,n), inverse(n, n)
double precision :: L(n,n), U(n,n), b(n), d(n), x(n)
double precision :: coeff
integer :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
inverse = c
end function inverse

end subroutine MultiSummaryLoop


