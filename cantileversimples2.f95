program cantilever1forca
   
   implicit none
    
    double precision, dimension (1,2) :: cc
    double precision, dimension (2,11) ::  mat_transf
    double precision, dimension (1,3) :: prop
    double precision, dimension (11,2) :: transmat
    double precision :: a, e, l, a1,a2
    integer :: i, ngdl, nele
    double precision, dimension (2,2) :: kop, kele
    integer, dimension (10,2) :: fa
    double precision, dimension (10,10) :: invkll, kll
    integer :: n
    double precision :: coef
    double precision, dimension (11,11) :: kglobal, kglobal1
    integer, dimension (11,11) :: ident
    integer, dimension (10,2) ::  mat_conect
    double precision, dimension (10) :: u_liv
    double precision, dimension (10) :: flivres
    integer, dimension (10) :: gdl_livres
    integer, dimension (1) :: gdl_impostos
    
    a1 = 0.05
    a2 = 0.1
    a = 0.05
    e = 10*10**9
    l = 0.2 !mm
    ngdl = 11
    nele = 10
    
    cc(1,:) = (/1,0/)
    
    fa(1,:) = (/2,0/)
    fa(2,:) = (/3,0/)
    fa(3,:) = (/4,0/)
    fa(4,:) = (/5,0/)
    fa(5,:) = (/6,0/)
    fa(6,:) = (/7,0/)
    fa(7,:) = (/8,0/)
    fa(8,:) = (/9,0/)
    fa(9,:) = (/10,0/)
    fa(10,:) = (/11,-10000/)
    
    prop(1,:) = (/a, e, l/)
    
    
    mat_conect(1,:) = (/1,2/)
    mat_conect(2,:) = (/2,3/)
    mat_conect(3,:) = (/3,4/)
    mat_conect(4,:) = (/4,5/)
    mat_conect(5,:) = (/5,6/)
    mat_conect(6,:) = (/6,7/)
    mat_conect(7,:) = (/7,8/)
    mat_conect(8,:) = (/8,9/)
    mat_conect(9,:) = (/9,10/)
    mat_conect(10,:) = (/10,11/)
    
    kop (1,:) = (/1,-1/)
    kop (2,:) = (/-1,1/)
    
    ident = 0
    do i=1,ngdl
        ident(i,i) = 1
    end do
    kglobal = 0
    kglobal1 = 0
    kele = 0
    
    do i=1,nele
                coef = prop(1,1)*prop(1,2)/prop(1,3)
                kele = coef*kop
                mat_transf (1,:) = ident(mat_conect(i,1),:)
                mat_transf (2,:) = ident(mat_conect(i,2),:)
                transmat = transpose(mat_transf)
                call matmu(transmat,kele,mat_transf,kglobal1)
                kglobal = kglobal+kglobal1
            end do
       
    gdl_livres = fa(:,1)
   
     gdl_impostos = cc(:,1)
    
     kll = kglobal(gdl_livres,gdl_livres)
!      kli = kglobal(gdl_livres,gdl_impostos)
     
     
     n = 10
     call inv(kll, invkll, n) 
     flivres = fa(:,2)
     u_liv = matmul(invkll,fa(:,2))
    
     print*, u_liv
        
     
    end program
    
    subroutine matmu (A,B,C,D)
      !programa para multiplicacao de 3 matrizes
  
  double precision, dimension (11,2) :: A, mul
  double precision, dimension (2,2) :: B
  double precision, dimension (2,11) :: C
  double precision, dimension (11,11) :: D
  
  mul = matmul(A,B)
  D = matmul (mul,C)
  
    end subroutine matmu
    subroutine inv(a,c,n)
!============================================================
! Inverse matrix
! Method: Based LU decomposition for Ax=b
! Marouf Abderahmane - Strasbourg University 2016 
!-----------------------------------------------------------
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

L=0.0
U=0.0
b=0.0

do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

do i=1,n
  L(i,i) = 1.0
end do
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inv
             
    
