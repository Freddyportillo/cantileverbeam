program cantilever1forca
   
   implicit none
    
    double precision, dimension (2,2) :: cc
    double precision, dimension (4,22) ::  mat_transf
    double precision, dimension (1,4) :: prop
    double precision, dimension (22,4) :: transmat
    double precision :: a, e, l, i
    integer :: ii, ngdl, nele, nnodes
    double precision, dimension (4,4) :: kop, kele
    integer, dimension (20,2) :: fa
    double precision, dimension (20,20) :: invkll, kll
    integer :: n
    double precision :: coef
    double precision, dimension (22,22) :: kglobal, kglobal1
    integer, dimension (22,22) :: ident
    integer, dimension (20,4) ::  mat_conect
    double precision, dimension (20) :: u_liv
    double precision, dimension (20) :: flivres
    integer, dimension (20) :: gdl_livres
    integer, dimension (2) :: gdl_impostos
    
!     a1 = 0.05
!     a2 = 0.1
    a = 0.05
    e = 2.1*10**11
    l = 0.2 !m
    i = 5/(12*10**8)
    
    nele = 10
    nnodes = nele+1
    ngdl = nnodes*2
    
    cc(1,:) = (/1,0/)
    cc(2,:) = (/2,0/)
    
    fa(1,:) = (/3,0/)
    fa(2,:) = (/4,0/)
    fa(3,:) = (/5,0/)
    fa(4,:) = (/6,0/)
    fa(5,:) = (/7,0/)
    fa(6,:) = (/8,0/)
    fa(7,:) = (/9,0/)
    fa(8,:) = (/10,0/)
    fa(9,:) = (/11,0/)
    fa(10,:) = (/12,0/)
    fa(11,:) = (/13,0/)
    fa(12,:) = (/14,0/)
    fa(13,:) = (/15,0/)
    fa(14,:) = (/16,0/)
    fa(15,:) = (/17,0/)
    fa(16,:) = (/18,0/)
    fa(17,:) = (/19,0/)
    fa(18,:) = (/20,0/)
    fa(19,:) = (/21,0/)
    fa(20,:) = (/22,-200/)
    
    
    prop(1,:) = (/a, e, l, i/)
    
    
    mat_conect(1,:) = (/1,2,3,4/)
    mat_conect(2,:) = (/3,4,5,6/)
    mat_conect(3,:) = (/5,6,7,8/)
    mat_conect(4,:) = (/7,8,9,10/)
    mat_conect(5,:) = (/9,10,11,12/)
    mat_conect(6,:) = (/11,12,13,14/)
    mat_conect(7,:) = (/13,14,15,16/)
    mat_conect(8,:) = (/15,16,17,18/)
    mat_conect(9,:) = (/17,18,19,20/)
    mat_conect(10,:) = (/19,20,21,22/)
    kop = 0
    kop(1,:) =     (/12.0d0,  6.0d0*l,   -12.0d0,  6.0d0*l/)
        kop(2,:) = (/6.0d0*l, 4.0d0*l*2, -6.0d0*l, 2.0d0*l**2/)
        kop(3,:) = (/-12.0d0, -6.0d0*l,      12.0d0, -6.0d0*l/)
        kop(4,:) = (/6.0d0*l, 2.0d0*l**2,-6.0d0*l, 4.0d0*l**2/)
    
    ident = 0
    do ii=1,ngdl
        ident(ii,ii) = 1
    end do
    kglobal = 0
    kglobal1 = 0
    kele = 0
    
    coef = 109375.0d0
    do ii=1,nele
                
                kele = coef*kop
                mat_transf (1,:) = ident(mat_conect(ii,1),:)
                mat_transf (2,:) = ident(mat_conect(ii,2),:)
                mat_transf (3,:) = ident(mat_conect(ii,3),:)
                mat_transf (4,:) = ident(mat_conect(ii,4),:)
                transmat = transpose(mat_transf)
                call matmu(transmat,kele,mat_transf,kglobal1)
                kglobal = kglobal+kglobal1
            end do
       print*, kglobal
!     gdl_livres = fa(:,1)
!    
!      gdl_impostos = cc(:,1)
!     
!      kll = kglobal(gdl_livres,gdl_livres)
! !      kli = kglobal(gdl_livres,gdl_impostos)
!      
!      
!      n = 10
!      call inv(kll, invkll, n) 
!      flivres = fa(:,2)
!      u_liv = matmul(invkll,fa(:,2))
!     
!      print*, u_liv
!         
     
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
!     subroutine kop (k)
!         
!         double precision, dimension (4,4) :: k
!         
!         k(1,:) = (/12, 6*l, -12, 6*l/)
!         k(2,:) = (/0, 4*l*2, -6*l, 2*l**2/)
!         k(3,:) = (/0, 0, 12, -6*l/)
!         k(4,:) = (/0, 0, 0, 4*l**2/)
!         
!         end subroutine kop
!     
