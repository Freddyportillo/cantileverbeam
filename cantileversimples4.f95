program cantilever5ele

implicit none

    
    
    integer, dimension (2,2) :: cc
    double precision, dimension (4,12) ::  mat_transf
    double precision, dimension (1,4) :: prop
    double precision, dimension (12,4) :: transmat
    double precision :: a, e, l, i
    integer :: ii, ngdl, nele, nnodes
    double precision, dimension (4,4) :: kop, kele
    integer, dimension (10,2) :: fa
    double precision, dimension (10,10) :: invkll, kll
    integer :: n
    double precision :: coef
    double precision, dimension (12,12) :: kglobal, kglobal1
    double precision, dimension (12) :: des_global
    integer, dimension (12,12) :: ident
    integer, dimension (5,4) ::  mat_conect
    double precision, dimension (10) :: u_liv
    double precision, dimension (10) :: flivres
    integer, dimension (10) :: gdl_livres
    integer, dimension (2) :: gdl_impostos
    double precision, dimension (12) :: desloc_trans
    

    a = 0.05
    e = 2.1*10**11
    l = 0.4 !m
    i = 4.16666/10**9
    
    nele = 5
    nnodes = nele+1
    ngdl = nnodes*2
    
    cc(1,:) = (/1,0/)
    cc(2,:) = (/2,0/)
    
    do ii = 1,8
        fa(ii,:) = (/ii+2,0/)
    end do

    fa(9,:) = (/11,-30/) !forca aplicada no extremo livre da viga
    fa(10,:) = (/12, 0/)
    
    prop(1,:) = (/a, e, l, i/)
    
    
    mat_conect(1,:) = (/1,2,3,4/)
    mat_conect(2,:) = (/3,4,5,6/)
    mat_conect(3,:) = (/5,6,7,8/)
    mat_conect(4,:) = (/7,8,9,10/)
    mat_conect(5,:) = (/9,10,11,12/)
    
    kop = 0
    kop(1,:) =     (/12.0d0,  6.0d0*l,   -12.0d0,  6.0d0*l/)
        kop(2,:) = (/6.0d0*l, 4.0d0*l**2, -6.0d0*l, 2.0d0*l**2/)
        kop(3,:) = (/-12.0d0, -6.0d0*l,      12.0d0, -6.0d0*l/)
        kop(4,:) = (/6.0d0*l, 2.0d0*l**2,-6.0d0*l, 4.0d0*l**2/)
    
    ident = 0
    do ii=1,ngdl
        ident(ii,ii) = 1
    end do
    kglobal = 0
    kglobal1 = 0
    kele = 0
    
    coef = e*i/(l**3) 
    do ii=1,nele
                
                kele = coef*kop
                mat_transf (1,:) = ident(mat_conect(ii,1),:)
                mat_transf (2,:) = ident(mat_conect(ii,2),:)
                mat_transf (3,:) = ident(mat_conect(ii,3),:)
                mat_transf (4,:) = ident(mat_conect(ii,4),:)
                transmat = transpose(mat_transf)
                call matmu(transmat,kele,mat_transf,kglobal1)
                kglobal = kglobal+kglobal1
!                 

                
            end do
       
      
       
     gdl_livres = fa(:,1)
    
    gdl_impostos = cc(:,1)
  
      kll = kglobal(gdl_livres,gdl_livres)
!       kli = kglobal(gdl_livres,gdl_impostos)
     
      
      n = 10
      call inv(kll, invkll, n) 
      flivres = fa(:,2)
      
      u_liv = matmul(invkll,fa(:,2))
      
        des_global= 0.0d0
        desloc_trans = 0.0d0
        
         des_global(cc(:,1))= cc(:,2)
         des_global(fa(:,1)) = u_liv
        print*,'Deslocamentos (em metros)'
         desloc_trans(:)=des_global(1:2:nele)
         
      print*,
      
      print*, des_global
         
     
    end program
    
    subroutine matmu (A,B,C,D)
      !programa para multiplicacao de 3 matrizes
  
  double precision, dimension (12,4) :: A, mul
  double precision, dimension (4,4) :: B
  double precision, dimension (4,12) :: C
  double precision, dimension (12,12) :: D
  
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
