program cantilever5ele
  use lapack_parcer

implicit none

    
    
    integer, dimension (2,2) :: cc
    double precision, dimension (4,12) ::  mat_transf
    double precision, dimension (1,4) :: prop
    double precision, dimension (12,4) :: transmat
    double precision :: a, e, l, i, ymax, x, ro
    integer :: ii, ngdl, nele, nnodes
    double precision, dimension (4,4) :: kop, kele, kdin, mele
    integer, dimension (10,2) :: fa
    double precision, dimension (10,10) :: invkll
    double precision, allocatable, dimension (:,:) :: kll, mll, eigvet
    double precision, allocatable, dimension (:) :: EigVal
    double precision, dimension (10,2) :: kli
    double precision, dimension (2,10) :: transkli
    double precision, dimension (2,2) :: kii
    integer :: n, nmodos
    double precision :: coef, coef_massa
    double precision, dimension (12,12) :: kglobal, kglobal1, mglobal, mglobal1
    double precision, dimension (12) :: des_global
    integer, dimension (12,12) :: ident
    integer, dimension (5,4) ::  mat_conect
    double precision, dimension (10) :: u_liv, freq, flivres
    double precision, dimension (5) :: freqs
    integer, dimension (10) :: gdl_livres
    integer, dimension (2) :: gdl_impostos
    double precision, dimension (2,1) :: u_imp
    double precision, dimension (2) :: f_imp
    double precision, dimension (6) :: desloc_trans, rot
    double precision, dimension (4,1) :: u_e
    double precision, dimension (1,4) :: psi2_0, psi2_l
    double precision, dimension (5,2) :: mat_sigma
    double precision, dimension (1,1) :: psi_ue1, psi_ue2
    double precision :: pi = 3.14159265359
   

    

  a = 0.05
  e = 2.1*10**11
  l = 0.4 !m
  i = 4.16666/10**9
  ro = 3.9    !kg/m
  ymax = 0.005
    
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
  end do
  n = 10    
  gdl_livres = fa(:,1)
    
  gdl_impostos = cc(:,1)
  allocate (kll(n,n))
  kll = kglobal(gdl_livres,gdl_livres)
  kli = kglobal(gdl_livres, gdl_impostos)
  kii = kglobal(gdl_impostos,gdl_impostos)
  call inv(kll, invkll, n) 
  flivres = fa(:,2)
      
  u_liv = matmul(invkll,fa(:,2))
! u_imp = cc(:,2)
  des_global= 0.0d0
  desloc_trans = 0.0d0
  des_global(cc(:,1))= cc(:,2)
  des_global(fa(:,1)) = u_liv
    print*, 'Deslocamentos (em metros)'
  desloc_trans=des_global(1:nele:2)
    print*, desloc_trans
    print*, '----------------------------'
    print*, 'os valores da rotação da seção transversal é (rad):'
  rot = des_global(2:nele+1:2)
    print*, rot
    print*, '----------------------------'
  transkli = transpose(kli)
  f_imp = (matmul(transkli,u_liv)) 
    print*, 'as reacoes no apoio são (N):'
    print*, f_imp
         
!     CÁLCULO DAS TENSÕES NORMAIS DE FLEXÃO MÁXIMAS NOS NÓS
        
  do ii=1,nele
    u_e(:,1) = des_global(mat_conect(ii, :))
    x = 0.0d0
    psi2_0(1,:) = (/-6/l**2+12*x/l**3, -4/l+6*x/l**2, 6/l**2-12*x/l**3, -2/l+6*x/l**2/)
    x = l
    psi2_l(1,:) = (/-6/l**2+12*x/l**3, -4/l+6*x/l**2, 6/l**2-12*x/l**3, -2/l+6*x/l**2/)
    psi_ue1 = matmul(psi2_0,u_e)
    mat_sigma(ii,1) = -E*psi_ue1(1,1)*ymax
    psi_ue2 = matmul(psi2_l,u_e)
    mat_sigma(ii,2) = -E*psi_ue2(1,1)*ymax
  end do
    
        
        !ANALISE DINAMICO DA VIGA
        
kdin(1,:) = (/156.0d0, 22.0d0*l, 54.0d0 ,-13.0d0*l/)
kdin(2,:) = (/22.0d0*l, 4.0d0*l**2, 13.0d0*l, -3.0d0*l**2/)
kdin(3,:) = (/54.0d0, -13.0d0*l, 156.0d0, -22.0d0*l/)
kdin(4,:) =(/-13.0d0*l, -3.0d0*l**2, -22.0d0*l, 4.0d0*l**2/)
        
mglobal = 0
mglobal1 = 0
        
coef_massa = ro*l/420
do ii = 1,nele
        
  mele = coef_massa*kdin
  mat_transf (1,:) = ident(mat_conect(ii,1),:)
  mat_transf (2,:) = ident(mat_conect(ii,2),:)
  mat_transf (3,:) = ident(mat_conect(ii,3),:)
  mat_transf (4,:) = ident(mat_conect(ii,4),:)
  transmat = transpose(mat_transf)
  call matmu(transmat,mele,mat_transf,mglobal1)
  mglobal = mglobal + mglobal1
    end do
  allocate (mll(n,n), eigval(n), eigvet(n,n))
       
  mll = mglobal(gdl_livres,gdl_livres)
  kll = kglobal(gdl_livres,gdl_livres)     
    nmodos = 5
  call lapack_eig (n,kll,mll,EigVal,EigVet)
  
  
  
    print*, '-----------------------'
    print*, 'As frequencias naturais sao (Hz): '
     freq = sqrt(eigval)*pi/2  
     
do ii = 1, nmodos
    freqs(ii) = minval(freq)/10
    freq(minloc(freq)) = 999999.0d0
end do
     !freqs(:) = freq(1:nmodos)
    print*, freqs

    
    deallocate (kll,mll,eigval,eigvet) 

    
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
