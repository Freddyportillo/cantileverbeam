program cantilever_ex

    implicit none
    
        !Foi gerada a submatriz de rigidez kll
        real*8, dimension (3,3) :: prop, kll, invkll
        real, dimension (3,1) ::  kli, kii
        real, dimension (1,3) :: transkli
        real, dimension (1) :: f_imp
        real, dimension (1,2) :: cc
        real, dimension (3,2) :: fa
        integer, dimension (3,2) :: mat_conect
        integer :: ii, nele, nnodes, ngdl,n
        real, dimension (3) :: coef, u_liv
        integer, dimension (3) :: gdl_livres
        real, dimension (2,2) :: kop, kele
        integer, dimension (4,4) :: ident
        real, dimension (2,4) :: mat_transf
        real, dimension (4,2) :: transmat
        real, dimension (4,4) :: kglobal, kglobal1
        real :: a,e
        integer,dimension (1) :: gdl_impostos
            nele = 3
            nnodes = 4
            ngdl = 4
            a = 25.3/10**4
            e = 9*10**10
            prop(1,:) = (/a, e, 1.0/)
            prop(2,:) = (/a, e, 0.5/)
            prop(3,:) = (/a, e, 3.5/)
            
            cc (1,:) = (/1,0/)
            
            fa (1,:) = (/2, -1000/)
            fa (2,:) = (/3, 1000/)
            fa (3,:) = (/4, -2000/)
            
            mat_conect (1,:) = (/1,2/)
            mat_conect (2,:) = (/2,3/)
            mat_conect (3,:) = (/3,4/)
            
        
            !construção de matriz de rigidez global
            
            kop (1,:) = (/1,-1/)
            kop (2,:) = (/-1,1/)
            ident = 0
            kglobal = 0
            kglobal1 = 0
            kele = 0
            
            do ii = 1,ngdl
                ident(ii,ii) = 1
            end do
            do ii=1,nele
                coef = prop(ii,1)*prop(ii,2)/prop(ii,3)
                kele = coef(ii)*kop
                mat_transf (1,:) = ident(mat_conect(ii,1),:)
                mat_transf (2,:) = ident(mat_conect(ii,2),:)
                transmat = transpose(mat_transf)
                call matmu(transmat,kele,mat_transf,kglobal1)
                kglobal = kglobal+kglobal1
            end do
            
            print*, kglobal
            
            gdl_livres = fa(:,1)
            gdl_impostos = cc(:,1)
            
            kll = kglobal(gdl_livres,gdl_livres)
            kli = kglobal(gdl_livres,gdl_impostos)
            print*, '----------'
            print*, kll
            print*, '----------'
            print*, kli
            !kii = kglobal(gdl_impostos,gdl_impostos)
            
            !deslocamentos no gdl gdl_livres
            n = 3.0
            call inv(kll, invkll, n)
            
            u_liv = matmul(invkll,fa(:,2))*1000
            print*, 'os valosres dos deloscamentos no gdl livres sao:'
            print*,u_liv
            print*, 'os valores das reacoes sao:'
            transkli = transpose(kli)
            u_liv = u_liv/1000 !em m
            f_imp = matmul(transkli,u_liv)/1000 !em kN
            print*, f_imp
end program
                    
        subroutine matmu (A,B,C,D)
      !programa para multiplicacao de 3 matrizes
  
  real, dimension (4,2) :: A, mul
  real, dimension (2,2) :: B
  real, dimension (2,4) :: C
  real, dimension (4,4) :: D
  
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
             
                
