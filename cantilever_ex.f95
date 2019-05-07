program cantilever_ex

    implicit none
    
    
        real, dimension (3,3) :: prop, kll
        real, dimension (3,1) ::  kli, kii
        real, dimension (1,2) :: cc
        real, dimension (3,2) :: fa
        integer, dimension (3,2) :: mat_conect
        integer :: ii, nele, nnodes, ngdl
        real, dimension (3) :: coef
        integer, dimension (3) :: gdl_livres
        real, dimension (2,2) :: kop, kele
        integer, dimension (4,4) :: ident
        real, dimension (2,4) :: mat_transf
        real, dimension (4,2) :: transmat
        real, dimension (4,4) :: kglobal, kglobal1
        real :: a,e
        integer,dimension (1,1) :: gdl_impostos
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
            kii = kglobal(gdl_impostos,gdl_impostos)

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
             
                
