gfortran -c pacotes/lapack_tools.f pacotes/lapack.f pacotes/lapack_parcer.f08 cantileversimples4.f95
gfortran lapack_tools.o lapack.o lapack_parcer.o cantileversimples4.o -o simples4
