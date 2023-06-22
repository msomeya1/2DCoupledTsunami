FC = ifort

a.out: sub_fdm4.o sub_program.o Lagrangian.f90
	$(FC) -xHOST -parallel Lagrangian.f90 sub_fdm4.o sub_program.o

# a.out: sub_fdm4.o sub_program.o MaedaFurumura.f90
# 	$(FC) -xHOST -parallel MaedaFurumura.f90 sub_fdm4.o sub_program.o

# a.out: sub_fdm4.o sub_program.o Eulerian.f90
# 	$(FC) -xHOST -parallel Eulerian.f90 sub_fdm4.o sub_program.o


sub_fdm4.o: sub_fdm4.f90
	$(FC) -c sub_fdm4.f90

sub_program.o: sub_program.f90
	$(FC) -c sub_program.f90