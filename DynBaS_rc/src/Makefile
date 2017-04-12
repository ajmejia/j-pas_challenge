.SUFFIXES:
.SUFFIXES: .f90 .o

RM       = /bin/rm -f

#FC       = gfortran -ffree-form -fcheck=all -Wunused -Wuninitialized -fopenmp
FC       = gfortran -ffree-form -fopenmp

%.o: %.f90
	$(FC) $< -c

dynbas: dynbas.f90 dynbas_pieces_module.o nnls.o brent.o locate.o
	$(FC) dynbas.f90 dynbas_pieces_module.o nnls.o  brent.o locate.o -o dynbas
	mv dynbas ../bin/.

brent.o: brent.f
	gfortran -std=legacy -ffixed-line-length-none brent.f -c

locate.o: locate.f
	gfortran -std=legacy -ffixed-line-length-none locate.f -c

nnls.o: nnls.f
	gfortran -std=legacy -ffixed-line-length-none nnls.f -c

dynbas_pieces_module.o: dynbas_pieces_module.f90 nnls.o

clean:
	$(RM) *.o *.mod $(DYNBAS)/bin/*
