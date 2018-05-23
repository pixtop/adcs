include Make.inc

# Machines salle TP N7
MATLAB_BIN=/applications/matlab/bin
# MATLAB_BIN=/opt/MATLAB/R2017a/bin

all: mexfile

matlabsetup:
	${MATLAB_BIN}/mex -setup
	${MATLAB_BIN}/mex -setup FORTRAN
	cp -f mex_FORTRAN_glnxa64.xml mex_C_glnxa64.xml $(HOME)/.matlab/R2017a/

mexfile:
	${MATLAB_BIN}/mex -Dmex -compatibleArrayDims -c subspace_iter_ev.F90 tools.f90 
	${MATLAB_BIN}/mex -Dmex -compatibleArrayDims mex_subspace_iter_ev.c tools.o subspace_iter_ev.o $(PLIBS) -lgfortran

clean:
	(rm -f *.o *.mod main *.mexa64)

%.o: %.f90
	$(FC) -c $<  

%.o: %.F90
	$(FC) -c $<  

%.o: %.f
	$(FC) -c $<  

