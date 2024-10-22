FC=gfortran
FFLAGS=-O3 
LIBS=-ldftbplus -fopenmp -lblas -llapack
INCLUDE_DIR=/home/bustamac/Programs_and_libraries/dftbplus_dev/_build/_install/include/dftbplus/modfiles
LIB_DIR=/home/bustamac/Programs_and_libraries/dftbplus_dev/_build/_install/lib

SRC:=constants_mod.F90
SRC+=grid_mod.F90
SRC+=classical_medium_mod.F90
SRC+=io_mod.F90
SRC+=testhelpers.F90
SRC+=q_medium_mod.F90
SRC+=td_propagator_mod.F90
MAIN=Mxll_Mxim.F90
MOD=${SRC:.F90=.mod}
OBJ=${SRC:.F90=.o}
OBJ+=Mxll_Mxim.o
EXC=program.e

%.o: %.F90 #wildcard rule, creation of *.o depends on *.f90
	$(FC) $(FFLAGS) -o $@ -c $< -I$(INCLUDE_DIR)

$(EXC): $(OBJ)
			$(FC) -o $@ $^ ${FFLAGS} -L$(LIB_DIR) $(LIBS)

.PHONY: clean
clean:
	      rm $(OBJ) $(MOD) $(EXC)
