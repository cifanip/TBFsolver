# Makefile
EXE = exe
POIS_OPT = FAST_MODE
MEM_OPT = MEM_SAVE
COMPILER = mpiifort


ifeq ($(POIS_OPT),MG_MODE)
OBJS = kinds.o formats.o interpolation.o initMpiGvar.o ompRoutines.o errorHandler.o allocateArrays.o mpiControl.o \
grid.o dictionary.o scalarField.o rbgs.o poissMat.o multiGrid.o vectorField.o momentumEqn.o \
poissonEqn.o time.o auxiliaryRoutines.o initialConditions.o VOF.o statistics.o pcg.o rampUpProp.o main.o
else ifeq ($(POIS_OPT),FAST_MODE)
OBJS = kinds.o formats.o interpolation.o initMpiGvar.o ompRoutines.o errorHandler.o allocateArrays.o \
mpiControl.o grid.o dictionary.o scalarField.o vectorField.o momentumEqn.o poissonEqn.o time.o \
auxiliaryRoutines.o initialConditions.o statistics.o pencilDec.o fastPoissonSolver.o VOF.o rampUpProp.o main.o
endif

MUMPS_FLAGS = -ldmumps -lmumps_common -lpord  -lmkl_blas95_lp64 \
-lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl
MKL_FLAGS = -lmkl_blas95_lp64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
-lmkl_blacs_intelmpi_lp64
FFTW_FLAGS = -lfftw3_omp -lfftw3 -lm

#######  compilation flags
#extended debug
#FLAGS = -g -check all -warn -traceback -debug extended

#debug
#FLAGS = -g -check all -traceback

#optimal 
FLAGS = -O2 -traceback

ifeq ($(POIS_OPT),MG_MODE)
$(EXE) : $(OBJS)
	$(COMPILER) -o $(EXE) -I$(SURFSARA_MUMPS_INCLUDE) -L$(SURFSARA_MUMPS_LIBS) $(OBJS) $(MUMPS_FLAGS) -qopenmp
else ifeq ($(POIS_OPT),FAST_MODE)	
$(EXE) : $(OBJS)
	$(COMPILER) $(OBJS) $(FFTW_FLAGS) $(MKL_FLAGS) -qopenmp -o $(EXE)
endif

kinds.o: kinds/kinds.f90
	$(COMPILER) -c $(FLAGS) kinds/kinds.f90

formats.o: IOobjects/formats.f90
	$(COMPILER) -c $(FLAGS) IOobjects/formats.f90

interpolation.o: kinds.o interpolation/interpolation.f90
	$(COMPILER) -D$(POIS_OPT) -c $(FLAGS) interpolation/interpolation.f90 -qopenmp

initMpiGvar.o: initMpiGvar/initMpiGvar.f90
	$(COMPILER) -c $(FLAGS) initMpiGvar/initMpiGvar.f90 -qopenmp

ompRoutines.o: ompRoutines/ompRoutines.f90
	$(COMPILER) -c $(FLAGS) ompRoutines/ompRoutines.f90 -qopenmp

errorHandler.o: initMpiGvar.o errorHandler/errorHandler.f90
	$(COMPILER) -c $(FLAGS) errorHandler/errorHandler.f90

allocateArrays.o: kinds.o errorHandler.o allocateArrays/allocateArrays.f90
	$(COMPILER) -D$(POIS_OPT) -c $(FLAGS) allocateArrays/allocateArrays.f90

dictionary.o: kinds.o errorHandler.o formats.o IOobjects/dictionary.f90
	$(COMPILER) -c $(FLAGS) IOobjects/dictionary.f90

scalarField.o: interpolation.o ompRoutines.o grid.o IOobjects/scalarField/scalarField.f90 \
IOobjects/scalarField/scalarBoundaryField_H.f90 IOobjects/scalarField/scalarBoundaryField_S.f90 
	$(COMPILER) -D$(POIS_OPT) -D$(MEM_OPT) -c $(FLAGS) IOobjects/scalarField/scalarField.f90 -qopenmp

grid.o: mpiControl.o grid/grid.f90
	$(COMPILER) -D$(POIS_OPT) -c $(FLAGS) grid/grid.f90 -qopenmp

rbgs.o: scalarField.o rbgs/rbgs.f90
	$(COMPILER) -c $(FLAGS) rbgs/rbgs.f90 -qopenmp

multiGrid.o: rbgs.o poissMat.o multiGrid/multiGrid.f90
	$(COMPILER) -c $(FLAGS) multiGrid/multiGrid.f90 -qopenmp
	
pencilDec.o: grid.o pencilDec/pencilDec.f90
	$(COMPILER) -c $(FLAGS) pencilDec/pencilDec.f90

fastPoissonSolver.o: scalarField.o initialConditions.o pencilDec.o fastPoissonSolver/fastPoissonSolver.f90
	$(COMPILER) -c $(FLAGS) -qopenmp fastPoissonSolver/fastPoissonSolver.f90

poissMat.o: scalarField.o poissMat/poissMat.f90
	$(COMPILER) -c $(FLAGS) -I$(SURFSARA_MUMPS_INCLUDE) poissMat/poissMat.f90 -qopenmp

mpiControl.o: allocateArrays.o dictionary.o mpiControl/mpiControl.f90
	$(COMPILER) -D$(POIS_OPT) -c $(FLAGS) mpiControl/mpiControl.f90

vectorField.o: scalarField.o IOobjects/vectorField/vectorField.f90
	$(COMPILER) -c $(FLAGS) IOobjects/vectorField/vectorField.f90

momentumEqn.o: time.o momentumEqn/momentumEqn.f90
	$(COMPILER) -D$(POIS_OPT) -c $(FLAGS) momentumEqn/momentumEqn.f90 -qopenmp

time.o: auxiliaryRoutines.o time/time.f90
	$(COMPILER) -c $(FLAGS) time/time.f90

ifeq ($(POIS_OPT),MG_MODE)
poissonEqn.o: pcg.o time.o poissonEqn/poissonEqn.f90
	$(COMPILER) -D$(POIS_OPT) -c $(FLAGS) poissonEqn/poissonEqn.f90 -qopenmp
else ifeq ($(POIS_OPT),FAST_MODE)	
poissonEqn.o: fastPoissonSolver.o time.o poissonEqn/poissonEqn.f90
	$(COMPILER) -D$(POIS_OPT) -c $(FLAGS) \
poissonEqn/poissonEqn.f90 -qopenmp
endif

auxiliaryRoutines.o: vectorField.o auxiliaryRoutines/auxiliaryRoutines.f90
	$(COMPILER) -c $(FLAGS) auxiliaryRoutines/auxiliaryRoutines.f90 -qopenmp

initialConditions.o: vectorField.o initialConditions/initialConditions.f90
	$(COMPILER) -c $(FLAGS) initialConditions/initialConditions.f90 -qopenmp

VOF.o: time.o VOF/VOF.f90
	$(COMPILER) -c $(FLAGS) VOF/VOF.f90 -qopenmp

statistics.o: vectorField.o statistics/statistics.f90
	$(COMPILER) -c $(FLAGS) statistics/statistics.f90 -qopenmp	

pcg.o: multiGrid.o pcg/pcg.f90
	$(COMPILER) -c $(FLAGS) pcg/pcg.f90 -qopenmp

rampUpProp.o: dictionary.o rampUpProp/rampUpProp.f90
	$(COMPILER) -c $(FLAGS) rampUpProp/rampUpProp.f90

ifeq ($(POIS_OPT),MG_MODE)
main.o: momentumEqn.o poissonEqn.o VOF.o initialConditions.o statistics.o rampUpProp.o \
createFields_H.f90 writeFields_H.f90 main.f90
	$(COMPILER) -D$(POIS_OPT) -D$(MEM_OPT) -c $(FLAGS) main.f90
else ifeq ($(POIS_OPT),FAST_MODE)
main.o: momentumEqn.o poissonEqn.o VOF.o initialConditions.o statistics.o rampUpProp.o \
createFields_H.f90 writeFields_H.f90 main.f90
	$(COMPILER) -D$(POIS_OPT) -D$(MEM_OPT) -c $(FLAGS) main.f90
endif

clean : 
	rm *.o *.mod $(EXE)
