# Makefile
EXE = TBFsolver
POIS_OPT = FAST_MODE
MEM_OPT = MEM_SAVE
#select the compiler
COMPILER = 

#if POIS_OPT=MG_MODE add paths to MUMPS lib
MUMPS_INCL = 
MUMPS_LIB_DIR = 

MUMPS_FLAGS = -ldmumps -lmumps_common -lpord -lmetis -lesmumps \
-lscotch -lscotcherr
FFTW_FLAGS = -lfftw3_omp -lfftw3 -lm
LAPACK_FLAGS = -lopenblas -lscalapack

ifeq ($(POIS_OPT),MG_MODE)
LIBS = $(MUMPS_FLAGS) $(LAPACK_FLAGS)
LIB_DIRS = $(MUMPS_LIB_DIR)
else ifeq ($(POIS_OPT),FAST_MODE)
LIBS = $(FFTW_FLAGS) $(LAPACK_FLAGS)
#path to above LIBS may be necessary. if so, add below
LIB_DIRS = 
endif

#modify the flags according to the chosen compiler
OMP_FLAG = -qopenmp

STACK_OPT =
#STACK_OPT = -heap-arrays

#debug
#FLAGS = -g -check all -traceback #-fpe0

#optimal 
FLAGS = -O2 -traceback

ifeq ($(POIS_OPT),MG_MODE)
OBJS = solverTypes.o kinds.o formats.o interpolation.o initMpiGvar.o ompRoutines.o errorHandler.o allocateArrays.o \
mpiControl.o grid.o parFile.o field.o rbgs.o poissMat.o multiGrid.o vfield.o momentumEqn.o vofBlocks.o \
poissonEqn.o time.o auxiliaryRoutines.o initialConditions.o VOF.o statistics.o pcg.o rampUpProp.o main.o
else ifeq ($(POIS_OPT),FAST_MODE)
OBJS = solverTypes.o kinds.o formats.o interpolation.o initMpiGvar.o ompRoutines.o errorHandler.o allocateArrays.o \
mpiControl.o grid.o parFile.o field.o vfield.o momentumEqn.o poissonEqn.o time.o vofBlocks.o \
auxiliaryRoutines.o initialConditions.o statistics.o pencilDec.o fastPoissonSolver.o VOF.o rampUpProp.o main.o
endif

ifeq ($(POIS_OPT),MG_MODE)
$(EXE) : $(OBJS)
	$(COMPILER) -o $(EXE) -I$(MUMPS_INCL) -L$(LIB_DIRS) $(OBJS) $(LIBS) $(OMP_FLAG)
else ifeq ($(POIS_OPT),FAST_MODE)
$(EXE) : $(OBJS)
	$(COMPILER) -o $(EXE) $(OBJS) $(LIBS) $(OMP_FLAG)
endif

solverTypes.o: solverTypes/solverTypes.f90
	$(COMPILER) -c $(FLAGS) solverTypes/solverTypes.f90

kinds.o: kinds/kinds.f90
	$(COMPILER) -c $(FLAGS) kinds/kinds.f90

formats.o: fields/formats.f90
	$(COMPILER) -c $(FLAGS) fields/formats.f90

interpolation.o: kinds.o interpolation/interpolation.f90
	$(COMPILER) -D$(POIS_OPT) -c -cpp $(FLAGS) interpolation/interpolation.f90 $(OMP_FLAG)

initMpiGvar.o: initMpiGvar/initMpiGvar.f90
	$(COMPILER) -c $(FLAGS) initMpiGvar/initMpiGvar.f90 $(OMP_FLAG)

ompRoutines.o: ompRoutines/ompRoutines.f90
	$(COMPILER) -c $(FLAGS) ompRoutines/ompRoutines.f90 $(OMP_FLAG)

errorHandler.o: initMpiGvar.o errorHandler/errorHandler.f90
	$(COMPILER) -c $(FLAGS) errorHandler/errorHandler.f90

allocateArrays.o: kinds.o errorHandler.o allocateArrays/allocateArrays.f90
	$(COMPILER) -D$(POIS_OPT) -c -cpp $(FLAGS) allocateArrays/allocateArrays.f90

parFile.o: kinds.o errorHandler.o formats.o fields/parFile.f90
	$(COMPILER) -c $(FLAGS) fields/parFile.f90

field.o: interpolation.o ompRoutines.o grid.o fields/field/field.f90 \
fields/field/boundaryField_H.f90 fields/field/boundaryField_S.f90 
	$(COMPILER) -D$(POIS_OPT) -D$(MEM_OPT) -c -cpp $(FLAGS) fields/field/field.f90 $(OMP_FLAG)

grid.o: mpiControl.o grid/grid.f90
	$(COMPILER) -D$(POIS_OPT) -c -cpp $(FLAGS) grid/grid.f90 $(OMP_FLAG)

rbgs.o: field.o rbgs/rbgs.f90
	$(COMPILER) -c $(FLAGS) rbgs/rbgs.f90 $(OMP_FLAG)

multiGrid.o: rbgs.o poissMat.o multiGrid/multiGrid.f90
	$(COMPILER) -c $(FLAGS) multiGrid/multiGrid.f90 $(OMP_FLAG)
	
pencilDec.o: grid.o pencilDec/pencilDec.f90
	$(COMPILER) -c $(FLAGS) pencilDec/pencilDec.f90

fastPoissonSolver.o: field.o initialConditions.o pencilDec.o fastPoissonSolver/fastPoissonSolver.f90
	$(COMPILER) -c $(FLAGS) $(OMP_FLAG) fastPoissonSolver/fastPoissonSolver.f90

poissMat.o: field.o poissMat/poissMat.f90
	$(COMPILER) -c $(FLAGS) -I$(MUMPS_INCL) poissMat/poissMat.f90 $(OMP_FLAG)

mpiControl.o: allocateArrays.o parFile.o mpiControl/mpiControl.f90
	$(COMPILER) -D$(POIS_OPT) -c -cpp $(FLAGS) mpiControl/mpiControl.f90

vfield.o: field.o fields/vfield/vfield.f90
	$(COMPILER) -c $(FLAGS) fields/vfield/vfield.f90

momentumEqn.o: time.o momentumEqn/momentumEqn.f90
	$(COMPILER) -D$(POIS_OPT) -c -cpp $(FLAGS) momentumEqn/momentumEqn.f90 $(OMP_FLAG)

time.o: solverTypes.o auxiliaryRoutines.o initialConditions.o time/time.f90
	$(COMPILER) -c $(FLAGS) time/time.f90

ifeq ($(POIS_OPT),MG_MODE)
poissonEqn.o: pcg.o time.o poissonEqn/poissonEqn.f90
	$(COMPILER) -D$(POIS_OPT) -c -cpp $(FLAGS) poissonEqn/poissonEqn.f90 $(OMP_FLAG)
else ifeq ($(POIS_OPT),FAST_MODE)	
poissonEqn.o: fastPoissonSolver.o time.o poissonEqn/poissonEqn.f90
	$(COMPILER) -D$(POIS_OPT) -c -cpp $(FLAGS) \
poissonEqn/poissonEqn.f90 $(OMP_FLAG)
endif

auxiliaryRoutines.o: vfield.o auxiliaryRoutines/auxiliaryRoutines.f90
	$(COMPILER) -c $(FLAGS) auxiliaryRoutines/auxiliaryRoutines.f90 $(OMP_FLAG)

initialConditions.o: vfield.o initialConditions/initialConditions.f90
	$(COMPILER) -c $(FLAGS) initialConditions/initialConditions.f90 $(OMP_FLAG)

vofBlocks.o: initialConditions.o VOF/vofBlocks.f90
	$(COMPILER) -c $(FLAGS) $(STACK_OPT) VOF/vofBlocks.f90 $(OMP_FLAG)

VOF.o: vofBlocks.o time.o VOF/VOF.f90
	$(COMPILER) -c $(FLAGS) VOF/VOF.f90 $(OMP_FLAG)

statistics.o: time.o statistics/statistics.f90
	$(COMPILER) -c $(FLAGS) statistics/statistics.f90 $(OMP_FLAG)	

pcg.o: multiGrid.o pcg/pcg.f90
	$(COMPILER) -c $(FLAGS) pcg/pcg.f90 $(OMP_FLAG)

rampUpProp.o: parFile.o rampUpProp/rampUpProp.f90
	$(COMPILER) -c $(FLAGS) rampUpProp/rampUpProp.f90

ifeq ($(POIS_OPT),MG_MODE)
main.o: momentumEqn.o poissonEqn.o VOF.o initialConditions.o statistics.o rampUpProp.o \
createFields_H.f90 writeFields_H.f90 main.f90
	$(COMPILER) -D$(POIS_OPT) -D$(MEM_OPT) -c -cpp $(FLAGS) main.f90
else ifeq ($(POIS_OPT),FAST_MODE)
main.o: momentumEqn.o poissonEqn.o VOF.o initialConditions.o statistics.o rampUpProp.o \
createFields_H.f90 writeFields_H.f90 main.f90
	$(COMPILER) -D$(POIS_OPT) -D$(MEM_OPT) -c -cpp $(FLAGS) main.f90
endif

clean : 
	rm -f *.o *.mod $(EXE)
