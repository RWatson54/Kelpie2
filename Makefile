# Set up the project directory
PROJECT_DIR = /home/raw54/Source/Kelpie_2p0

# Where to hide things
BIN_DIR = $(PROJECT_DIR)/bin
SRC_DIR = $(PROJECT_DIR)/src
LIB_DIR = $(PROJECT_DIR)/lib
INC_DIR = $(PROJECT_DIR)/include

# Set up the compiler
CC=gcc
FC=gfortran

# Set the optimisation and error checking options
OPTIMISE = -O3
DEBUGGING = -g -Wall -Wno-unused-dummy-argument -fcheck=all -pedantic -fbacktrace -Wextra
INCLUDE = -I$(INC_DIR)

# Bundle these together
FCFLAGS= $(OPTIMISE) $(DEBUGGING) $(INCLUDE)
CCFLAGS= $(OPTIMISE) $(DEBUGGING) $(INCLUDE)

# Set up any libraries needed
LFLAGS= -L$(LIB_DIR) 

# Set up the list of source code files
EXEC_FILES = $(SRC_DIR)/precision.f90 \
             $(SRC_DIR)/hello.f90 \
             $(SRC_DIR)/parameters.f90 \
             $(SRC_DIR)/input.f90 \
             $(SRC_DIR)/readalloc.f90 \
             $(SRC_DIR)/main.f90 

MESH_FILES = $(SRC_DIR)/precision.f90 \
             $(SRC_DIR)/mesh.f90 


# Set up the list of object files from the source code files
OBJ_FILES1 = $(EXEC_FILES:%.f90=%.o)
OBJ_FILES2 = $(MESH_FILES:%.f90=%.o)

# Set up the hooks
all: clean meshinit.exe kelpie.exe

kelpie.exe: $(OBJ_FILES1)
	@echo ""
	@echo "Building Kelpie executable"
	$(FC) $(FCFLAGS) $(OBJ_FILES1) -o $@ $(LFLAGS)
	$(MKDIR) $(BIN_DIR)
	$(MV) $@ $(BIN_DIR)/$@

meshinit.exe: $(OBJ_FILES2)
	@echo ""
	@echo "Building Mesh and initialisation executable"
	$(FC) $(FCFLAGS) $(OBJ_FILES2) -o $@ $(LFLAGS)
	$(MKDIR) $(BIN_DIR)
	$(MV) $@ $(BIN_DIR)/$@

%.o: %.F90
	@echo ""
	@echo "Compiling F90 file " $*.F90
	$(FC) -c $(FCFLAGS) -J$(SRC_DIR) -o $@ -c $<

%.o: %.f90
	@echo ""
	@echo "Compiling f90 file " $*.f90
	$(FC) -c $(FCFLAGS) -J$(SRC_DIR) -o $@ -c $<

%.o: %.f
	@echo ""
	@echo "Compiling f77 file " $*.f
	$(FC) -c $(FCFLAGS) -J$(SRC_DIR) -o $@ -c $<

%.o: %.c
	@echo ""
	@echo "Compiling C file " $*.c
	$(CC) -c $(CCFLAGS) -o $@ -c $<

clean:
	$(RM) *.o *.mod *~ *# fort.*
	${RM} $(SRC_DIR)/*.o $(SRC_DIR)/*~ $(SRC_DIR)/*.mod $(SRC_DIR)/*#
	${RM} ${INC_DIR}/*.mod
	${RM} ${BIN_DIR}/*.exe

#
# Utilities
#

CD     = cd
MV     = mv
RM     = \rm -f
CP     = \cp -f
RMDIR  = \rm -rf
MKDIR  = \mkdir -p
AR     = \ar -crs
