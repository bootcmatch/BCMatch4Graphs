# Eventually it will be generated with a configure script.
# 
.SUFFIXES: .f90 .F90 .f .F .c .o


#
# Optional packages
#

# SuperLU: sparse factorization
#  Versions 4 and 5 have a slightly different interface, need
#  to define which one. 
#
#SUPERLU_DIR=/home/users/pasqua/NUMERICAL/LIB/SLU50
#SUPERLU_INCDIR=$(SUPERLU_DIR)/include
#SUPERLU_LIBDIR=$(SUPERLU_DIR)/lib 
#SUPERLU_LIBS=-lsuperlu -L$(SUPERLU_LIBDIR)
#SUPERLU_FLAGS=-DHAVE_SUPERLU -DSLU_VERSION_5 -I$(SUPERLU_INCDIR)

# HSL library (for MC64)
#HSL_DIR=/home/users/pasqua/BootCMATCH/hsl
#HSL_INCDIR=$(HSL_DIR)/include
#HSL_LIBDIR=$(HSL_DIR)/lib
#HSL_LIBS=-lhsl_mc64 -L$(HSL_LIBDIR)
#HSL_FLAGS= -DHAVE_HSL -I$(HSL_INCDIR)

# SPRAL package for auction algorithm
#SPRAL_DIR=/home/users/pasqua/BootCMATCH/spral
#SPRAL_INCDIR=$(SPRAL_DIR)/include
#SPRAL_LIBDIR=$(SPRAL_DIR)/lib
#SPRAL_LIBS=-lspral -L$(SPRAL_LIBDIR)
#SPRAL_FLAGS=-DHAVE_SPRAL -I$(SPRAL_INCDIR)


#
# LAPACK and BLAS are mandatory for the third party libraries
#
#BLAS_LIBS=-lcblas
#
# Compilers & such
#
FC=gfortran 
CC=gcc
FCOPT=-O3  
#FCOPT=-g  
CCOPT=-O3  
#CCOPT=-g  
#
# Linking is done with CC by default; if using SPRAL and HSL you
#  may need to link in the Fortran runtime. 
#
LINKER=$(CC)
#FORTRAN_LIBS= -lgfortran


#
# Putting it all together
#
LIBS= $(SUPERLU_LIBS) $(HSL_LIBS) $(SPRAL_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(FORTRAN_LIBS)
CFLAGS=$(SUPERLU_FLAGS) $(HSL_FLAGS) $(SPRAL_FLAGS)


#
# Compile rules
#
.c.o:
	$(CC) $(CCOPT) $(CFLAGS) $(CINCLUDES) -c $< -o $@
.f.o:
	$(FC) $(FCOPT) $(FINCLUDES)  -c $< -o $@
.f90.o:
	$(FC) $(FCOPT) $(FINCLUDES) -c $< -o $@
.F.o:
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $< -o $@
.F90.o:
	$(FC) $(FCOPT) $(FINCLUDES) $(FDEFINES) -c $< -o $@

