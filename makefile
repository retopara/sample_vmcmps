OPTS = -O3 -w 
#OPTS = -O3 -msse2 -static -s -pipe -fprefetch-loop-arrays -frename-registers -ftracer 

#MYLQCC
CXX = /usr/bin/mpicxx

#INSPUR
#CXX = /opt/intel/impi/openmpi/bin/mpicxx
# LAPACK_DIR = /opt/intel/Compiler/11.1/072/mkl/lib/em64t/
# BLAS_DIR = /opt/intel/Compiler/11.1/072/mkl/lib/em64t/
#F2C_DIR = /opt/intel/Compiler/11.1/072/mkl/lib/em64t/

#Lenovo_GPU
#CXX = /opt/mvapich-1.2/bin/mpicxx
#LAPACK_DIR = /usr/lib/
#BLAS_DIR = /usr/lib/
#F2C_DIR = /share/data2/lib/

#Dirac
# CXX = /opt/openmpi/bin/mpicxx
#LAPACK_DIR = /share/data1/lib/
#BLAS_DIR = /share/data1/lib/
#F2C_DIR = /share/data2/lib/

#Einstein
#LAPACK_DIR = /share/apps/lib
#BLAS_DIR = /share/apps/lib
#F2C_DIR = /share/apps/lib

INCLUDE = -I ./include \
					-I ./headers

#LIBS = -L${LAPACK_DIR} -lmkl_lapack  \
	-L${BLAS_DIR} -lmkl_blas95_lp64 \

OBJ = .adjusttemper.o \
.target.o \
.csa.o \
.esite.o \
.exchtemper.o \
.normalize.o \
.main.o \
.random.o \
.remc.o \
.toolbox.o \
.xgenerator.o \
.fileio.o \
.flip.o \
.screenio.o \
.share.o \

#compile
vmps.exe : ${OBJ}
	${CXX} ${OPTS} ${OBJ} ${INCLUDE} -o vmps.exe

#connect
.%.o :		%.cpp
	${CXX} ${OPTS} ${INCLUDE} -c $< -o $@

#clean
clean:
	rm .*.o vmps.exe
