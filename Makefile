ifeq ($(DEBUG),yes)
    OPT_FLAG = -g3
    PETSC_ARCH = arch-linux-debug
    BUILD_DIR = debug
else
    OPT_FLAG = -g0 -O3 -march=native -mtune=native
    PETSC_ARCH = arch-linux-opt
    BUILD_DIR = opt
endif

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

COMPILE = mpicxx -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas ${OPT_FLAG} \
	  -I ${MYLIB_DIR}/Eigen -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include \
	  -I ./include
LINK = mpicxx -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas ${OPT_FLAG}

output: src/main.o src/mesh.o
	${LINK} main.o mesh.o -o output ${PETSC_KSP_LIB}
src/main.o: src/main.cc
	${COMPILE} -c src/main.cc
src/mesh.o: src/mesh.cc
	${COMPILE} -c src/mesh.cc
