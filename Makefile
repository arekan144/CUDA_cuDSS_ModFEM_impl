COMP = g++
NCOMP = nvcc
ARH = ar
LINK = g++
NLINK = nvcc

ARCH = -arch=sm_60
#OPS = -O3
NOPS = --fdevice-time-trace $(shell od -An -N2 -i /dev/random | tr -d ' ')trace 
OPS = -g -D _DEBUG=1

NAME = "defaultnametoreplace"


LIB_EIGEN = -I./Eigen5 -I./Eigen5/Eigen
CTK_DIR=/usr/local/cuda-12.9/
CUDSS_DIR=/usr/lib/x86_64-linux-gnu/libcudss/12
CUDSS_INCLUDE=/usr/include/libcudss/12

NINCL = -I$(CTK_DIR)include -I$(CUDSS_INCLUDE)
NLIBS = -L$(CUDSS_DIR) -L$(CTK_DIR)lib64 -lcudart -lcublas -lcublasLt -lculibos -lcudss 
#NLIBS = -lcudart -lcublas -lcublasLt -lculibos -lcudss 

SNLIBINCL = -I${CUDSS_INCLUDE} \
-Xlinker=${CUDSS_DIR}/libcudss_static.a \
-Xlinker=${CTK_DIR}/lib64/libcublas_static.a \
-Xlinker=${CTK_DIR}/lib64/libcublasLt_static.a \
-Xlinker=${CTK_DIR}/lib64/libculibos.a

all: libEigenWrapper.a libSparseStructures.a libCUDSSWrapper.a main.exe

clean:
	rm *.o

remove:
	rm *.a
	rm main.exe

rebuild: clean remove all


libEigenWrapper.a: eigen_wrapper.cpp eigen_wrapper.h
	$(COMP) -c eigen_wrapper.cpp $(LIB_EIGEN) -o eigen_wrapper.o $(OPS)
	$(ARH) crs $@ eigen_wrapper.o

libSparseStructures.a: ./matrixes/CSR_class.cpp ./matrixes/SparseStructures.h
	$(COMP) $(OPS)  -o SparseStructures.o -c ./matrixes/CSR_class.cpp
	$(ARH) crs $@ SparseStructures.o

libCUDSSWrapper.a: cudss_wrapper.cu cudss_wrapper.h
	$(NCOMP) -Xcompiler $(OPS) -o cudss_wrapper.o -c cudss_wrapper.cu $(NOPS)
	$(NLINK) -Xcompiler $(OPS) -dlink -o cudss_wrapper.dlink.o cudss_wrapper.o $(NOPS) 
	$(NLINK) -Xcompiler $(OPS) -lib -o libCUDSSWrapper.a cudss_wrapper.o cudss_wrapper.dlink.o $(NOPS) 

main.o: main.cpp
	$(COMP) $(OPS) -o $@ -c $<

main.exe: main.o libEigenWrapper.a libSparseStructures.a libCUDSSWrapper.a
	$(NLINK) $(SNLIBINCL) $^ -o $@ $(NOPS)

#$(NLINK) $(ARCH) $^ $(NLIBS) $(NINCL) -o $@
#$(LINK) $^ -I${CUDSS_INCLUDE} -Xlinker=${CUDSS_DIR}/libcudss_static.a 
#$(LINK) $(SNLIBINCL) $^ -o $@

#Correct usege: ./main.exe <file_name> <t for text file, b for binary, OPT def= t> <matrix type: 0 (General), 1 (Symetric), 2 (Hermitian), OPT def= 0><y - skip Eigen, s - save the resoults, n - don't skip, l - load solution from file (adv skip), OPT def= n>  <file_name_with_prev_result, OPT, ignored if prev != l>

test: main.exe
	./testing.sh $(NAME)
# \
-L. -lSparseStructures -lEigenWrapper -lCUDSSWrapper \
-cudart static \
-I${CUDSS_INCLUDE} \
-I${CTK_DIR} \
-lcuda \
-lcudss_static \
-lcublas_static \
-lcublasLt_static \
-lculibos \
-Xcompiler -O3 -o main.exe

# \
-I${CUDSS_INCLUDE} \
-Xlinker=${CUDSS_DIR}/libcudss_static.a \
-Xlinker=${CTK_DIR}/lib64/libcublas_static.a \
-Xlinker=${CTK_DIR}/lib64/libcublasLt_static.a \
-Xlinker=${CTK_DIR}/lib64/libculibos.a \
end of comment


