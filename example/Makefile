#DB Add MaCh3 methodology of finding CUDA version and compute capability
CUDAVER=$(shell nvcc --version | grep -o 'V[0-9].*' |  cut -d. -f1 |sed 's:V::')

ifeq ($(shell expr $(CUDAVER) \< 9), 1)
ARCH=	-gencode arch=compute_30,code=sm_30 \
	-gencode arch=compute_32,code=sm_32 \
	-gencode arch=compute_35,code=sm_35 \
	-gencode arch=compute_35,code=compute_35
else
ifeq ($(shell expr $(CUDAVER) \>= 11), 1)
ARCH= 	-gencode arch=compute_52,code=sm_52 \
	-gencode arch=compute_60,code=sm_60 \
	-gencode arch=compute_61,code=sm_61 \
	-gencode arch=compute_70,code=sm_70 \
	-gencode arch=compute_75,code=sm_75 \
	-gencode arch=compute_80,code=sm_80 \
	-gencode arch=compute_86,code=sm_86 \
	-gencode arch=compute_86,code=compute_86
else
ARCH=	-gencode arch=compute_30,code=sm_30 \
	-gencode arch=compute_32,code=sm_32 \
	-gencode arch=compute_35,code=sm_35 \
	-gencode arch=compute_35,code=compute_35
endif
endif

# Had to modify ptx assembler option
all:
	nvcc -g -O0 -x cu $(ARCH) -lineinfo -std=c++11 -Xcompiler="-fopenmp -Wall" -I.. main.cpp -o maingpu

cpu:
	g++ -g -O0 -std=c++11 -fopenmp -Wall -I.. main.cpp -o maincpu

clean:
	rm -f maingpu maincpu
