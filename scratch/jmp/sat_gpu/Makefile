
INCPATH=-I/scratch/jmp/install/include \
				-I/scratch/jmp/install/cuda-5.0/include \
				-I/scratch/jmp/cuda/cudpp-master/include

LIBPATH=-L/scratch/jmp/install/lib \
				-L/scratch/jmp/install/cuda-5.0/lib

LIBS=-lcudpp

default:
	nvcc $(INCPATH) $(LIBPATH) $(LIBS) main.cu -o main
