all: program

program: cudacode.o
    g++ -o program /usr/local/cuda/lib64 -lcuda -lcudart matrixMul.cpp  cudacode.o 

cudacode.o:
    nvcc -c -arch=sm_20 matrixMul.cu 

clean: rm -rf *o program
