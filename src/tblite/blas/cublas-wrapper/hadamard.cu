#include <cuda_runtime.h>
#include "cublas_v2.h"
#include <math.h>

// this scales the element of the first vector by the element of the second vector 
__global__ static void dElementwiseScale(double* vecinout, const double* vecin,size_t ndim2)
{
  // get the global id (in vector)
  size_t id = blockIdx.x*blockDim.x+threadIdx.x; 
  if (id < ndim2) vecinout[id] *= vecin[id];
}

// this scales the element of the first vector by the element of the second vector 
__global__ static void sElementwiseScale(float* vecinout, const float* vecin,size_t ndim2)
{
  // get the global id (in vector)
  size_t id = blockIdx.x*blockDim.x+threadIdx.x; 
  if (id < ndim2) vecinout[id] *= vecin[id];
}
extern "C" void sHadamard(float* vecinout, const float* vecin, size_t ndim2){
    int blockSize, gridSize;
    blockSize = 512;
    // Number of thread blocks in grid
    gridSize = (int)ceil((float)ndim2/blockSize);

    sElementwiseScale<<<gridSize,blockSize>>>(vecinout, vecin, ndim2);
}

extern "C" void dHadamard(double* vecinout, const double* vecin, size_t ndim2){
    int blockSize, gridSize;
    blockSize = 512;
    // Number of thread blocks in grid
    gridSize = (int)ceil((float)ndim2/blockSize);

    dElementwiseScale<<<gridSize,blockSize>>>(vecinout, vecin, ndim2);
}