// includes, system

#define _DER_ADJ

#define _DISTURB
#define _GLOBAL_DISTURB

#define _DISTURB_FACTOR 0.00015f

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <curand.h>
#include <vector_types.h>

#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "kernel.h"
#include "kernel_b.cu"
#include "BuildGrid.cu"

//#define _NON_WRAP

extern __global__ void MoveAdapt(float *data,float *data_new,size_t size,float DIS,int RESOLUTION,float CONSTANT);
extern __global__ void DDAAdapt(float *data, float *resolution,size_t size,float DIS,int RESOLUTION);

float *dev_data = 0;
float *dev_res = 0;
float *dev_target=0;
float *dev_diff=0;
float *dev_cov=0;
float *dev_cov_x=0;
float *dev_cov_y=0;
float *dev_datanew=0;
float *dev_randNum=0;
float *dev_gauss_coe = 0;
float *dev_gauss_der_coe = 0;
cudaArray* cu_array_x;
cudaArray* cu_array_y;
bool gbIsAdapt = false;
int gInnerLoop = 10;
curandGenerator_t gGen;
float gDelta = 1.5f;
int lRows;
int lCols;
float gAdaptAspect = 1.0f;
bool gbIsChangeStep = true;
bool gIsUseCell = true;
float2 *dev_force=0;
float *dev_force_len=0;
float gGridRelativeSize = 0.1f;

void InitKernel(float* data,
				 float* target,
				 size_t sampleSize,
				 size_t Resolution,
				 int seed,
				 float iDDASize)
{
	cudaMalloc((void**)&dev_datanew, sampleSize*2 * sizeof(float));
	cudaMemcpy(dev_datanew, data, sampleSize*2 * sizeof(float), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&dev_data, sampleSize*2 * sizeof(float));
	cudaMemcpy(dev_data, data, sampleSize*2 * sizeof(float), cudaMemcpyHostToDevice);
	cudaMalloc((void**)&dev_res,  Resolution*Resolution * sizeof(float));

	cudaMalloc((void**)&dev_target,  Resolution*Resolution * sizeof(float));
	cudaMemcpy(dev_target, target, Resolution*Resolution* sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc((void**)&dev_cov,  Resolution*Resolution * sizeof(float));
	cudaMalloc((void**)&dev_diff,  Resolution*Resolution * sizeof(float));

	cudaMalloc((void**)&dev_cov_x,  Resolution*Resolution* sizeof(float));
	cudaMalloc((void**)&dev_cov_y,  Resolution*Resolution* sizeof(float));

	cudaMalloc(&dev_randNum, sampleSize*2*sizeof(float));
	curandCreateGenerator(&gGen, CURAND_RNG_PSEUDO_DEFAULT);
	curandSetPseudoRandomGeneratorSeed(gGen, seed);

    cudaChannelFormatDesc channelDesc_x = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaChannelFormatDesc channelDesc_y = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaMallocArray( &cu_array_x, &channelDesc_x, Resolution,Resolution );
	cudaMallocArray( &cu_array_y, &channelDesc_y, Resolution,Resolution );
	tex_x.addressMode[0] = cudaAddressModeClamp;
	tex_x.addressMode[1] = cudaAddressModeClamp;
	tex_x.filterMode = cudaFilterModeLinear;
	tex_x.normalized = true;  
	tex_y.addressMode[0] = cudaAddressModeClamp;
	tex_y.addressMode[1] = cudaAddressModeClamp;
	tex_y.filterMode = cudaFilterModeLinear;
	tex_y.normalized = true;  

	// init gaussian coefficients
	int kernelSize = int(gDelta*3.0f)*2+1;
	cudaMalloc((void**)&dev_gauss_coe, kernelSize * kernelSize * sizeof(float));
	cudaMalloc((void**)&dev_gauss_der_coe, kernelSize * kernelSize * sizeof(float) * 2);
	_CHECK_ERROR

	const int blockSize = 8;
	dim3 grid ((kernelSize-1)/blockSize+1, (kernelSize-1)/blockSize+1);
	dim3 threads (blockSize, blockSize, 1);
	InitGaussCoe<blockSize><<<grid, threads>>>(dev_gauss_coe, dev_gauss_der_coe, kernelSize, gDelta);
	_CHECK_ERROR

	cudaMalloc((void**)&dev_force_len, sampleSize * sizeof(float));
	cudaMalloc((void**)&dev_force, sampleSize * sizeof(float2));
	_CHECK_ERROR

	if (gIsUseCell)
	{
		gCellReso = int(1.0f/(gGridRelativeSize*iDDASize));
		if (gCellReso<1)
			gCellReso = 1;
		InitGrid(sampleSize);
		_CHECK_ERROR
		printf("Uniform Grid Resolution %d\n", gCellReso);
	}
}

struct square
{
	__host__ __device__
	float operator()(float x)
	{
		return x * x;
	}
};

float snrm2_fast(float* x, int size)
{
	// with fusion
	return sqrt( thrust::transform_reduce(
		thrust::device_ptr<float>(x),
		thrust::device_ptr<float>(x+size),
		square(),
		0.0f,
		thrust::plus<float>()));
}

float Sampling(float* data,
						 float* res,
						 float* cov,
						 size_t sampleSize,
						 size_t Resolution,
						 float DDA_dis,
						 float step)
{
	dim3 grid ((sampleSize-1)/256+1, 1);
	dim3 threads (256, 1, 1);
	dim3 grid2(( Resolution*Resolution-1)/256+1,1);

	cudaChannelFormatDesc channelDesc_x = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaChannelFormatDesc channelDesc_y = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	
	float sumErr = 0;
	for(int iter=0;iter<gInnerLoop;iter++)
	{
		if (gIsUseCell && !gbIsAdapt)
			BuildGrid((float2*)dev_data, gSampleKeys, sampleSize, gCellBeg, gCellEnd);
		float error = 0;
		cudaMemset(dev_res, 0, Resolution*Resolution*sizeof(float));

		_CHECK_ERROR
		if (gbIsAdapt)
			DDAAdapt<<<grid, threads>>>(dev_data, dev_res, sampleSize,DDA_dis,Resolution, gAdaptAspect);
		else
		{
			if (gIsUseCell)
				DDAGrid<<<grid, threads>>>((float2*)dev_data,
					sampleSize,
					gCellReso,
					gCellBeg,
					gCellEnd,
					dev_res,
					Resolution,
					DDA_dis
					);
			else
				DDA<<<grid, threads>>>(dev_data, dev_res, sampleSize,DDA_dis,Resolution);
		}
		_CHECK_ERROR


		Gaussian<<<grid2, threads>>>(dev_res, dev_target,dev_cov, dev_diff,Resolution*Resolution,Resolution,int(gDelta*3.0f),gDelta,dev_gauss_coe);
		_CHECK_ERROR
		error = snrm2_fast(dev_cov, Resolution*Resolution);
		float disturbStep=snrm2_fast(dev_diff, Resolution*Resolution)/error;
		sumErr += error;

        Gaussian2<<<grid2, threads>>>(dev_cov, dev_cov_x,dev_cov_y, Resolution*Resolution,Resolution,int(gDelta*3.0f),gDelta, dev_gauss_der_coe);
		_CHECK_ERROR      
		cudaMemcpyToArray(cu_array_x, 0, 0,dev_cov_x, Resolution*Resolution*4, cudaMemcpyDeviceToDevice);
	    
		_CHECK_ERROR
		cudaMemcpyToArray(cu_array_y, 0, 0,dev_cov_y, Resolution*Resolution*4, cudaMemcpyDeviceToDevice);

		curandGenerateUniform(gGen, dev_randNum, sampleSize*2);

		cudaBindTextureToArray( tex_x, cu_array_x, channelDesc_x);
		cudaBindTextureToArray( tex_y, cu_array_y, channelDesc_y);
		if (gbIsAdapt)
			MoveAdapt<<<grid, threads>>>(dev_data, dev_datanew, dev_randNum, sampleSize,DDA_dis, Resolution, step, lRows, lCols, gAdaptAspect);
		else
		{
			if (gIsUseCell)
				CalcForceGrid<<<grid, threads>>>((float2*)dev_data,
													sampleSize,
													dev_force,
													dev_force_len,
													DDA_dis,
													gCellBeg,
													gCellEnd,
													gCellReso);
			else
				CalcForce<<<grid, threads>>>(dev_data,
					dev_force,
					dev_force_len,
					sampleSize,
					DDA_dis,
					Resolution);
			_CHECK_ERROR
			float maxPow = thrust::reduce(thrust::device_ptr<float>(dev_force_len),
				thrust::device_ptr<float>(dev_force_len+sampleSize),
				0.0f,
				thrust::maximum<float>());
			_CHECK_ERROR
			MoveForce<<<grid, threads>>>
				(dev_data,
				dev_datanew,
				sampleSize,
				thrust::raw_pointer_cast(&dev_force[0]),
				step/sqrt(maxPow),
				step*disturbStep*0.02f,
				dev_randNum);
			_CHECK_ERROR
		}
		_CHECK_ERROR
		Change<<<grid, threads>>>(dev_data, dev_datanew,sampleSize);
		_CHECK_ERROR
	}
	cudaMemcpy(res, dev_res, Resolution*Resolution * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(cov, dev_cov, Resolution*Resolution * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(data, dev_data, sampleSize*2 * sizeof(float), cudaMemcpyDeviceToHost);
	return sumErr/gInnerLoop;
}

void SetAdaptive(float* iAdapt,
				 float* iDAdapt,
				 int iRows,
				 int iCols)
{
	lRows = iRows;
	lCols = iCols;

    cudaChannelFormatDesc channelDesc_adapt = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaChannelFormatDesc channelDesc_dadapt = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
	cudaArray* cu_array_adapt;
	cudaArray* cu_array_dadapt;

	cudaMallocArray( &cu_array_adapt, &channelDesc_adapt, iCols, iRows);
	cudaMallocArray( &cu_array_dadapt, &channelDesc_dadapt, iCols, iRows);
	cudaMemcpyToArray(cu_array_adapt, 0, 0, iAdapt, iCols*iRows*4, cudaMemcpyHostToDevice);
	cudaMemcpyToArray(cu_array_dadapt, 0, 0, iDAdapt, iCols*iRows*4*2, cudaMemcpyHostToDevice);

	_CHECK_ERROR
	gAdapt.addressMode[0] = cudaAddressModeWrap;
	gAdapt.addressMode[1] = cudaAddressModeWrap;
	gAdapt.filterMode = cudaFilterModeLinear;
	gAdapt.normalized = true;

	gDAdapt.addressMode[0] = cudaAddressModeWrap;
	gDAdapt.addressMode[1] = cudaAddressModeWrap;
#ifdef _DER_ADJ
	gDAdapt.filterMode = cudaFilterModePoint;
#else // _DER_ADJ
	gDAdapt.filterMode = cudaFilterModeLinear;
#endif // _DER_ADJ
	gDAdapt.normalized = true;

	cudaBindTextureToArray( gAdapt, cu_array_adapt, channelDesc_adapt);
	cudaBindTextureToArray( gDAdapt, cu_array_dadapt, channelDesc_dadapt);
	_CHECK_ERROR

}


