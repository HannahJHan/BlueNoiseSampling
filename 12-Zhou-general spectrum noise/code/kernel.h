#pragma once

#define _NON_WRAP

#include "stdio.h"

#ifdef _DEBUG
#define _CHECK_ERROR printf("%s\n", cudaGetErrorString(cudaGetLastError()));
#else
#define _CHECK_ERROR
#endif

extern void InitKernel(float* data,
				 float* target,
				 size_t sampleSize,
				 size_t Resolution,
				 int seed,
				 float iDDASize);

extern float Sampling(float* data,
						 float* res,
						 float*cov,
						 size_t sampleSize,
						 size_t Resolution,
						 float DDA_dis,
						 float step);
extern void SetAdaptive(float* iAdapt,
				 float* iDAdapt,
				 int iRows,
				 int iCols);
extern bool gbIsAdapt;
extern int gInnerLoop;
extern float gDelta;
extern float gAdaptAspect;
extern bool gbIsChangeStep;
extern int gCellReso;
extern float gGridRelativeSize;