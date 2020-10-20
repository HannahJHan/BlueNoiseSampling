
#include "DiffGauss.h"

using namespace SSSLib;

DiffGauss::DiffGauss(void)
{
}

DiffGauss::~DiffGauss(void)
{
}

void DiffGauss::SetFilter(float iDelta)
{
	float boarder = BoarderDens();
	mGaussKernel.SetKernel(iDelta);
	mFilterVal = mGaussKernel.FilterMat(mImg, GaussKernel::CLAMP);
}