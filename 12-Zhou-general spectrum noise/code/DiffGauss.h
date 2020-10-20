#pragma once

#include "DiffFunc.h"
#include "SSSLib\SSSSmallVec.h"
#include "SSSLib\SSSGaussKernel.h"


class DiffGauss :
	public DiffFunc
{
public:
	DiffGauss(void);
	~DiffGauss(void);
	void SetFilter(float iDelta);
	SSSLib::GaussKernel mGaussKernel;
	cv::Mat mFilterVal;
	cv::Mat mFilterDev;
	cv::Mat mFilterWeighted;
};
