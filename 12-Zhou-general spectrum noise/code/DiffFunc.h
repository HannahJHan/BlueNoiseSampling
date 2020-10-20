#pragma once

#include <opencv\cxcore.h>
#include <vector>
#include "SSSLib\SSSSmallVec.h"
#include "SSSLib\SSSGaussKernel.h"

class DiffFunc
{
public:

	DiffFunc(void);
	~DiffFunc(void);

	void Load(const char* iFileName, float iD = -1, float iAddTerm = 0.0f); // the length must be the power of 2
	void ScaleDensity(float iTargetD);
	void ScaleR2Density(float iTargetD);
	float BoarderDens();
	cv::Mat SetAdaptve(cv::Mat iAdapt);
	cv::Mat SetAdaptveOld(cv::Mat iAdapt);
	cv::Mat mImg;
	float mD;
	float mDensity;
	int GetSize(){return mImg.rows;}
	float GetGridSize(){return mD/GetSize();}
	float mLoadScale;
	SSSLib::GaussKernel::TBoarderFilter mKernelWarp;
};
