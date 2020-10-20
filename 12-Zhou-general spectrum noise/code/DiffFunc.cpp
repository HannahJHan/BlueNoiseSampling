
#include "DiffFunc.h"
#include <opencv\highgui.h>
#include "SSSLib\SSSFunc.h"
#include "SSSLib\SSSPfm.h"

using namespace SSSLib;
using namespace std;

DiffFunc::DiffFunc(void)
: mD(0.1f)
, mLoadScale(1.0f)
, mKernelWarp(GaussKernel::CLAMP_TO_BOARDER)
{
}

DiffFunc::~DiffFunc(void)
{
}

void DiffFunc::Load(const char *iFileName, float iD, float iAddTerm)
{
	if (iD>0)
		mD = iD;
	mImg = PFPfRead(iFileName);

	mDensity = float(cv::sum(mImg)[0])/(mD*mD);
	float cmpDens = BoarderDens();
	float avgDens = float(cv::sum(mImg)[0])/(mImg.rows*mImg.cols);
	cout<<"Avg Density: "<<avgDens
		<<"\tBoarder Density: "<<cmpDens<<endl;
	float scale = sqr(mD)/mImg.rows/mImg.cols;
	cmpDens /= scale;
	avgDens /= scale;

}

float DiffFunc::BoarderDens()
{
	float val = 0;
	for (int i=0; i<mImg.rows-1; i++)
		val += mImg.at<float>(i, 0)
			+ mImg.at<float>(mImg.rows-1-i, mImg.cols-1)
			+ mImg.at<float>(0, i)
			+ mImg.at<float>(mImg.rows-1, mImg.cols-1-i);
	return val/(4.0f*(mImg.rows-1));
}

void DiffFunc::ScaleDensity(float iTargetD)
{
	mImg *= iTargetD/mDensity;
	mDensity = iTargetD;
}

void DiffFunc::ScaleR2Density(float iTargetD)
{
	mD = mD*sqrt(mDensity/iTargetD);
	mDensity = iTargetD;
}

cv::Mat DiffFunc::SetAdaptve(cv::Mat iAdapt)
{
	mD *= sqrt(float(iAdapt.rows)/iAdapt.cols);
	cv::Mat intense = iAdapt.clone();
	intense *= mDensity/float(cv::sum(intense)[0]);
	cv::Mat initGuess = SetAdaptveOld(intense);
	cv::Mat curGuess = initGuess.clone();
	return curGuess;
	
}

cv::Mat DiffFunc::SetAdaptveOld(cv::Mat iAdapt) // input shall be the density matrix
{
	cv::Mat turnRSqr = iAdapt.clone();
	int nonZeroCount = 0;
	for (int r=0; r<iAdapt.rows; r++)
	for (int c=0; c<iAdapt.cols; c++)
		if (iAdapt.at<float>(r, c) != 0)
			nonZeroCount ++;
	float scale = 1.0f/float(cv::sum(turnRSqr)[0])*turnRSqr.rows*turnRSqr.cols;
	turnRSqr = turnRSqr*scale;
	cv::Mat retMat = turnRSqr.clone();
	for (int r=0; r<retMat.rows; r++)
	for (int c=0; c<retMat.cols; c++)
	{

		retMat.at<float>(r, c) = sqrt(1.0f/retMat.at<float>(r, c));
	}
	return retMat;
}
