#include "DDAUtil.h"

using namespace SSSLib;
using namespace std;

DiffGauss gDiff;
cv::Mat gAdaptImgOgn;
cv::Mat gAdaptImg;
cv::Mat gScaleImg;
cv::Mat gBoundaryImg;
float gRange = 1.0f;
int gN = -1;
string gInputDDA = "dda.pfm";
bool gContinue = false;
string gAdaptImgName;
bool gElecInit = false;
bool gGridInit = false;
bool gClusterInit = false;
bool gStop = true;
int gLoop = 1000;
float gT = 0.02f;
float gRScale = 1.0f;
float gAddTerm = 0;
int gClusterN = 0;
int gTerminateIterations = -1;
vector<Vec2f> gSamples;
int gCritLoop = 200;

char gDDAFileName[100];
char gSampleImageName[100];
char gSampleFileName[100];
char gLogFileName[100];

void LoadConfig(char* iFileName)
{
	ifstream inFile(iFileName);
	if (!inFile.is_open())
		return;
	while (!inFile.eof())
	{
		string name;
		inFile>>name;
		if (name == "Range")
			inFile>>gRange;
		else if (name == "Points")
			inFile>>gN;
		else if (name == "Delta")
			inFile>>gDelta;
		else if (name == "DDA")
			inFile>>gInputDDA;
		else if (name == "Loop")
			inFile>>gLoop;
		else if (name == "DeltaT")//temperature
			inFile>>gT;
		else if (name == "Continue")
		{
			string val;
			inFile>>val;
			if (val=="y")
				gContinue = true;
			else
				gContinue = false;
		}
		else if (name == "Adaptive")
		{
			inFile>>gAdaptImgName;
			gbIsAdapt = true;
		}
		else if (name == "ElecInit")
		{
			string val;
			inFile>>val;
			if (val == "y")
				gElecInit = true;
			else
				gElecInit = false;
		}
		else if (name == "GridInit")
		{
			string val;
			inFile>>val;
			if (val == "y")
				gGridInit = true;
			else
				gGridInit = false;
		}
		else if (name == "ClusterInit")
		{
			string val;
			inFile>>val;
			if (val == "y")
				gClusterInit = true;
			else
				gClusterInit = false;
		}
		else if (name == "ClusterN")
			inFile>>gClusterN;
		else if (name == "InnerLoop")
			inFile>>gInnerLoop;
		else if (name == "RScale")
			inFile>>gRScale;
		else if (name == "AddTerm")
			inFile>>gAddTerm;
		else if (name == "GridSize")
			inFile>>gGridRelativeSize;
		else if (name == "AutoStop")
		{
			string val;
			inFile>>val;
			if (val == "y")
				gStop = true;
			else
				gStop = false;
		}
		else if (name == "IteCount")
			inFile>>gTerminateIterations;
		else if (name == "CritLoop")
			inFile>>gCritLoop;
	}
}

float CalcError(float* iData, float* iTarget, int iSize)
{
	using namespace SSSLib;
	float retVal = 0;
	float sum = 0;
	int len = iSize;
	for (int i=0; i<iSize*iSize; i++)
	{
		int r = i/len;
		int c = i%len;
		if (
			sqr((r+0.5f)/len-1.0f)+sqr((c+0.5f)/len-1.0f)
			<sqr(0.45f))
		{
			retVal += SSSLib::sqr(iData[i]-iTarget[i]);
			sum += SSSLib::sqr(iTarget[i]);
		}
	}
	return sqrt(retVal/sum);
}

cv::Mat GetDerivative(cv::Mat iMat)
{
	cv::Mat retMat(iMat.size(), CV_32FC2);
	for (int r=0; r<retMat.rows; r++)
	for (int c=0; c<retMat.cols; c++)
		retMat.at<cv::Vec2f>(r, c) = cv::Vec2f(
			(iMat.at<float>(r, warp(c, iMat.cols))
				-iMat.at<float>(r, warp(c-1, iMat.cols)))
				*(float)iMat.cols,
			(iMat.at<float>(warp(r, iMat.rows), c)
				-iMat.at<float>(warp(r-1, iMat.rows), c))
				*(float)iMat.rows);
	return retMat;
}

vector<Vec2f> LoadSamples(const char* iFileName)
{
	vector<Vec2f> retVec;
	int num;
	ifstream inputFile(iFileName, ios_base::binary);
	inputFile.read((char*)&num, 4);
	retVec.resize(num);
	inputFile.read((char*)&retVec[0], 8*num);
	return retVec;
}

void SaveSamples(const char* iFileName, const std::vector<SSSLib::Vec2f>& iSamples)
{
	ofstream outputFile(iFileName, ios_base::binary);
	int num = iSamples.size();
	outputFile.write((char*)&num, 4);
	outputFile.write((char*)&iSamples[0], 8*iSamples.size());
}

std::vector<SSSLib::Vec2f> GridSamples(int iN)
{
	vector<Vec2f> retVec;
	retVec.resize(iN);
	int count = int(ceil(sqrt((double)iN)));
	float step = 1.0f/count;
	float x = step*0.5f;
	float y = step*0.5f;
	for (int i=0; i<(int)retVec.size(); i++)
	{
		retVec[i] = Vec2f(x, y);
		x += step;
		if (x>1.0f)
		{
			x = step*0.5f;
			y += step;
		}
	}
	return retVec;
}

std::vector<SSSLib::Vec2f> ClusterSamples(int iN, int iClusterCount)
{
	vector<Vec2f> retVec = RandomSamples(iN);
	int numClusters = SSSLib::sqr(iClusterCount);
#define _GETGRID(k) ((k)%iClusterCount+0.5f)*(1.0f/iClusterCount)
	for (int i=0; i<(int)retVec.size(); i++)
		retVec[i] = retVec[i]/50.0f
			+SSSLib::Vec2f(_GETGRID(i/numClusters), _GETGRID(i%numClusters));
#undef _GETGRID
	return retVec;
}

std::vector<SSSLib::Vec2f> RandomSamples(int iN)
{
	vector<Vec2f> retVec(iN);
	for (int i=0; i<iN; i++)
		retVec[i] = Vec2f(randf<float>(), randf<float>());
	return retVec;
}

std::vector<SSSLib::Vec2f> ImportanceSamples(cv::Mat iAdapt, int iN)
{
	vector<Vec2f> retVec(iN);
	vector<double> importance(iAdapt.rows*iAdapt.cols);
	vector<double> cum(importance.size());
	for (int i=0; i<(int)importance.size(); i++)
		importance[i] = (double)iAdapt.at<float>(i/iAdapt.cols, i%iAdapt.cols);
	getCDF<double>(&importance[0], importance.size(), &cum[0]);
	for (int i=0; i<iN; i++)
	{
		Vec2f pt(randf<float>(), randf<float>());
		int num = randCDF<double>(&cum[0], importance.size());
		pt[0] += float(num%iAdapt.cols);
		pt[1] += float(num/iAdapt.cols);
		pt[0] /= iAdapt.cols;
		pt[1] /= iAdapt.rows;
		retVec[i] = pt;
	}
	return retVec;
}

cv::Mat PlotSamples2Mat(const std::vector<SSSLib::Vec2f>& iSamples,
		int iReso,
		float iAspect
		)
{
	cv::Mat img = cv::Mat::zeros(int(iReso*iAspect), iReso, CV_8U);
	for (int i=0; i<(int)iSamples.size(); i++)
	{
		int x = clamp(int(iSamples[i][0]*img.cols), 0, img.cols-1);
		int y = clamp(int(iSamples[i][1]*img.rows), 0, img.rows-1);
		img.at<unsigned char>(y, x) = 255;
	}
	return img;
}

void SetFileNames(const char* iTitle)
{
	sprintf(gDDAFileName, "%s_DDA.pfm", iTitle);
	sprintf(gSampleImageName, "%s_samples.png", iTitle);
	sprintf(gSampleFileName, "%s.data", iTitle);
	sprintf(gLogFileName, "%s.log", iTitle);
}

float L2Sum(float* iData, int iSize)
{
	float retVal = 0;
	for (int i=0; i<iSize; i++)
		retVal += SSSLib::sqr(iData[i]);
	return sqrt(retVal/iSize);
}

vector<Vec2f> CropAdapt(const cv::Mat& iAdapt, const std::vector<SSSLib::Vec2f>& iSamples)
{
	vector<Vec2f> retVec;
	for (int i=0; i<(int)iSamples.size(); i++)
	{
		int r = int(iSamples[i][1]*iAdapt.rows);
		int c = int(iSamples[i][0]*iAdapt.cols);
		if (iAdapt.at<unsigned char>(r, c) != 255)
			retVec.push_back(iSamples[i]);
	}
	return retVec;
}