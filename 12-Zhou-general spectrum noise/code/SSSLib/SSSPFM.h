#pragma once

#include <iostream>
#include <fstream>
#include <opencv\cxcore.h>

namespace SSSLib
{
	inline void PFPfWrite(int iWidth, int iHeight, const float* iData, const char* oFileName, bool iPForPf)
	{
		FILE* oFile = fopen(oFileName, "wb");
		using namespace std;
		if (iPForPf)
			fprintf(oFile, "PF\n");
		else
			fprintf(oFile, "Pf\n");
		fprintf(oFile, "%d %d\n-1.00\n", iWidth, iHeight);
		fwrite((const char*)iData, sizeof(float)*(iPForPf?3:1), iWidth*iHeight, oFile);
		fclose(oFile);
	}

	inline void PFWrite(int iWidth, int iHeight, const float* iData, const char* oFileName)
	{
		PFPfWrite(iWidth, iHeight, iData, oFileName, true);
	}
	inline void PfWrite(int iWidth, int iHeight, const float* iData, const char* oFileName)
	{
		PFPfWrite(iWidth, iHeight, iData, oFileName, false);
	}
	inline void PfWrite(const cv::Mat& iImg, const char* oFileName)
	{
		PfWrite(iImg.cols, iImg.rows, &iImg.at<float>(0, 0), oFileName);
	}
	inline void PFWrite(const cv::Mat& iImg, const char* oFileName)
	{
		PFWrite(iImg.cols, iImg.rows, iImg.at<cv::Vec3f>(0, 0).val, oFileName);
	}

	inline cv::Mat PFPfRead(const char* iFileName)
	{
		using namespace std;
		FILE* inFile = fopen(iFileName, "rb");
		char name[100];
		int h,w;
		float f;
		fscanf(inFile, "%s\n%d %d %f\n", name, &w, &h, &f);
		cv::Mat retMat(w, h, name[1]=='F'? CV_32FC3: CV_32F);
		fread((char*)retMat.data, retMat.elemSize(), w*h, inFile);
		fclose(inFile);
		return retMat;
	}
}