#pragma once

#include "DiffGauss.h"
#include <iostream>
#include <stdlib.h>
#include<iostream>
#include<string>
#include <sstream>
#include<fstream>
#include<vector>
#include <opencv\highgui.h>
#include <opencv\cxcore.h>

#include "SSSLib\SSSSmallVec.h"
#include "SSSLib\SSSFunc.h"
#include "SSSLib\SSSPFM.h"

#include "kernel.h"

extern void LoadConfig(char* iFileName);
extern float CalcError(float* iData, float* iTarget, int iSize);
extern cv::Mat GetDerivative(cv::Mat iMat);

extern DiffGauss gDiff;
extern cv::Mat gAdaptImgOgn;
extern cv::Mat gAdaptImg;
extern cv::Mat gScaleImg;
extern cv::Mat gBoundaryImg;
extern float gRange;
extern int gN;
extern std::string gInputDDA;
extern bool gContinue;
extern std::string gAdaptImgName;
extern bool gElecInit;
extern bool gGridInit;
extern bool gClusterInit;
extern bool gStop;
extern int gLoop;
extern float gT;
extern float gRScale;
extern float gAddTerm;
extern int gClusterN;
extern int gTerminateIterations;
extern std::vector<SSSLib::Vec2f> gSamples;
extern int gCritLoop;

extern std::vector<SSSLib::Vec2f> LoadSamples(const char* iFileName);
extern void SaveSamples(const char* iFileName, const std::vector<SSSLib::Vec2f>& iSamples);
	
extern std::vector<SSSLib::Vec2f> GridSamples(int iN);
extern std::vector<SSSLib::Vec2f> ClusterSamples(int iN, int iClusterCount = 1);
extern std::vector<SSSLib::Vec2f> RandomSamples(int iN);
extern std::vector<SSSLib::Vec2f> ImportanceSamples(cv::Mat iAdapt, int iN);

extern cv::Mat PlotSamples2Mat(const std::vector<SSSLib::Vec2f>& iSamples, int iReso, float iAspect);

extern char gDDAFileName[100];
extern char gSampleImageName[100];
extern char gSampleFileName[100];
extern char gLogFileName[100];
extern void SetFileNames(const char* iTitle);

extern float L2Sum(float* iData, int iSize);

extern std::vector<SSSLib::Vec2f> CropAdapt(const cv::Mat& iAdapt, const std::vector<SSSLib::Vec2f>& iSamples);