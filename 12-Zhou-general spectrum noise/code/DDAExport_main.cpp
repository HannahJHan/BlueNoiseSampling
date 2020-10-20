#include "DDAUtil.h"

using namespace std;

int main(int argc, char** argv)
{
	srand(cvGetTickCount());
	LoadConfig(argv[1]);
	gDiff.mLoadScale = 1.0f;
	if (gRange>0)
		gDiff.Load(gInputDDA.c_str(), gRange, gAddTerm);
	else
		gDiff.Load(gInputDDA.c_str(), 1.0f, gAddTerm);

	if (gRange<=0)
	{
		gDiff.ScaleR2Density((float)gN);
		gRange = gDiff.mD;
	}
	if (gN != -1)
	{
		gDiff.ScaleDensity((float)gN);
		gN = (int)gDiff.mDensity;
	}
	gDiff.mD *= gRScale;

	gRange = gDiff.mD;
	char generalDesc[100];
	sprintf(generalDesc, "delta_%d_point_%d_range_%03d", int(gDelta*10), int(gN), int(gRange*100));
	SetFileNames(generalDesc);
	gDiff.SetFilter(gDelta);

	if (gAdaptImgName.length() != 0)
	{
		gAdaptImgOgn = cv::imread(gAdaptImgName.c_str(), 0);
		if (gAdaptImgOgn.rows == 0)
		{
			cout<<"Failed to read image "<<gAdaptImgName;
			return -1;
		}
		gAdaptImg = cv::Mat(gAdaptImgOgn.size(), CV_32F);
		for (int r=0; r<(int)gAdaptImgOgn.rows; r++)
		for (int c=0; c<(int)gAdaptImgOgn.cols; c++)
			gAdaptImg.at<float>(r, c) = SSSLib::Max(5.0f, 255.0f-float(gAdaptImgOgn.at<unsigned char>(r, c)));

		gScaleImg = gDiff.SetAdaptve(gAdaptImg);

		SetAdaptive(&gScaleImg.at<float>(0, 0),
			(float*)&(GetDerivative(gScaleImg).at<cv::Vec2f>(0, 0)),
			gScaleImg.rows,
			gScaleImg.cols);
		gAdaptAspect = float(gAdaptImg.rows)/gAdaptImg.cols;
	}

	if (gContinue)
		gSamples = LoadSamples("tmp.data");
	if (!gContinue || gSamples.size() == 0)
	{
		if (gAdaptImg.rows!=0)
			gSamples = ImportanceSamples(gAdaptImg, gN);
		else if (gGridInit)
			gSamples = GridSamples(gN);
		else if (gClusterInit)
			gSamples = ClusterSamples(gN, 1);
		else if (gClusterN>0)
			gSamples = ClusterSamples(gN, gClusterN);
		else
			gSamples = RandomSamples(gN);
	}

	cv::imwrite("init.png", PlotSamples2Mat(gSamples, 512, gAdaptAspect));

	int initTime = GetTickCount();

	InitKernel((float*)&gSamples[0], &gDiff.mImg.at<float>(0, 0), gSamples.size(), gDiff.mImg.rows, cvGetTickCount(), gDiff.mD);

	float showScale = gDiff.mImg.rows*gDiff.mImg.cols/(float)cv::sum(gDiff.mImg)[0]*0.5f;

	cv::Mat curDDA = gDiff.mImg.clone();
	cv::Mat gaussDDA = gDiff.mImg.clone();

	float error = 0.1f;
	float lasAvgError = 1.0f;
	float curSumError = 0.0f;

	int i=0;
	bool bInit = true;
	int iteCount = 0;
	int sumCount = 0;
	float ddaL2Sum = L2Sum(&gDiff.mFilterVal.at<float>(0, 0), gaussDDA.rows*gaussDDA.cols)*sqrt((float)gaussDDA.rows*gaussDDA.cols);
	if (gTerminateIterations!=0)
	while (1)
	{
		int begTime = GetTickCount();
		float kernelAvgErr = 0;
		kernelAvgErr =
			Sampling((float*)&gSamples[0],
				&curDDA.at<float>(0, 0),
				&gaussDDA.at<float>(0, 0),
				gSamples.size(),
				curDDA.rows,
				gRange/2,
				gT*gRange/2);
		sumCount++;
		cout<<"Ite "<<iteCount
			<<"\tPt "<<gSamples.size()
			<<"\tTime "<<GetTickCount()-begTime
			<<"\t"<<kernelAvgErr/ddaL2Sum
			<<endl;

		ofstream logFile(gLogFileName, ios_base::app);
		logFile<<"Pt "<<gSamples.size()
			<<"\tLoop "<<iteCount
			<<"\tError "<<kernelAvgErr*100/ddaL2Sum<<"%"
			<<"\tTime "<<(GetTickCount()-initTime)
			<<endl;

		SaveSamples("tmp.data", gSamples);
		cv::Mat samplesImg = PlotSamples2Mat(gSamples, 512, gAdaptAspect);
		cv::imshow("samples", samplesImg);

		SSSLib::PfWrite(curDDA, "tmpdda.pfm");
		cv::imshow("dda", curDDA*showScale);
		cv::imshow("target", gDiff.mImg*showScale);
		cv::waitKey(10);

		curSumError += kernelAvgErr;
		i++;
		iteCount += gInnerLoop;
		if (iteCount/gCritLoop!=(iteCount-gInnerLoop)/gCritLoop)
		{
			if (curSumError/sumCount>0.98*lasAvgError && !bInit && gStop)
				break;
			bInit = false;
			lasAvgError = curSumError/sumCount;
			curSumError = 0;
			logFile<<"Avg Error "<<curSumError/sumCount/ddaL2Sum<<endl;
			sumCount = 0;
		}
		if (gTerminateIterations>0 && iteCount>=gTerminateIterations)
			break;
	}
	if (gbIsAdapt)
		gSamples = CropAdapt(gAdaptImgOgn, gSamples);
	SaveSamples(gSampleFileName, gSamples);
	SSSLib::PfWrite(curDDA, gDDAFileName);
	cv::imwrite(gSampleImageName, PlotSamples2Mat(gSamples, 512, gAdaptAspect));
	cout<<"finish"<<endl;

	return 0;
}
