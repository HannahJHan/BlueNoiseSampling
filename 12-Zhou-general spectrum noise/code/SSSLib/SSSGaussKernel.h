#pragma once
#include <opencv\cxcore.h>
#include <math.h>
#include "SSSFunc.h"
#include "SSSSmallVec.h"

namespace SSSLib
{
	class GaussKernel
	{
	public:
		enum TBoarderFilter
		{WARP, CLAMP, CLAMP_TO_BOARDER};
		GaussKernel(){}
		void SetKernel(float iDelta, int iReso = -1)
		{
			if (iReso == -1)
				iReso = int(ceil(iDelta*6.0f));
			if (iReso%2==0)
				iReso += 1;
			mDelta = iDelta;
			mKernel = cv::Mat(iReso, iReso, CV_32F);
			for (int r=0; r<mKernel.rows; r++)
			for (int c=0; c<mKernel.cols; c++)
			{
				float x = c-0.5f*mKernel.rows+0.5f;
				float y = r-0.5f*mKernel.rows+0.5f;
				mKernel.at<float>(r, c) = exp(-(sqr(x)+sqr(y))/2.0f/sqr(mDelta));
			}
			mScale = (float)(cv::sum(mKernel)[0]);
			mKernel /= mScale;
		}
		float GetGauss(SSSLib::Vec2f iPst)
		{
			return exp(-iPst.SqrLength()/2.0f/sqr(mDelta))/mScale;
		}
		cv::Mat FilterMat(const cv::Mat& iMat,
			TBoarderFilter iBoarderFlag,
			float iBoarder = -1.0f)
		{
			return FilterMat<float>(iMat, mKernel, iBoarderFlag, iBoarder);
		}
		template <class _TKernelEle>
		static cv::Mat FilterMat(const cv::Mat& iMat,
			const cv::Mat& iKernel,
			TBoarderFilter iBoarderFlag,
			float iBoarder)
		{
			cv::Mat retMat = cv::Mat::zeros(iMat.size(), iKernel.type());
			for (int r=0; r<retMat.rows; r++)
			for (int c=0; c<retMat.cols; c++)
			{
				for (int kr=0; kr<iKernel.rows; kr++)
				for (int kc=0; kc<iKernel.cols; kc++)
				{
					int dr = r-iKernel.rows/2+kr;
					int dc = c-iKernel.cols/2+kc;
					int coordr;
					int coordc;
					if (iBoarderFlag == WARP)
					{
						coordr = warp(dr, retMat.rows);
						coordc = warp(dc, retMat.cols);
					}
					else if (iBoarderFlag == CLAMP)
					{
						coordr = clamp(dr, 0, retMat.rows-1);
						coordc = clamp(dc, 0, retMat.cols-1);
					}

					if (iBoarderFlag == CLAMP_TO_BOARDER)
					{
						if (within(dr, 0, retMat.rows-1)
							&& within(dc, 0, retMat.cols-1))
							retMat.at<_TKernelEle>(r, c) += iMat.at<float>(dr, dc)
								*iKernel.at<_TKernelEle>(kr, kc);
						else
							retMat.at<_TKernelEle>(r, c) += iBoarder
								*iKernel.at<_TKernelEle>(kr, kc);
					}
					else
						retMat.at<_TKernelEle>(r, c) += iMat.at<float>(coordr, coordc)
							*iKernel.at<_TKernelEle>(kr, kc);
				} 
			}
			return retMat;
		}
		Vec2f GetDerivative(Vec2f iDeltaPst)
		{
			float commonscale = exp(-(sqr(iDeltaPst[0])+sqr(iDeltaPst[1]))/(2*sqr(mDelta)))
				/(-sqr(mDelta)*mScale);
			return iDeltaPst*commonscale;
		}
		cv::Mat FilterDerivative(const cv::Mat& iMat,
			TBoarderFilter iBoarderFlag,
			float iBoarder = -1.0f,
			int iReso = -1)
		{
			if (iReso == -1)
				iReso = mKernel.rows;
			cv::Mat kernel(iReso, iReso, CV_32FC2);
			for (int r=0; r<kernel.rows; r++)
			for (int c=0; c<kernel.cols; c++)
			{
				Vec2f der = GetDerivative(
					Vec2f(float(c-kernel.cols*0.5f+0.5f), float(r-kernel.rows*0.5f+0.5f)));
				kernel.at<cv::Vec2f>(r, c) = cv::Vec2f(der[0], der[1]);
			}
			return FilterMat<cv::Vec2f>(iMat, kernel, iBoarderFlag, iBoarder);
		}
		float mDelta;
		float mScale;
		cv::Mat mKernel;
	};

}