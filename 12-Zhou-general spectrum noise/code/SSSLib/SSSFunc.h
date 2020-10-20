#pragma once
#include <vector>
#include <algorithm>
#include <math.h>
//#include "RandomNumber.h"
#include "SSSSmallVec.h"

namespace SSSLib
{
	static const double Pi = 3.1415926535897932384626433832795;

	template <typename _TVal>
	int sign(_TVal val1)
	{
		if (val1==0)
			return 0;
		else if (val1>0)
			return 1;
		else
			return -1;
	}

	template <typename _TVal, typename _TFac>
	_TVal interpolate(_TVal val1, _TVal val2, _TFac intPosi) {return val1 + (val2-val1) * intPosi;}

	template <typename _TVal>
	_TVal Min(_TVal val1, _TVal val2){return val1<val2?val1:val2;}

	template <typename _TVal>
	_TVal Max(_TVal val1, _TVal val2){return val1>val2?val1:val2;}

	template <typename _TVal>
	_TVal clamp(_TVal val, _TVal lowBound, _TVal highBound) {return Min(Max(val, lowBound), highBound);}

	template <typename _TVal>
	bool within(_TVal val, _TVal lowBound, _TVal highBound) {return val>=lowBound && val<=highBound;}

	template <typename _TVal>
	_TVal sqr(_TVal val){return val*val;}

	inline int warp(int iVal, int iHighBound)
	{
		int retVal = iVal - iHighBound*(iVal/iHighBound);
		if (retVal<0)
			retVal += iHighBound;
		return retVal;
	}

	template <typename _TVal>
	_TVal warp(_TVal iVal, _TVal iHighBound)
	{
		return (_TVal)(iVal - iHighBound*floor(((double)iVal)/iHighBound));
	}

	template <typename _TVal>
	void getCDF(const _TVal* iPDF, int iN, _TVal* oCDF)
	{
		oCDF[0] = iPDF[0];
		for (int i=1; i<iN; i++)
			oCDF[i] = oCDF[i-1] + iPDF[i];
		_TVal vmax = oCDF[iN-1];
		if (vmax != 0)
		for (int i=0; i<iN; i++)
			oCDF[i] /= vmax;
	}

	template <typename _TVal>
	_TVal randf()
	{return ((_TVal)((((rand()&32767) << 15) + (rand()&32767)))) / 1073741824;}
	//{return ((_TVal)((((GetRand(omp_get_thread_num())&32767) << 15) + GetRand(omp_get_thread_num())) & 65535)) / 65536;}

	template <typename _TVal>
	int randCDF(const _TVal* iCDF, int iN)
	{
		_TVal r = randf<_TVal>();
		return upper_bound(iCDF, iCDF+iN, r) - iCDF;
	}

	template <typename _TVal>
	int randPDF(const _TVal* iPDF, int iN)
	{
		// a random floating point number
		using namespace std;
		vector<_TVal> cdf(iN);
		getCDF(iPDF, iN, (_TVal*)&cdf[0]);
		return randCDF((_TVal*)&cdf[0], iN);
	}

	inline void randGauss(float& oX, float& oY)
	{
		float u = randf<float>();
		float v = randf<float>();
		float c = sqrt(-2.0f*log(1.0f-u));
		oX = c*cos(2.0f*3.1415926535897932384626433832795f*v);
		oY = c*sin(2.0f*3.1415926535897932384626433832795f*v);
	}

	inline Vec2f randGauss()
	{
		Vec2f retVec;
		randGauss(retVec[0], retVec[1]);
		return retVec;
	}
}