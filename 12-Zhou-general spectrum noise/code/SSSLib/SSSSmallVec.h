#pragma once

// Developed by Yahan Zhou
// Email: yahan.zhou@gmail.com

#include <assert.h>

namespace SSSLib
{
	template <typename _TVal, int _LENGTH>
	struct SmallVec;

	template <typename _TVal>
	struct Type
	{ 
		typedef SmallVec<_TVal, 2> Vec2;
		typedef SmallVec<_TVal, 3> Vec3;
		typedef SmallVec<_TVal, 4> Vec4;
	};

	typedef Type<float>::Vec2 Vec2f;
	typedef Type<float>::Vec3 Vec3f;
	typedef Type<float>::Vec4 Vec4f;

	typedef Type<int>::Vec2 Vec2i;
	typedef Type<int>::Vec3 Vec3i;
	typedef Type<int>::Vec4 Vec4i;

	typedef Type<double>::Vec2 Vec2lf;
	typedef Type<double>::Vec3 Vec3lf;
	typedef Type<double>::Vec4 Vec4lf;

	typedef Type<double>::Vec2 Vec2d;
	typedef Type<double>::Vec3 Vec3d;
	typedef Type<double>::Vec4 Vec4d;

	typedef Type<unsigned int>::Vec2 Vec2ui;
	typedef Type<unsigned int>::Vec3 Vec3ui;
	typedef Type<unsigned int>::Vec4 Vec4ui;

	typedef Type<char>::Vec2 Vec2c;
	typedef Type<char>::Vec3 Vec3c;
	typedef Type<char>::Vec4 Vec4c;

	typedef Type<unsigned char>::Vec2 Vec2uc;
	typedef Type<unsigned char>::Vec3 Vec3uc;
	typedef Type<unsigned char>::Vec4 Vec4uc;

	typedef Type<short>::Vec2 Vec2s;
	typedef Type<short>::Vec3 Vec3s;
	typedef Type<short>::Vec4 Vec4s;

	typedef Type<unsigned short>::Vec2 Vec2us;
	typedef Type<unsigned short>::Vec3 Vec3us;
	typedef Type<unsigned short>::Vec4 Vec4us;

	template <typename _TVal, int _LENGTH>
	struct _VecData
	{
		union
		{
			_TVal val[_LENGTH];
			struct
			{
				_TVal x;
				_TVal y;
				_TVal z;
				_TVal w;
			};
			struct
			{
				_TVal r;
				_TVal g;
				_TVal b;
				_TVal a;
			};
		};
	};

	template <typename _TVal>
	struct _VecData<_TVal, 3>
	{
		union
		{
			_TVal val[3];
			struct
			{
				_TVal x;
				_TVal y;
				_TVal z;
			};
			struct
			{
				_TVal r;
				_TVal g;
				_TVal b;
			};
		};
	};

	template <typename _TVal>
	struct _VecData<_TVal, 2>
	{
		union
		{
			_TVal val[2];
			struct
			{
				_TVal x;
				_TVal y;
			};
			struct
			{
				_TVal r;
				_TVal g;
			};
		};
	};

	template <typename _TVal, int _LENGTH>
	struct SmallVec//: _VecData<_TVal, _LENGTH>
	{
		_TVal val[_LENGTH];
		typedef _TVal t_Ele;
		typedef SmallVec<_TVal, _LENGTH> t_Vec;
		static const int size = _LENGTH;

		_TVal& operator[](int i){assert(i<size); return val[i];}
		const _TVal& operator[](int i) const{assert(i<size); return val[i];}

		template <typename _iTVal>
		SmallVec(const SmallVec<_iTVal, _LENGTH>& iVec)
		{
			for (int i=0; i<_LENGTH; i++)
				val[i] = (_TVal)iVec[i];
		}
		SmallVec()
		{}
		SmallVec(const _TVal& ix,
			const _TVal& iy)
		{val[0] = ix; val[1] = iy;}
		SmallVec(const _TVal& ix,
			const _TVal& iy,
			const _TVal& iz)
		{val[0] = ix; val[1] = iy; val[2] = iz;}
		SmallVec(const _TVal& ix,
			const _TVal& iy,
			const _TVal& iz,
			const _TVal& iw)
		{val[0] = ix; val[1] = iy; val[2] = iz; val[3] = iw;}
		template <typename _TVec>
		SmallVec(const _TVec& iVec)
		{
			for (int i=0; i<size; i++)
				val[i] = iVec[i];
		}
		_TVal SqrLength() const
		{
			_TVal retVal = 0;
			for (int i=0; i<_LENGTH; i++)
				retVal += val[i] * val[i];
			return retVal;
		}
		_TVal Length() const
		{
			return sqrt(SqrLength());
		}

		template <typename _TVec>
		SmallVec<_TVal, _LENGTH>& operator=(const _TVec& iVec)
		{
			for (int i=0; i<_LENGTH; i++)
				val[i] = iVec[i];
			return *this;
		}
		template <typename _TVec>
		bool operator==(const _TVec& iVec) const
		{
			for (int i=0; i<_LENGTH; i++)
			if (val[i]!=iVec[i])
				return false;
			return true;
		}
		template <typename _TVec>
		bool operator<(const _TVec& iVec) const
		{
			for (int i=0; i<_LENGTH; i++)
			if (val[i]<iVec[i])
				return true;
			else if (val[i]>iVec[i])
				return false;
			return false;
		}
		template <typename _TVec>
		bool operator>(const _TVec& iVec) const
		{
			for (int i=0; i<_LENGTH; i++)
			if (val[i]>iVec[i])
				return true;
			else if (val[i]<iVec[i])
				return false;
			return false;
		}
		template <typename _TVec>
		bool operator!=(const _TVec& iVec)
		{
			return !((*this)==iVec);
		}
		SmallVec<_TVal, _LENGTH> operator-() const
		{
			SmallVec<_TVal, _LENGTH> retVec;
			for (int i=0; i<_LENGTH; i++)
				retVec[i] = -val[i];
			return retVec;
		}
		t_Vec Norm()
		{
			_TVal len = sqrt(SqrLength());
			if (len!=0)
				return (*this)/len;
			else
				return *this;
		}

#define _OPE_IN(ope)\
		SmallVec<_TVal, _LENGTH>& operator##ope##=(const t_Ele& iVal)\
		{\
			for (int i=0; i<_LENGTH; i++)\
				val[i] ##ope##= iVal;\
			return *this;\
		}\
		template <typename _TPVal>\
		SmallVec<_TVal, _LENGTH>& operator##ope##=(const _TPVal* iVal)\
		{\
			for (int i=0; i<_LENGTH; i++)\
				val[i] ##ope##= iVal[i];\
			return *this;\
		}\
		template <typename _TVal2>\
		SmallVec<_TVal, _LENGTH>& operator##ope##=(const SmallVec<_TVal2, _LENGTH>& iVec)\
		{\
			for (int i=0; i<_LENGTH; i++)\
			val[i] ##ope##= iVec[i];\
			return *this;\
		}\
		template <typename _TVal2>\
		SmallVec<_TVal, _LENGTH> operator##ope(const SmallVec<_TVal2, _LENGTH>& iVec) const\
		{\
			SmallVec<_TVal, _LENGTH> retVec;\
			for (int i=0; i<_LENGTH; i++)\
				retVec[i] = val[i] ##ope iVec[i];\
			return retVec;\
		}\
		template <typename _TPVal>\
		SmallVec<_TVal, _LENGTH> operator##ope(const _TPVal* iVec) const\
		{\
			SmallVec<_TVal, _LENGTH> retVec;\
			for (int i=0; i<_LENGTH; i++)\
				retVec[i] = val[i] ##ope iVec[i];\
			return retVec;\
		}\
		SmallVec<_TVal, _LENGTH> operator##ope(const _TVal& iVal) const\
		{\
			SmallVec<_TVal, _LENGTH> retVec;\
			for (int i=0; i<_LENGTH; i++)\
				retVec[i] = val[i] ##ope iVal;\
			return retVec;\
		}

		_OPE_IN(+);
		_OPE_IN(-);
		_OPE_IN(*);
		_OPE_IN(/);
		_OPE_IN(%);
#undef _OPE_IN

	};
 
#define _OPE_OUT(ope)\
	template <typename _TVal, int _LENGTH>\
	SmallVec<_TVal, _LENGTH> operator##ope(const typename SmallVec<_TVal, _LENGTH>::t_Ele& iVal,\
		const SmallVec<_TVal, _LENGTH>& iVec)\
	{\
		SmallVec<_TVal, _LENGTH> retVec;\
		for (int i=0; i<iVec.size; i++)\
			retVec[i] = iVec[i] ##ope iVal;\
		return retVec;\
	}

	_OPE_OUT(+);
	_OPE_OUT(-);
	_OPE_OUT(*);
	_OPE_OUT(/);
	_OPE_OUT(%);
#undef _OPE_OUT

	template <typename _TResult, typename _TVecA, typename _TVecB>
	_TResult dot(const _TVecA& iVeca, const _TVecB& iVecb)
	{
		assert(iVeca.size == iVecb.size);
		_TResult retVal = 0;
		for (int i=0; i<iVeca.size; i++) retVal += iVeca[i]*iVecb[i];
		return retVal;
	}

	template <typename _TResult, typename _TVecA, typename _TVecB>
	_TResult dis(const _TVecA& iVeca, const _TVecB& iVecb)
	{
		assert(iVeca.size == iVecb.size);
		return sqrt((iVeca-iVecb).SqrLength());
	}

	template <typename _TResult, typename _TVecA, typename _TVecB>
	_TResult cross2(const _TVecA& iVeca, const _TVecB& iVecb)
	{return iVeca[0]*iVecb[1] - iVeca[1]*iVecb[0];}

	template <typename _TVecA, typename _TVecB>
	_TVecA cross3(const _TVecA& iVeca, const _TVecB& iVecb)
	{
		_TVecA retVec;
		for (int i=0; i<3; i++)
			retVec[i] = iVeca[(i+1)%3]*iVecb[(i+2)%3] - iVeca[(i+2)%3]*iVecb[(i+1)%3];
		return retVec;
	}
}