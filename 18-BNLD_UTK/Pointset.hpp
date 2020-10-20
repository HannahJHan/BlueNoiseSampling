#ifndef POINTSET_H
#define POINTSET_H

#include <iostream>
#include <vector>
#include <typeinfo>
#include <bitset>
#include <array>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <cstring>

template <uint D, typename T>
class Vector{
private:
	
public:	
	
	T data[D];
	
	
	Vector(T f = 0) { for(uint i=0; i<D; i++) data[i] = f; }
	Vector(const T val[D]) { for(uint i=0; i<D; i++) data[i] = val[i]; }
	Vector(T i, T j, T k=0, T l=0) { data[0] = i; data[1] = j; if (D > 2) data[2] = k; if (D > 3) data[3] = l; for(uint i=4; i<D; i++) data[i] = 0; }
	Vector(const Vector<D, T>& arg_v) { memcpy(data, arg_v.data, D*sizeof(T)); }
	
	Vector<D,T> operator+(const Vector<D,T>& v) const { Vector<D,T> res; for(uint i=0; i<D; i++) res[i] = data[i]+v[i]; return res; }
	Vector<D,T>& operator+=(const Vector<D,T>& v) { for(uint i=0; i<D; i++) data[i]+=v[i]; return *this; }
	Vector<D,T> operator-(const Vector<D,T>& v) const { Vector<D,T> res; for(uint i=0; i<D; i++) res[i] = data[i]-v[i]; return res; }
	Vector<D,T>& operator-=(const Vector<D,T>& v) { for(uint i=0; i<D; i++) data[i]-=v[i]; return *this; }
	Vector<D,T> operator-() const { Vector<D,T> res; for(uint i=0; i<D; i++) res[i] = -data[i]; return res; }
	Vector<D,T> operator*(double f) const { Vector<D,T> res; for(uint i=0; i<D; i++) res[i] = f*data[i]; return res; }
	Vector<D,T>& operator*=(double f) { for(uint i=0; i<D; i++) data[i] *= f; return *this; }
	Vector<D,T> operator*(const Vector<D,T>& v) const { Vector<D,T> res; for(int i=0; i<D; i++) res[i] = v[i]*data[i]; return res; }
	Vector<D,T>& operator*=(const Vector<D,T>& v) {for(uint i=0; i<D; i++) data[i] *= v[i]; return *this; }
	Vector<D,T> operator/(double f) const { assert(f!=0); double inv = 1.f/f; return (*this)*inv; }
	Vector<D,T>& operator/=(double f) { assert(f!=0); double inv = 1.f/f; (*this)*=inv; return *this; }
	T operator[](uint i) const { assert(i>=0 && i<D); return data[i]; }
	T& operator[](uint i) { assert(i>=0 && i<D); return data[i]; }
	
	void setAll(T d) { for(uint i=0; i<D; i++) data[i] = d; }
	void set(uint i, T d) { assert(i>=0 && i<D); data[i] = d; }

	bool operator>(const Vector<D,T>& v) { for(uint i=0; i<D;i++) { if(v[i] >= data[i]) return false; } return true; }
	bool operator<(const Vector<D,T>& v) { for(uint i=0; i<D;i++) { if(v[i] <= data[i]) return false; } return true; }
	bool operator>=(const Vector<D,T>& v) { for(uint i=0; i<D;i++) { if(v[i] > data[i]) return false; } return true; }
	bool operator<=(const Vector<D,T>& v) { for(uint i=0; i<D;i++) { if(v[i] < data[i]) return false; } return true; }
	
	bool operator==(const Vector<D,T>& v) const { 
		for(uint i=0; i<D;i++) {
			if(v[i] != data[i]) 
				return false; 
		}
		return true; 
	}
	
	bool operator!=(const Vector<D,T>& v) { 
		for(uint i=0; i<D;i++) {
			if(v[i] != data[i]) 
				return true; 
		}
		return false; 
	}
	
	void operator=(const Vector<D, T>& arg_v) { memcpy(data, arg_v.data, D*sizeof(T)); }
	
	double lengthSquare() const { double res=0; for(uint i=0; i<D; i++) res += data[i]*data[i]; return res; }
	double length() const { return sqrt(lengthSquare()); }
	
	unsigned int nbCoord() { return D; }
};

template <unsigned int D, typename T>
std::ostream& operator<<(std::ostream& os, const Vector<D, T>& obj)
{
	for(uint i=0; i<D-1; i++) os << obj[i] << "\t"; 
    os << obj[D-1];
	return os;
}

template <unsigned int D, typename T>
std::istream& operator>>(std::istream& is, Vector<D, T>& obj)
{
	for(uint i=0; i<D; i++) {
		T d;
		is >> d;
		obj.set(i, d);
	}
	return is;
}

template <unsigned int D, typename T>
inline Vector<D,T> operator*(double f, const Vector<D,T>& v) { return v*f; }

template <unsigned int D, typename T>
inline double Dot(const Vector<D,T>& v1, const Vector<D,T>& v2) { double res = 0; for(uint i=0; i<D; i++) res += v1[i]*v2[i]; return res; }

template <unsigned int D, typename T>
inline double AbsDot(const Vector<D,T>& v1, const Vector<D,T>& v2) { return std::abs(Dot(v1, v2)); }

template <unsigned int D, typename T>																	
inline Vector<D, T> Normalize(const Vector<D, T>& v) { return v/(v.Length()); }




template <unsigned int D, typename T>
class Point
{
public:
	Point() : m_pos() {}
	Point(const T& f) : m_pos(f) {}
	Point(const T val[D]) : m_pos(val) {}
	Point(const T& i, const T& j) : m_pos(i, j) {}
	Point(const T& i, const T& j, const T& k) : m_pos(i, j, k) {}
	Point(const T& i, const T& j, const T& k, const T& l) : m_pos(i, j, k ,l) {}
	Point(const Point<D, T>& arg_pt) : m_pos(arg_pt.pos()) {}
	
	void pos(const Vector<D, T>& arg_pos) { m_pos = arg_pos; }
	Vector<D, T>& pos() { return m_pos; }
	Vector<D, T> pos() const { return m_pos; }
	//T* pos() { return m_pos; }
	//Vector<D, T> pos() const { return m_pos; }
	
	bool operator==(const Point<D,T> arg_pt) const { return m_pos == arg_pt.pos(); }
	
protected:
	Vector<D, T> m_pos;
	//T m_pos[D];
};

template <unsigned int D, typename T>
std::ostream& operator<<(std::ostream& os, const Point<D, T>& obj)
{
	os << std::setprecision(15) << std::fixed;
	os << obj.pos();
	return os;
}

template <unsigned int D, typename T>
std::istream& operator>>(std::istream& is, Point<D, T>& obj)
{
	is >> obj.pos();
	return is;
}

template <typename T>
struct Domain
{
//public:
	//Domain() {}
	T pMin;
	T pMax;
};


template <unsigned int D, typename T, typename P>
class Pointset : public std::vector< P >
{
public:
	Pointset() 
	{
		for(uint d=0; d<D; d++)
		{
			domain.pMin.pos()[d] = 0;
			domain.pMax.pos()[d] = 1;
		}
		toricity=-1;
	}

	Domain<P> domain;
	int toricity;
	
protected:
};

#endif
