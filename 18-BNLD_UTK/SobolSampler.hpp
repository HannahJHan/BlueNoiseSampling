#include "Pointset.hpp"
#include <fstream>

#define DMAX 10
#define a058947_LENGTH 32

	
class SamplerSobol
{
public:
	
	SamplerSobol() {
		
		for(int i=0; i<2; i++)
			m_indices[i] = i+1;
		
		buildMSobol(m_indices, 2);
	}
	
	void setIndices(int* arg_indices, uint arg_dimension)
	{
		if(arg_dimension > DMAX)
			std::cout << "SamplerSobol::setIndices Can only generate Sobol Samples in dimensions < 10, only taking the 10 firsts indices\n" << std::endl;
		
		std::cout << "Sobol is using indices ";
		for(uint i=0; i<arg_dimension; i++)
		{
			m_indices[i] = arg_indices[i];
			std::cout << arg_indices[i] << ",";
		}
		std::cout << std::endl;
		
		buildMSobol(m_indices, arg_dimension);
		
		
		
	}
	
	Point<2, double> generateIthPoint(uint i)
	{
		Point<2, double> pt;
		for(uint d=0; d<2; d++)
		{
			pt.pos()[d] = sobol1d(i, d);
		}
		return pt;
	}
	
	bool generateSamples(Pointset<2, uint, Point<2,uint> >& arg_pts, unsigned long long int arg_points)
	{
		if (!arg_pts.empty())
		{
			std::cout << "SamplerSobol::generateSamples the pointset to fill is not empty, clearing it ..." << std::endl;
			arg_pts.clear();
		}
		
		arg_pts.resize(arg_points);
		arg_pts.toricity = 1;
		
		unsigned long long int smax = ceil(log2(arg_points));
		smax = pow(2, smax);
		
		for(uint d=0; d<2; d++){
			arg_pts.domain.pMin.pos()[d] = 0;
			arg_pts.domain.pMax.pos()[d] = smax;
		}
		
		for(unsigned long long int i=0; i<arg_points; i++)
		{
			Point<2, double> pt = generateIthPoint(i);
			for(uint d=0; d<2; d++)
				arg_pts[i].pos()[d] = pt.pos()[d]*smax;
		}

		return true;
	}
	
protected:
	std::string m_filename_direction_numbers;

	int m_indices[DMAX] = {0};
	
	uint m_polynomial[a058947_LENGTH] = {1,11,111,1011,1101,
	10011,11001,100101,101001,101111,
	110111,111011,111101,1000011,1011011,
	1100001,1100111,1101101,1110011,10000011,
	10001001,10001111,10010001,10011101,10100111,
	10101011,10111001,10111111,11000001,11001011,
	11010011,11010101};

	//columns of the sobol matrix from diagonal to top,
	//least significant bit on the diagonal and most significant bit on the upper row
	//if m_values is always 1, we have an identity matrix.
	uint m_mvalues[DMAX][a058947_LENGTH] = { {1} };
	
	uint getValueFromMatrix(uint row, uint column, uint dimension)
	{
		if(row > column)
			return 0;
		return ithBit(m_mvalues[dimension][column], column-row);
	}

	template<uint B>
	uint fromBase10toBase(int n)
	{
		uint res = 0;
		uint power = 0;
		while(n != 0)
		{
			uint tmp = n%10;
			tmp *= B;
			res += pow(tmp, power);
			n /= 10;
			power++;
		}
		return res;
	}
	uint ithBit(uint val, uint bitid)
	{
		return (val >> bitid) & 1u;
	}
	
	void buildMSobol1D(uint ID, uint d)
	{
		assert(ID > 0);
		uint polynomial = fromBase10toBase<2>(m_polynomial[ID-1]);
		uint nbbits_polynomials = log2(polynomial)+1;

		for(uint i=0; i<a058947_LENGTH; i++)
			m_mvalues[d][i]=1;
		for(uint i=2; i<nbbits_polynomials; i++)
			m_mvalues[d][i-1] = (pow(2, (i-1)) - 1)*2 + 1;
		
		uint* val_xor_0 = (uint*) malloc(sizeof(uint)*nbbits_polynomials);
		for(uint i=nbbits_polynomials; i<=a058947_LENGTH; i++)
		{
			for(uint j=0; j<nbbits_polynomials-1; j++)
			{
				val_xor_0[j] = pow(2, (j+1)) * ithBit(polynomial,nbbits_polynomials-j-2) * m_mvalues[d][i-j-2];
			}
			
			uint val_xor_1 = m_mvalues[d][i-nbbits_polynomials];
			//std::cout << "xor1 "<< val_xor_1 << std::endl;
			
			m_mvalues[d][i-1] = val_xor_1;
			for(uint j=0; j<nbbits_polynomials-1; j++)
				m_mvalues[d][i-1] ^= val_xor_0[j];
		}
			
		free(val_xor_0);
	}
	void buildMSobol(int* IDs, int nb)
	{
		for(int i=0; i<nb; i++)
			buildMSobol1D(IDs[i], i);
	}

	double sobol1d(uint n, int d)
	{
		uint nbbits_polynomials = log2(n)+1;
		if(n==0)
			nbbits_polynomials = 1;
		
		uint* val_xor = (uint*) malloc(sizeof(uint)*nbbits_polynomials);

		for(uint i=0; i<nbbits_polynomials; i++)
			val_xor[i] = ithBit(n, i) * m_mvalues[d][i] * pow(2,nbbits_polynomials-1-i);
		
		uint res = val_xor[0];
		for(uint i=1; i<nbbits_polynomials; i++)
			res ^= val_xor[i];
		
		free(val_xor);

		return (double)res/pow(2, nbbits_polynomials);
	}
};
