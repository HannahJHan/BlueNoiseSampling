#include "Pointset.hpp"
#include "SobolSampler.hpp"
#include <random>
#include <chrono>
	
class IncrOwenSampler
{
public:

	IncrOwenSampler() {
		leveldepth = 8;
		setRandomSeedTime();
		
		int indices[2];
		for(int i=0; i<2; i++)
			indices[i] = i+1;
		sobol.setIndices(indices, 2);
	}
	
	void setDepthIncrement(int arg_depth) { leveldepth = arg_depth; }
	
	void setRandomSeed( long unsigned int arg_seed ) { m_mersenneTwister.seed(arg_seed); }
	void setRandomSeedTime() { m_mersenneTwister.seed(std::chrono::system_clock::now().time_since_epoch().count()); }
	
	void setIndices(int* arg_indices, uint arg_dimension){ sobol.setIndices(arg_indices, arg_dimension); }

	bool generateSamples(Pointset< 2, double, Point<2, double> >& arg_pts, unsigned long long int arg_points)
	{
		if(arg_points == 0)
			return true;

		if (!arg_pts.empty())
		{
			std::cout << " warning: the pointset to fill is not empty, clearing it ..." << std::endl;
			arg_pts.clear();
		}
		
		treedepth = ceil( log2(arg_points) );
		treedepth = ceil( (double)treedepth / (double)leveldepth ) * leveldepth;
		std::cout << "Increasing tree " << ceil(treedepth/leveldepth) << " times to generate " << arg_points << " samples ..." << std::endl;
		treenbflags = pow(2, treedepth);

		Vector<2, bool*> permutations;
		
		computeTree<2>(permutations);

		Pointset<2, uint, Point<2, uint> > int_pts;
		sobol.generateSamples(int_pts, treenbflags);
		int_pts.resize(arg_points);
		
		owenscrambling<2, uint, Point<2, uint> >(int_pts, permutations);
		
		arg_pts.resize(arg_points);
		for(uint d=0; d<2; d++){
			arg_pts.domain.pMin.pos()[d] = 0;
			arg_pts.domain.pMax.pos()[d] = 1;
		}
		
		for(uint i=0; i<arg_points; i++)
		{
			for(uint d=0; d<2; d++)
				arg_pts[i].pos()[d] = ((double)int_pts[i].pos()[d] / (double)treenbflags);
		}
		
		return true;
	}

private:
	
	template<unsigned int D>
	void computeTree(Vector<D, bool*>& permutations)
	{
		for(int d=0; d<D; d++)
			permutations[d] = new bool[treenbflags];
	
		for(int i=0; i<treenbflags; i++)
			for(int d=0; d<D; d++)
				permutations[d][i] = (uint)(getRandom01()*1000)%2;
		
		for(int d=0; d<D; d++)
		{
			for(int l=0; l<treedepth; l++)
			{
				uint begin_level = pow(2, l)-1;
				uint end_level = pow(2,l+1)-1;
				uint increm = pow(2, l%leveldepth);
				
				for(int i=begin_level; i<end_level; i+=increm)
					permutations[d][i] = 0;
			}
		}
	}
	
	template<unsigned int D, typename T, typename P>
	void owenscrambling(Pointset<D, T, P>& arg_pts, Vector<D, bool*> permutations)
	{
		int ind;
		int nbBits = treedepth;

		bool* pt_digits = new bool[nbBits];
		bool* res_digits = new bool[nbBits];
	  
		for(uint d=0; d<D; d++)
		{
			bool* permut = permutations[d];

			for(uint ipt=0; ipt<arg_pts.size(); ipt++) 
			{
				uint pt = arg_pts.at(ipt).pos()[d];

				for(int i=0; i<nbBits; i++)
					pt_digits[i] = (pt >> i) & 1;

				for(int ilevel=0; ilevel<nbBits; ilevel++) {
					int begin_level = pow(2, ilevel)-1;
					ind = (pt >> (nbBits-ilevel)) + begin_level;
					res_digits[ilevel] = permut[ind] ^ pt_digits[nbBits-ilevel-1];
				}

				uint tmp = 0;
				for(int i=0; i<nbBits; i++) {
					tmp <<= 1;
					tmp |= res_digits[i];
				}

				arg_pts[ipt].pos()[d] = tmp;
			}
		}

		delete [] pt_digits;
		delete [] res_digits;
	}
	
	uint treedepth;
	uint treenbflags;
	uint leveldepth;
	SamplerSobol sobol;
	
	std::mt19937 m_mersenneTwister;
    double getRandom01()
    {
        return (double)m_mersenneTwister()/(double)m_mersenneTwister.max();
    }
};