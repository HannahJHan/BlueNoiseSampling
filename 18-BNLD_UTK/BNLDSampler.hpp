#include "Pointset.hpp"
#include "SobolSampler.hpp"
#include "BNLDLut.hpp"
#include <stack>
#include <queue>
#include <chrono>
#include <iostream>
#include <bitset>
#include <algorithm>

#define LEVEL_MAX 15

bool sortBitsetX_2D(const Point<2, std::bitset<BNLDSampler_BITSETSIZE> >& a, const Point<2, std::bitset<BNLDSampler_BITSETSIZE> >& b)
{
	return a.pos()[0].to_ullong() < b.pos()[0].to_ullong();
}
bool sortBitsetY_2D(const Point<2, std::bitset<BNLDSampler_BITSETSIZE> >& a, const Point<2, std::bitset<BNLDSampler_BITSETSIZE> >& b)
{
	return a.pos()[1].to_ullong() < b.pos()[1].to_ullong();
}

class BNLDSampler
{
public:

	BNLDSampler() {
		subdiv_factor_K = 4;
		tilesize = 16;
		lookuptable.setK(subdiv_factor_K);
		
		write_lut=false;
		generates_new_patterns=false;
		
		setRandomSeedTime();
		
		sobol_ids.resize(4);
		sobol_ids[0] = 1;
		sobol_ids[1] = 2;
		sobol_ids[2] = 3;
		sobol_ids[3] = 7;
		
		optimize_patterns = true;
		null_patterns = false;
		
		silent = false;
		
		//srand48(std::chrono::system_clock::now().time_since_epoch().count());
		//srand(std::chrono::system_clock::now().time_since_epoch().count());
	}
	
	void setSilent(bool arg_silent) { silent = arg_silent; }
	
	void useNullPermuts(bool arg_null) {  null_patterns = arg_null; }
	void useOptimizedPermuts(bool arg_optim) {  optimize_patterns = arg_optim; }
	
	void setSobolIds(const std::vector<uint>& ids) { sobol_ids.clear(); sobol_ids.resize( ids.size() ); for(uint i=0; i<ids.size(); i++) sobol_ids[i] = ids[i]; }
	void generatesNewPatterns(bool arg_generates_new_patterns) { generates_new_patterns = arg_generates_new_patterns; }
	void setRandomSeed( long unsigned int arg_seed ) { srand48(arg_seed); }
	void setRandomSeedTime() { srand48(std::chrono::system_clock::now().time_since_epoch().count()); }
	
	void setLookupTableInFile(std::string arg_lookuptable_file) { lookuptable_infile=arg_lookuptable_file; }
	void setLookupTableOutFile(std::string arg_lookuptable_file) { write_lut=true; lookuptable_outfile=arg_lookuptable_file; }
	void writeLUT(bool arg_write_lut) { write_lut=arg_write_lut; }
	
	void setSubdivFactorK(uint arg_K) { subdiv_factor_K=arg_K; tilesize=pow(arg_K, 2); lookuptable.setK(arg_K); }

	bool generateSamples(Pointset<2, double, Point<2, double> >& arg_pts, unsigned long long int arg_points)
	{
		if(arg_points == 0)
			return true;

		if (!arg_pts.empty())
		{
			std::cout << "warning: BNLD::generateSamples the pointset to fill is not empty, clearing it ..." << std::endl;
			arg_pts.clear();
		}
		
		if (generates_new_patterns && !write_lut)
		{
			std::cout << "BNLD::generateSamples will compute new LUT entries but will not write them ..." << std::endl;
		}
		if (!generates_new_patterns && write_lut)
		{
			std::cout << "BNLD::generateSamples will NOT compute new LUT entries but will write a new LUT anyway ..." << std::endl;
		}
		if (write_lut && !optimize_patterns)
		{ 
			std::cout << "BNLD::generateSamples will fill a LUT with random patterns, Is this the desired behaviour ? y/n" << std::endl;
			char c = ' ';
			while(1)
			{
				std::cin >> c;
				if (c == 'n')
					return false;
				if (c == 'y')
					break;
				std::cout << "press y or n" << std::endl;
			}
		}
		
		if(!lookuptable.fill(lookuptable_infile.c_str()))
		{
			std::cout << "Failed to open LUT";
			return false;
		}
		
		uint level_max = ceil( log(arg_points)/log(tilesize) );
		unsigned long long int nextpoweroftilesize = pow(tilesize, level_max);
	
		Pointset< 2, uint, Point<2, uint> > input_pts;
		
		int indices[2];
		indices[0] = sobol_ids[0];
		indices[1] = sobol_ids[1];
		sobol.setIndices(indices, 2);

		if(!silent)
			std::cout << "Generating " << nextpoweroftilesize << " Sobol samples ..." << std::flush;
		if(!sobol.generateSamples(input_pts, nextpoweroftilesize))
		{
			std::cout << "BNLD::generateSamples Failed to generates sobol samples ..."<< std::endl;
			return false;
		}
		if(!silent)
			std::cout << "Done" << std::endl;
		
		Pointset< 2, uint, Point<2, uint> > output_pts;
		output_pts.resize(nextpoweroftilesize);

		
		if(!initializeStrata(input_pts, level_max))
			return false;
		if(!recursivePermutations(level_max, input_pts, output_pts))
			return false;

		arg_pts.resize(arg_points);
		arg_pts.toricity = 1;
		
		for(uint d=0; d<2; d++)
		{
			arg_pts.domain.pMin.pos()[d] = 0;
			arg_pts.domain.pMax.pos()[d] = 1;
		}

		for(unsigned long long int i=0; i<arg_points; i++)
		{
			for(int d=0; d<2; d++)
			{
				arg_pts[i].pos()[0+d] = output_pts[i].pos()[d]  / (double)nextpoweroftilesize;
			}
		}
	      
	      if(!silent)
		      std::cout << "Done" << std::endl;
		      
	      lookuptable.clear();
	      
	      return true;
	}

	
private:
    
	struct Recursive
	{
		uint level;
		Vector<2, uint> cornerStrata;
		Vector<2, uint> digital_shift;
	};

	bool recursivePermutations(uint levelmax, Pointset< 2, uint, Point<2, uint> >& input_pts, Pointset< 2, uint, Point<2, uint> >& output_pts)
	{
		/***
		 * Variables
		 ***/
		uint kd = tilesize;
		uint logkd = log2(kd);
		uint logk = log2(subdiv_factor_K);
		uint bitmax = logkd*levelmax;
		unsigned long long int nb_pts_total = 1;//pow(k, levelmax);
		nb_pts_total <<= (logkd*levelmax);
		
		std::queue<Recursive> stack;
		
		Recursive r0;
		r0.level = 0;
		r0.cornerStrata = input_pts[0].pos();
		r0.digital_shift[0] = 0;
		r0.digital_shift[1] = 0;
		stack.push(r0);
		
		while(!stack.empty())
		{
			Recursive r = stack.front();
			stack.pop();
		
			uint level = r.level;
			
			//If I just do nb_strata_axis = pow(sqrt(k),level); I will face an int overflow when generating 4G samples
			unsigned long long int nb_strata_axis = 1;
			nb_strata_axis <<= (logk*level);

			//Used for retrieving the LUT key from the digits of the selected points
			uint bit1 = (logkd*levelmax)-(level*logk);
			
			/***
			* Select all samples within the current strata
			* Also saves their sobol indices to be able to properly reorder them
			***/
			std::vector<unsigned long long int> selected_ids;
			selected_ids.resize(kd);
			Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point< 2, std::bitset<BNLDSampler_BITSETSIZE> > > selected_pts; // = keyPattern16BITSTRING
			selected_pts.resize(kd);
		
			double x = r.cornerStrata[0];
			double y = r.cornerStrata[1];
			x /= nb_pts_total;
			y /= nb_pts_total;
			unsigned long long int stratum_x = x * nb_strata_axis;
			unsigned long long int stratum_y = y * nb_strata_axis;
			for(uint i=0; i<kd; i++)
			{
				selected_ids[i] = strata[level][stratum_x*nb_strata_axis+stratum_y][i];
				selected_pts[i].pos()[0] = input_pts[ selected_ids[i] ].pos()[0];
				selected_pts[i].pos()[1] = input_pts[ selected_ids[i] ].pos()[1];
			}

			/**
			* We apply the accumulated XOR values to reposition the 0 point as it would have been premutted in the previous level
			***/
			Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point< 2, std::bitset<BNLDSampler_BITSETSIZE> > > xored_selected_pts;
			xored_selected_pts.resize(kd);
			for(uint i=0; i<selected_pts.size(); i++)
			{
				xored_selected_pts[i].pos()[0] = selected_pts[i].pos()[0].to_ulong() ^ r.digital_shift[0];
				xored_selected_pts[i].pos()[1] = selected_pts[i].pos()[1].to_ulong() ^ r.digital_shift[1];
			}

			
			/**
			* Retrieving the LUT key from the digits of the xored points
			***/
			Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point< 2, std::bitset<BNLDSampler_BITSETSIZE> > > key_pts;
			key_pts.resize(kd);
			for(uint i=0; i<xored_selected_pts.size(); i++)
				for(uint b=0; b<logkd; b++)
				{
					key_pts[i].pos()[0][logkd-1-b] = xored_selected_pts[i].pos()[0][bit1-b-1]; //we read the logk most significant bits from the position (bitmax-bit0) to (bitmax-bit0-logk)
					key_pts[i].pos()[1][logkd-1-b] = xored_selected_pts[i].pos()[1][bit1-b-1];
				}
			
			/**
			* Retrieving the permutation from the LUT
			***/
			Vector<2, unsigned long long int> permuts;
			BNLDLut_Key key(key_pts);

			if(!lookuptable.getAnOptimizedPermut(key, permuts))
			{
				if(!generates_new_patterns)
				{
					std::string str;
					str += "{";
					for(uint i=0; i<key_pts.size(); i++)
					{
						str += "{";
						str += std::to_string(key_pts[i].pos()[0].to_ulong());
						str += ",";
						str += std::to_string(key_pts[i].pos()[1].to_ulong());
						str += "}";
					}
					str += "}";
					std::cout << "BNLD::recursivePermutations Couldn't find a LUT entry for pattern " + str << std::endl;
					return false;
				}
			}
			
			/**
			* Apply the permutation to the pointset
			***/
			Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point< 2, std::bitset<BNLDSampler_BITSETSIZE> > > permutted_selected_pts;
			permutted_selected_pts.resize(kd);
			applyPermut(xored_selected_pts, permutted_selected_pts, permuts, level, levelmax);
			
			if((level+1) == levelmax)
			{
				/**
				* Return the permutted pointset
				***/
				for(uint i=0; i<permutted_selected_pts.size(); i++)
				for(int d=0; d<2; d++)
					output_pts[selected_ids[i]].pos()[d] = permutted_selected_pts[i].pos()[d].to_ulong();
			}
			else
			{
				for(uint i=0; i<kd; i++)
				{
					Recursive rprime;
					rprime.level = level+1;
					
					/**
					* Computes the new digital shift
					***/
					rprime.digital_shift[0] = r.digital_shift[0] ^ ( permutted_selected_pts[i].pos()[0].to_ulong() ^ xored_selected_pts[i].pos()[0].to_ulong() );
					rprime.digital_shift[1] = r.digital_shift[1] ^ ( permutted_selected_pts[i].pos()[1].to_ulong() ^ xored_selected_pts[i].pos()[1].to_ulong() );
					
					/**
					* Computes the corner of the point's stratum
					***/
					Vector<2, uint> newcornerstrata;
					std::bitset<BNLDSampler_BITSETSIZE> pt;
					for(uint b=0; b<(level+1)*logk; b++)
						pt[(level+1)*logk - 1 -b] = selected_pts[i].pos()[0][bitmax-1-b];
					uint nb_pts_strata = nb_pts_total / pow(sqrt(kd), level+1);
					rprime.cornerStrata[0] = pt.to_ulong() * nb_pts_strata;
					for(uint b=0; b<(level+1)*logk; b++)
						pt[(level+1)*logk - 1 -b] = selected_pts[i].pos()[1][bitmax-1-b];
					rprime.cornerStrata[1] = pt.to_ulong() * nb_pts_strata;
					
					/**
					* Recurse
					***/
					stack.push(rprime);
				}
			}
		}
		
		return true;
		
	}
	
	void applyPermut(Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point< 2, std::bitset<BNLDSampler_BITSETSIZE> > > in_pts, 
					 Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point< 2, std::bitset<BNLDSampler_BITSETSIZE> > >& out_pts,
					 Vector<2, unsigned long long int> permuts,
					 uint level,
				     uint level_max
					)
	{
		uint N = in_pts.size();
		Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point< 2, std::bitset<BNLDSampler_BITSETSIZE> > > pts2;
		pts2.resize(N);
		
		uint kd = tilesize;
		uint logkd = log2(tilesize);
		
		Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point< 2, std::bitset<BNLDSampler_BITSETSIZE> > > in_pts_sortx;
		in_pts_sortx.resize(N);
		Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point< 2, std::bitset<BNLDSampler_BITSETSIZE> > > in_pts_sorty;
		in_pts_sorty.resize(N);
		for(uint i=0; i<N; i++)
		{
			in_pts_sortx[i] = in_pts[i];
			in_pts_sorty[i] = in_pts[i];
		}
		std::sort(in_pts_sortx.begin(), in_pts_sortx.end(), sortBitsetX_2D);
		std::sort(in_pts_sorty.begin(), in_pts_sorty.end(), sortBitsetY_2D);
		
		std::bitset<BNLDSampler_BITSETSIZE> xpermuts[2];
		xpermuts[0] = permuts[0];
		xpermuts[1] = permuts[1];
		
		Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point< 2, std::bitset<BNLDSampler_BITSETSIZE> > > pts_indices;
		pts_indices.resize(in_pts.size());
		for(uint ipt=0; ipt<in_pts.size(); ipt++) 
		{
			uint idx = 0;
			for(idx=0; idx<N && in_pts[ipt].pos()[0] != in_pts_sortx[idx].pos()[0]; idx++);
			uint idy = 0;
			for(idy=0; idy<N && in_pts[ipt].pos()[1] != in_pts_sorty[idy].pos()[1]; idy++);
			
			pts_indices[ipt].pos()[0] = idx;
			pts_indices[ipt].pos()[1] = idy;
		}
		
		int ind;
		
		uint NBITS = logkd;
		uint PERMUTS_NBITS = kd-1;
		
		for(uint ipt=0; ipt<tilesize; ipt++) 
		{
		  std::bitset<BNLDSampler_BITSETSIZE> xresdigits[2];
		  std::bitset<BNLDSampler_BITSETSIZE> thisSobolPt[2];
		  
		  for(uint d=0; d<2; d++)
		    thisSobolPt[d] = pts_indices[ipt].pos()[d];

		  for(uint ilevel=0; ilevel<NBITS; ilevel++) {
			  
		    for(uint d=0; d<2; d++)
		    {
		      ind = (pts_indices[ipt].pos()[d].to_ulong() >> (NBITS-ilevel)) + (pow(2,ilevel)-1);
		      xresdigits[d].set(NBITS-ilevel-1, (xpermuts[d][PERMUTS_NBITS-ind-1] ^ thisSobolPt[d][NBITS-ilevel-1]) );
		    }
			
		  }
		  
		  for(uint d=0; d<2; d++)
		    pts_indices[ipt].pos()[d] = xresdigits[d].to_ulong();
		}
		
		
		for(uint ipt=0; ipt<N; ipt++)
			out_pts[ipt].pos()[0] = in_pts_sortx[pts_indices[ipt].pos()[0].to_ullong()].pos()[0];
		
		
		for(uint ipt=0; ipt<N; ipt++)
			out_pts[ipt].pos()[1] = in_pts_sorty[pts_indices[ipt].pos()[1].to_ullong()].pos()[1];
	}
	
	bool initializeStrata(Pointset< 2, uint, Point<2, uint> >& input_pts, uint level_max)
	{
		
		strata.resize(level_max);
		uint logKd = log2(tilesize);
		uint logK = log2(subdiv_factor_K);
		unsigned long long int nb_pts_total = 1;
		nb_pts_total <<= (logKd*level_max);//pow(tilesize, level_max);
		
		for(uint l=0; l<level_max; l++)
		{	
			uint nb_strata = 1;
			nb_strata <<= (logKd*l);
			
			//unsigned long long int nb_strata_axis = sqrt( nb_strata );
			uint nb_strata_axis = 1;
			nb_strata_axis <<= (logK*l);
			
			strata[l].resize( nb_strata );
			for(uint i=0; i<nb_strata; i++)
				strata[l][i].resize(tilesize);
			
			std::vector< uint > cpt_strata;
			cpt_strata.resize(nb_strata);
			
			unsigned long long int nbpts = 1;//pow(tilesize,l+1);
			nbpts <<= (logKd*(l+1));//pow(tilesize, level_max);
			for(unsigned long long int i=0; i<nbpts; i++)
			{
				if(!silent && (i%100 == 0))
					std::cout << "Init stratas ... Level " << l << "/" << level_max << "\t" <<  (double)i/(double)nbpts * 100.0 << "%...\r" << std::flush;
				
				double x = input_pts[i].pos()[0];
				double y = input_pts[i].pos()[1];
				x /= nb_pts_total;
				y /= nb_pts_total;
				unsigned long long int stratum_x = x * nb_strata_axis;
				unsigned long long int stratum_y = y * nb_strata_axis;
				
				uint& n = cpt_strata[stratum_x*nb_strata_axis + stratum_y];
				if(n >= tilesize)
				{
					std::cout << "BNLD::initializeStrata non dyadic point repartition" << std::endl;
					return false;
				}
				strata[l][stratum_x*nb_strata_axis + stratum_y][n] = i;
				n++;
			}
			
			for(unsigned long long int i=0; i<nb_strata; i++)
				if( cpt_strata[i] != tilesize )
				{
					std::cout << "BNLD::initializeStrata non dyadic point repartition" << std::endl;
					return false;
				}
		}

		if(!silent)
				std::cout << "Done" << std::endl;
			
		return true;
	}
	
	uint subdiv_factor_K;
	uint tilesize;
	std::string lookuptable_infile;
	std::string lookuptable_outfile;
	bool write_lut;
	bool generates_new_patterns;
	bool optimize_patterns;
	bool silent;
	bool null_patterns;
	std::vector<uint> sobol_ids;
	
	SamplerSobol sobol;
	std::vector< std::vector< std::vector<uint> > > strata;
	BNLDLut lookuptable;
};
