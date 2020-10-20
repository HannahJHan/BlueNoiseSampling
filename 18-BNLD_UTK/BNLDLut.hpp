#include "Pointset.hpp"
#include <map>
#include <cstdlib>
#include <fstream>

const uint BNLDSampler_BITSETSIZE=64;

class BNLDLut_Key
{
  public:
	Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point<2, std::bitset<BNLDSampler_BITSETSIZE> > > pts;
	
	BNLDLut_Key() {}
	BNLDLut_Key(Pointset< 2, std::bitset<BNLDSampler_BITSETSIZE>, Point<2, std::bitset<BNLDSampler_BITSETSIZE> > >& arg_pts) { pts=arg_pts; }

	/**
	 * Needed by std::map
	 ***/
	bool operator<(const BNLDLut_Key& k) const
	{
		for(uint i=0; i<pts.size(); i++)
		{
			for(uint d=0; d<2; d++)
			{
				if(pts[i].pos()[d].to_ulong() < k.pts[i].pos()[d].to_ulong())
					return true;

				if(pts[i].pos()[d].to_ulong() > k.pts[i].pos()[d].to_ulong())
					return false;
			}
		}
		return false;
	}
};

class BNLDLut
{
private:
	std::map< BNLDLut_Key, Vector<2, unsigned long long int> >  data; //map linking key to permutation
	uint K; //subdiv factor
	
public:
	BNLDLut() { K=4; }
	BNLDLut(uint arg_K, uint arg_NBPERMUTS) { K=arg_K; }
	
	void setK(uint arg_K) { K=arg_K; }
	
	uint size() { return data.size(); }
	void clear() { data.clear(); }
	bool fill(const std::string& name)
	{
		std::ifstream file(name);
		if(!file.is_open())
		{
			std::ofstream newfile(name);
			newfile.close();
			file.open(name);
			if(!file.is_open())
				return false;
		}
		
		if(file.eof())
			return false;
		
		std::cout << "Filling LUT . . ." << std::endl;
		while(!file.eof())
		{
			uint tilesize = pow(K, 2);
			BNLDLut_Key key;
			key.pts.resize(tilesize);
			
			std::vector<uint> keyui;
			keyui.resize(tilesize*2);
			for(uint i=0; i<tilesize*2; i++)
				file >> keyui[i];
			
			for(uint i=0; i<tilesize; i++)
			{
				for(uint d=0; d<2; d++)
					key.pts[i].pos()[d] = keyui[2*i+d];
			}
			
			if(file.eof())
				break;
			
			Vector<2, unsigned long long int> p;
			for(uint d=0; d<2; d++)
				file >> p[d];
			data.insert(std::pair< BNLDLut_Key, Vector<2, unsigned long long int> >(key, p));
		}
		
		std::cout << "Filled LUT with " << data.size() << " patterns" << std::endl;
		file.close();
		return true;
	}
	
	bool getAnOptimizedPermut(BNLDLut_Key& key, Vector<2, unsigned long long int>& resulting_permut)
	{
		//std::map< BNLDLut_Key<2>, Vector<2, uint> >::iterator it;
		auto it = data.find(key);
		if(it == data.end())
			return false;
		resulting_permut = it->second;
		return true;
	}
};
