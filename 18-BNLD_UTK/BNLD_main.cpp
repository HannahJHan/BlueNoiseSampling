#include "BNLDSampler.hpp"

int main(int argc, char**argv)
{
  Pointset<2, double, Point<2, double> > pts;
  BNLDSampler sampler;
  
  if(argc < 5 || argc > 7)
  {
    std::cout << "Usage: " << argv[0] << " <LUT_file> <subdivision factor K> <nbpts> <output file> [sobol indice x=1] [sobol indice y=2]" << std::endl;
    return 0;
  }
  
  sampler.setLookupTableInFile(argv[1]);
  
  uint K = atoi(argv[2]);
  sampler.setSubdivFactorK(K);
  
  std::vector<uint> sobol_ids;
  sobol_ids.resize(2);
  sobol_ids[0] = 1;
  sobol_ids[1] = 2;
  if(argc >= 6)
	  sobol_ids[0] = atoi(argv[5]);
  if(argc >= 7)
	  sobol_ids[1] = atoi(argv[6]);
  sampler.setSobolIds(sobol_ids);
  
  uint nbpts_orig = atoi(argv[3]);
  uint tilesize = K*K;
  uint level = ceil(log(nbpts_orig)/log(tilesize));
  uint nbpts = pow(tilesize, level);
  if(nbpts != nbpts_orig)
	  std::cout << nbpts_orig << " is not a multiple of " << K << "^2 ... generating " << nbpts << " samples instead" << std::endl;
  sampler.generateSamples(pts, nbpts);
  
  std::ofstream file(argv[4]);
  for(int i=0; i<pts.size(); i++)
    file << pts[i] << std::endl;
  
  std::cout << "Generated " << pts.size() << " samples" << std::endl;
}