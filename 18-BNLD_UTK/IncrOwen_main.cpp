#include "IncrementalOwen.hpp"

int main(int argc, char**argv)
{
  Pointset<2, double, Point<2, double> > pts;
  IncrOwenSampler sampler;
  
  if(argc < 4 || argc > 7)
  {
    std::cout << "Usage: " << argv[0] << " <increase size tree> <nbpts> <output file> [sobol x indice=1] [sobol y indice=2] [seed]" << std::endl;
    return 0;
  }

  sampler.setDepthIncrement(atoi(argv[1]));
  
  if(argc >= 7)
	  sampler.setRandomSeed(atoi(argv[6]));
  
  int sobol_ids[2];
  sobol_ids[0] = 1;
  sobol_ids[1] = 2;
  if(argc >= 5)
	  sobol_ids[0] = atoi(argv[4]);
  if(argc >= 5)
	  sobol_ids[1] = atoi(argv[5]);
  sampler.setIndices(sobol_ids, 2);
  
  sampler.generateSamples(pts, atoi(argv[2]));
  
  std::ofstream file(argv[3]);
  for(int i=0; i<pts.size(); i++)
    file << pts[i] << std::endl;
  
  std::cout << "Generated " << pts.size() << " samples" << std::endl;
}