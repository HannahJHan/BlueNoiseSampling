#include "SamplerGradientDescent.hpp"
#include "../parameters/ParamParser_getopt.hpp"
#include "../io/fileIO.hpp"
#include <chrono>

#include "runSampler.hpp"

using namespace utk;

typedef double T;
#define D 2
typedef Point<D, T> P;
typedef SamplerGradientDescent S;

int main(int argc, char** argv)
{
	ParamParser_getopt parser;
	S sampler;
	
	std::string initpts_file="";
	
	//PARSE PARAM
	initParserSampler(&parser);
	parser.addShortOption('i', &initpts_file, 1, assignString, displayString, (char*)"[string]\t\tThe pointset to optimize", (char*)"Pointset");
	//PARSING
	parser.parse(argc, argv);
	
	if(!dealParamParserSampler(&parser))
		return 0;

	if(initpts_file.empty())
	{
		ERROR("Parameter -i mandatory");
		std::cout << parser.getHelp() << std::endl;
		return 0;
	}
	
	sampler.setInputPointset(initpts_file);
	
	PointsetWriter<D, T, P> writer;
	writer.open(param_output.c_str());
	while(param_nbrealisations>0)
	{
		//SAMPLE
		Pointset<D, T, P> pts;
		
		std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
		if(!sampler.generateSamples<D, T, P>(pts))
			return 1;
		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		if(param_verbose)
			std::cout << std::fixed << std::setprecision(5) << "Generated " << pts.size() << " samples in " << time_span.count() << " secs" << std::endl;

		//WRITE
		writer.writePointset(pts);
		param_nbrealisations--;
	}
	writer.close();
	
	return 0;
}
