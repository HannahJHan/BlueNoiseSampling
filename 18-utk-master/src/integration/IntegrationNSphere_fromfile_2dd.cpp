/*
 * Hélène Perrier helene.perrier@liris.cnrs.fr
 * and David Coeurjolly david.coeurjolly@liris.cnrs.fr
 *
 * Copyright (c) 2018 CNRS Université de Lyon
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the UTK project.
 */
#include "../parameters/ParamParser_getopt.hpp"
#include "../integration/IntegrationNSphere.hpp"
#include "../io/fileIO.hpp"
#include <chrono>
#include <fstream>

#include "runIntegration.hpp"

#ifdef CIMG_ON
#include "../../externals/CImg-1.6.9/CImg.h"
#endif

#define D 2

using namespace utk;

int main(int argc, char** argv)
{
	ParamParser_getopt parser;
	IntegrationNSphere disc;

	double param_radius = 0.25;
	bool param_output_func=false;
	parser.addLongOption((char*)"radius", &param_radius, 1, assignDouble, displayDouble, (char*)"<double>\tSets the radius of the n_sphere. Must be in ]0, 0.5[", (char*)"Radius");
	parser.addLongOption((char*)"outputfunc", &param_output_func, 0, assignBoolTrue, displayBool, (char*)"Outputs the integrated function in func.dat (plottable with gnuplot)", (char*)"Output Func");
#ifdef CIMG_ON
	using namespace cimg_library;
	bool param_output_img=false;
	parser.addLongOption((char*)"outputimg", &param_output_img, 0, assignBoolTrue, displayBool, (char*)"Outputs the integrated function in func.png", (char*)"Output Img");
#endif
	initParserIntegration(&parser);
	
	parser.parse(argc, argv);
	
	if(!dealParamParserIntegration(&parser))
		return 0;
	
	disc.setRadius(param_radius);

	PointsetReader< D, double, Point<D, double> > stream;
	stream.open(param_input);
	
	IntegrationStatistics stats;

	int size = 0;
	Pointset< D, double, Point<D, double> > pts;
	while(stream.readPointset(pts))
	{
		if(param_npts != 0)
			pts.resize(param_npts);
		
		if(stats.nbpts == 0)
			stats.nbpts = pts.size();
		else if(stats.nbpts != pts.size())
			WARNING("File contains pointsets of different sizes, averaging");
		
		double analytical = 0;
		double integration = 0;
		if(!disc.compute< D, double, Point<D, double> >(pts, integration, analytical))
			return 1;
		
		stats.setAnalytical(analytical);
					
		if(param_brute)
		{
			if(param_output.empty())
				std::cout << pts.size() << " " << integration << std::endl;
			else
			{
				std::ofstream file(param_output, std::ios::out);
				file << pts.size() << " " << integration << std::endl;
				file.close();
			}
		}
		else
			stats.addValue(integration);
		
		size+=pts.size();
		pts.clear();
	}

	if(!param_brute)
	{
		stats.computeStatistics();
		
		if(param_output.empty())
			std::cout << stats << std::endl;
		else
		{
			std::ofstream file(param_output, std::ios::out);
			file << stats << std::endl;
			file.close();
		}
	}
	
	if(param_output_func)
	{
		std::ofstream file("func.dat");
		
		for(double i=-0.5; i<=0.5; i+=1.0/256.0)
		{
			for(double j=-0.5; j<=0.5; j+=1.0/256.0)
			{
				Point<2, double> pt; pt.pos()[0] = i; pt.pos()[1] = j;
				file << i << " " << j << " " << disc.sampleIntegrand(pt) << std::endl;
			}
			file << std::endl;
		}
		
		file.close();
	}
#ifdef CIMG_ON
	if(param_output_img)
	{
		CImg<double> m_image_rgb(1024, 1024, 1, 3);
		
		for(double i=-0.5; i<=0.5; i+=1.0/1024.0)
		{
			for(double j=-0.5; j<=0.5; j+=1.0/1024.0)
			{
				Point<2, double> pt; pt.pos()[0] = i; pt.pos()[1] = j;
				m_image_rgb.atXYZC((i+0.5)*1024, (j+0.5)*1024, 0, 0) = disc.sampleIntegrand(pt)*255;
				m_image_rgb.atXYZC((i+0.5)*1024, (j+0.5)*1024, 0, 1) = disc.sampleIntegrand(pt)*255;
				m_image_rgb.atXYZC((i+0.5)*1024, (j+0.5)*1024, 0, 2) = disc.sampleIntegrand(pt)*255;
			}
		}
		
		m_image_rgb.save("func.png");
	}
#endif
	
	return 0;
}
