/*
 * Coded by Hélène Perrier helene.perrier@liris.cnrs.fr
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
 * either expressed or implied, of the Halton project.
 */
#ifndef _UTK_SAMPLER_STRATIFIED_
#define _UTK_SAMPLER_STRATIFIED_

#include "../pointsets/Pointset.hpp"
#include <cstdlib>
#include <random>
#include <chrono>
#include <sstream>

#include "../utils.hpp"

namespace utk
{

	class SamplerStratified
	{
	protected:
	public:

		SamplerStratified() { m_jittering_range = 1.0f; setRandomSeedTime(); }

		void setRandomSeed(long unsigned int arg_seed) { m_mersenneTwister.seed(arg_seed); }
		void setRandomSeedTime() { m_mersenneTwister.seed(std::chrono::system_clock::now().time_since_epoch().count()); }
		void setJitteringRange(double arg_jittering_range) { m_jittering_range = fabs(arg_jittering_range); }

		template<unsigned int D, typename T, typename P>
		bool generateSamples(Pointset<D, T, P>& arg_pts, unsigned int arg_points)
		{
			if (!arg_pts.empty())
			{
				WARNING("SamplerStratified::generateSamples the pointset to fill is not empty, clearing it ...");
				arg_pts.clear();
			}

			if (m_jittering_range < 1.0f)
			{
				WARNING("SamplerStratified::generateSamples jittering doesn't allow points to be placed inside the each whole strata");
			}
			else if (m_jittering_range > 1.0f)
			{
				WARNING("SamplerStratified::generateSamples jittering allows points to be placed outside of their strata");
			}

			if(!isFloatingType<T>())
			{
				ERROR("SamplerStratified::generateSamples generates samples in [0, 1]^d, can only run with points of type float or double");
				return false;
			}

			/// Number of points to generate needs to be a power of D
			double nthroot = pow(arg_points, 1.0 / D);
			bool is_whole_number = floor(nthroot) == nthroot;
			int near_nthroot = int(nthroot);



			if(!is_whole_number)
			{
				std::stringstream ss_dim;
				ss_dim << D;

				std::stringstream ss_given;
				ss_given << arg_points;

				std::stringstream ss_used;
				ss_used << pow(near_nthroot, D);

				WARNING("SamplerStratified::generateSamples needs " + ss_given.str() + " (the number of points to generate) to be a power of " + ss_dim.str() + ". Using " + ss_used.str() + " points instead.");

				arg_points = pow(near_nthroot, D);
			}

			//load domain & toricity
			for(uint d = 0; d<D; d++)
			{
				arg_pts.domain.pMin.pos()[d] = 0;
				arg_pts.domain.pMax.pos()[d] = 1;
			}
			arg_pts.toricity = 1;

			arg_pts.resize(arg_points);

			std::vector<int> counts;
			counts.resize(D+1, 0); /// D+1 to avoid overflow
			for(uint i = 0; i < arg_points; i++)
			{
				for(uint d = 0; d < D; d++)
				{
					double unit = (1.0 / near_nthroot);
					double jittered = getRandom01() * m_jittering_range * unit;

					arg_pts[i].pos()[d] = counts[d] * unit + jittered;
				}

				counts[0] += 1;
				for(uint j = 0; j < D; j++)
				{
					if(counts[j] >= near_nthroot)
					{
						counts[j] = 0;
						counts[j + 1] += 1;
					}
				}
			}

			return true;
		};

	protected:
		double m_jittering_range;

		std::mt19937 m_mersenneTwister;
		double getRandom01()
		{
			return (double)m_mersenneTwister() / (double)m_mersenneTwister.max();
		}

	};

} //end namespace utk


#endif
