This code is the code for the paper
"Multi-Dimensional Sequences with Low-Discrepancy Blue-Noise 2-D Projections"
submission Id 1019

It contains the code of our BNLD sampler in arbitray
dimension.

More precisely, it contains the LUT optimization based
sampler we propse for one or two 2-D projections and
the incremental Owen sampler for other dimensions.


########################################################
## To compile both sampler
########################################################

mkdir build
cd build
cmake ..
make

or simply:

c++ -o BNLD_main BNLD_main.cpp -I. --std=c++11
c++ -o IncrOwen_main IncrOwen_main.cpp -I. --std=c++11


########################################################
## Quick example: 6-D samples with optimized projections
## on the first 4 dimensions (using Sobol 1-2 and 3-7)
## and our Owen for the remaining dimensions.
########################################################

./BNLD_main ../LUT/LUT_K4_16Mpts_Sobol12.dat 4 4096 part0.dat 1 2
./BNLD_main ../LUT/LUT_K4_16Mpts_Sobol37.dat 4 4096 part1.dat 3 7
./IncrOwen_main 8 4096 part2.dat 4 5
../cat_files.sh part0.dat part1.dat pts_4D.dat
../cat_files.sh pts_4D.dat part3.dat pts_6D.dat

cat pts_6D.dat
0.000000000000000	0.000000000000000	0.000000000000000	0.000000000000000	0.000000000000000	0.000000000000000
0.750000000000000	0.875000000000000	0.687500000000000	0.687500000000000	0.953125000000000	0.500000000000000
0.250000000000000	0.625000000000000	0.875000000000000	0.812500000000000	0.679687500000000	0.835937500000000
0.562500000000000	0.250000000000000	0.312500000000000	0.375000000000000	0.253906250000000	0.398437500000000
0.125000000000000	0.812500000000000	0.375000000000000	0.875000000000000	0.582031250000000	0.878906250000000
...

########################################################
## To run BNLD
########################################################

cd build
./BNLD_main <LUT> <Factor K> <#Pts> <output> [sobol id x =1] [sobol id y =2]

This sampler generates 2D scrambling of the sobol sequence obtained
from the primitive polynomials given by sobol id x and sobol id y
(see http://oeis.org/A058947).
We give with this sampler the LUTs for sobol pairs 1,2 and 3,7 letting you generate up to 16*10^6 samples.


########################################################
## To run Incremental Owen
########################################################

cd build
./IncrOwen_main <Depth Increment> <#Pts> <output> [sobol id x =1] [sobol id y =2] [random seed]

This sampler generates our incremental Owen's scrambling of 2D sobol sequence obtained
from the primitive polynomials given by sobol id x and sobol id y
(see http://oeis.org/A058947).
The parameter depth increment controls the depth subtrees to concatene. If this parameter is
too low, the resulting scrambling will not be aliasing free. From our experiement,
we suggest you use depth increment = 8.

For nD sampling, you should generate yourself the 2D projections using different sobol ids
and join them using the cat_files.sh script

../cat_files <input1> <input2> <output>
careful, if output exists, it will be overwritten
