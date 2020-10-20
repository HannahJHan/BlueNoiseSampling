/*
 * Generation code for LDBN in the paper:
 *      Ahmed, Perrier, Coeurjolly, Ostromoukhov, Guo, Yan, Huang, and Deussen
 *      Low-Discrepancy Blue Noise Sampling
 *      SIGGRAPH Asia 2016
 *
 * Coded by Abdalla G. M. Ahmed, 2016-08-31.
 *
 * Copyright (c) 2016, Abdalla G. M. Ahmed
 * All rights reserved.
 *
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
 * either expressed or implied, of the LDBN project.
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <algorithm>

void abort(const char *message) {
    fprintf(stderr, message);
    exit(1);
}

unsigned mirror[256];
void tabulate() {
    for (unsigned i = 0; i < 256; i++) {
        mirror[i] = (i >> 7) + ((i >> 5) & 2) + ((i >> 3) & 4) + ((i >> 1) & 8)
        + ((i << 1) & 16) + ((i << 3) & 32) + ((i << 5) & 64) + ((i << 7) & 128);
    }
}

inline double phi(const unsigned &i) {
    const unsigned ONE = 0x1000000;                                             // 24 bits is considered sufficient
    const double scl = 1.0 / ONE;
    return scl * (
        mirror[(i >> 16) & 255] +
        (mirror[(i >> 8) & 255] << 8) +
        (mirror[i & 255] << 16)
    );
}

struct Point { double x, y; };
std::vector<Point> s;
unsigned char *lut;
unsigned t;                                                                     // Width of table (stored point set).

void generate(int n) {
    double inv = 1.0 / n;
    unsigned mask = t - 1;                                                      // t should be a power of 2.
    unsigned shift = log2(t);
    int i = 0;
    for (unsigned Y = 0; Y < n; Y++) {
        for (unsigned X = 0; X < n; X++) {
            unsigned index = ((Y & mask) << shift) + (X & mask);
            double u = phi( (Y & 0xfffffff0) + (lut[index] & 0xf) );
            double v = phi( (X & 0xfffffff0) + (lut[index] >> 4) );
            s[i].x = inv * (X + u);
            s[i++].y = inv * (Y + v);
        }
    }
}

void loadTable(const char *fileName) {
    FILE *file = fopen(fileName, "rb");
    int w, h, grades;
    fscanf(file, " P5 %d %d %d", &w, &h, &grades);
    if ((w != h) || (grades != 255) || ((w & (w-1)) != 0)) {
        abort("Invalid table format\n");
    }
    if (fgetc(file) != '\n') abort("Missing new line char");
    t = w;
    lut = new unsigned char [t*t];
    fread(lut, 1, t*t, file);
    fclose(file);
}

const char *USAGE_MESSAGE = "Usage: %s <data.pgm> <n>\n"
    "Produces a point set of size n x n\n";

int main(int argc,char **argv) {
    if (argc != 3) {
        fprintf(stderr, USAGE_MESSAGE, argv[0]);
        exit(1);
    }
    loadTable(argv[1]);
    unsigned n = atoi(argv[2]);
    tabulate();
    s.resize(n*n);
    clock_t t0 = clock();
    generate(n);
    clock_t t1 = clock();
    printf("%d\n", n * n);
    for (int i = 0; i < n*n; i++) printf("%0.14f %0.14f\n", s[i].x, s[i].y);
    double totalTime = (double)(t1 - t0) / CLOCKS_PER_SEC;
    fprintf(stderr,
            "generated %'d samples in %.12fs (%'d points per second)..\n",
            n*n, totalTime, (int)(n*n/totalTime)
    );
}
