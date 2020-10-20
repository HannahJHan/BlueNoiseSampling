/*
 * Reference-matching code for LDBN in the paper:
 *      Ahmed, Perrier, Coeurjolly, Ostromoukhov, Guo, Yan, Huang, and Deussen
 *      Low-Discrepancy Blue Noise Sampling
 *      SIGGRAPH Asia 2016
 *
 * Coded by Abdalla G. M. Ahmed, 2016-09-19.
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
#include <algorithm>
#include <string>


std::vector<float> u, v;
std::vector<unsigned> Xm, Ym;
std::vector<unsigned char> lut;
unsigned t;
const unsigned m = 16;                                                          // We hard-code this for simplicity
const unsigned mirror[] = {                                                        // Bit reversals over 4-bits.
    0, 8, 4, 12,
    2, 10, 6, 14,
    1, 9, 5, 13,
    3, 11, 7, 15
};

struct TCompare {
    float *p;
    bool operator() (int i1, int i2) { return (p[i1] < p[i2]); }
};

void sortIndexes(std::vector<float> &x, std::vector<unsigned> &order) {        // Returns the indexes in ascending order of entry values.
    int count = x.size();
    TCompare cmp;
    cmp.p = x.data();
    for (int i = 0; i < count; i++) order[i] = i;
    std::sort(order.begin(), order.end(), cmp);
}

void algorithm2() {
    std::vector<float> d(m);                                                   // Distances from respective edges
    std::vector<unsigned> order(m);
    for (int X = 0; X < t; X++) {                                               // Scan columns
        for (int Y0 = 0; Y0 < t; Y0 += m) {                                     // Scan the current chunk
            for (int i = 0; i < m; i++) {                                       // Extract the offsets in this chunk from the reference set
                int Y = Y0 + i;
                int index = Y * t + X;
                d[i] = u[index];
            }
            sortIndexes(d, order);
            for (int i = 0; i < m; i++) {
                unsigned j = order[i];                                          // The i'th element in the reference (ascending order)
                int index = (Y0 + j) * t + X;                                   // Translate to global index
                Xm[index] = mirror[i];                                          // place in this index the sequence number of the ith element of vdC in the chunk
            }
        }
    }
    for (int Y = 0; Y < t; Y++) {                                               // Repeat for rows
        for (int X0 = 0; X0 < t; X0 += m) {
            for (int i = 0; i < m; i++) {
                int X = X0 + i;
                int index = Y * t + X;
                d[i] = v[index];
            }
            sortIndexes(d, order);
            for (int i = 0; i < m; i++) {
                unsigned j = order[i];
                int index = Y * t + X0 + j;
                Ym[index] = mirror[i];
            }
        }
    }
    for (int index = 0; index < t*t; index++) {
        lut[index] = (Ym[index] << 4) | Xm[index];
    }
}

const char *USAGE_MESSAGE = "Usage: %s <offsets-file> <t>\n"
    "Offsets are lists of u's and v's separated by white space\n";

int main(int argc,char **argv) {
    if (argc != 3) {
        fprintf(stderr, USAGE_MESSAGE, argv[0]);
        exit(1);
    }
    std::string fileName = argv[1];
    t = atoi(argv[2]);
    u.resize(t*t);
    v.resize(t*t);
    Xm.resize(t*t);
    Ym.resize(t*t);
    lut.resize(t*t);
    FILE *file = fopen(fileName.c_str(), "r");
    for (int i = 0; i < t*t; i++) {
        fscanf(file, " %f %f", &u[i], &v[i]);
        if (feof(file)) {
            fprintf(stderr, "Error reading the reference file\n");
            exit(1);
        }
    }
    fclose(file);
    algorithm2();
    fileName.replace(fileName.length()-3, 3, "pgm");
    fprintf(stderr, "Saving table to %s\n", fileName.c_str());
    file = fopen(fileName.c_str(), "wb");
    int w, h, grades;
    fprintf(file, "P5 %d %d\n%d\n", t, t, 255);
    fwrite(lut.data(), 1, t*t, file);
    fclose(file);
}
