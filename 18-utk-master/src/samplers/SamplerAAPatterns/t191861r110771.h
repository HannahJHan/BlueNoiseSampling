/* A header file to define pattern parameters and dependent constants
 * for use with dttp.cpp
 * Note that T00 and T63 never need to be used.
 * Created by ?? ??
 * 2014-11-11
 */
#ifndef T191861R110771_H
#define T191861R110771_H

const unsigned long aa_t = 191861;
const unsigned long aa_r = 110771;
const unsigned long MAP_WIDTH = 63;

enum {                                                                          // These threshold are generated by algorithm 1 in the paper (see thresholds.cpp)
    T00 = 0,
    T01 = 2131,
    T02 = 5822,
    T03 = 7953,
    T04 = 13775,
    T05 = 15906,
    T06 = 19597,
    T07 = 21728,
    T08 = 23859,
    T09 = 27550,
    T10 = 29681,
    T11 = 31812,
    T12 = 35503,
    T13 = 37634,
    T14 = 43456,
    T15 = 45587,
    T16 = 49278,
    T17 = 51409,
    T18 = 53540,
    T19 = 57231,
    T20 = 59362,
    T21 = 61493,
    T22 = 65184,
    T23 = 67315,
    T24 = 73137,
    T25 = 75268,
    T26 = 78959,
    T27 = 81090,
    T28 = 83221,
    T29 = 86912,
    T30 = 89043,
    T31 = 94865,
    T32 = 96996,
    T33 = 102818,
    T34 = 104949,
    T35 = 108640,
    T36 = 110771,
    T37 = 112902,
    T38 = 116593,
    T39 = 118724,
    T40 = 124546,
    T41 = 126677,
    T42 = 130368,
    T43 = 132499,
    T44 = 134630,
    T45 = 138321,
    T46 = 140452,
    T47 = 142583,
    T48 = 146274,
    T49 = 148405,
    T50 = 154227,
    T51 = 156358,
    T52 = 160049,
    T53 = 162180,
    T54 = 164311,
    T55 = 168002,
    T56 = 170133,
    T57 = 172264,
    T58 = 175955,
    T59 = 178086,
    T60 = 183908,
    T61 = 186039,
    T62 = 189730,
    T63 = 191861
};

/* Binary search set-region indexes in vectors look-up table.
 * We hard-code thresholds for better performance: Yes, the algorithm is so
 * efficient that saving these few machine cycles gives substantial speedup
 * We put this in header file because the depth depends on the map, which in
 * turn depends on t and r.
 */
inline int bsearch (int x) {
    if (x < T32)
        if (x < T16)
            if (x < T08)
                if (x < T04)
                    if (x < T02)
                        if (x < T01)
                            return 0;
                        else
                            return 1;
                    else
                        if (x < T03)
                            return 2;
                        else
                            return 3;
                else
                    if (x < T06)
                        if (x < T05)
                            return 4;
                        else
                            return 5;
                    else
                        if (x < T07)
                            return 6;
                        else
                            return 7;
            else
                if (x < T12)
                    if (x < T10)
                        if (x < T09)
                            return 8;
                        else
                            return 9;
                    else
                        if (x < T11)
                            return 10;
                        else
                            return 11;
                else
                    if (x < T14)
                        if (x < T13)
                            return 12;
                        else
                            return 13;
                    else
                        if (x < T15)
                            return 14;
                        else
                            return 15;
        else
            if (x < T24)
                if (x < T20)
                    if (x < T18)
                        if (x < T17)
                            return 16;
                        else
                            return 17;
                    else
                        if (x < T19)
                            return 18;
                        else
                            return 19;
                else
                    if (x < T22)
                        if (x < T21)
                            return 20;
                        else
                            return 21;
                    else
                        if (x < T23)
                            return 22;
                        else
                            return 23;
            else
                if (x < T28)
                    if (x < T26)
                        if (x < T25)
                            return 24;
                        else
                            return 25;
                    else
                        if (x < T27)
                            return 26;
                        else
                            return 27;
                else
                    if (x < T30)
                        if (x < T29)
                            return 28;
                        else
                            return 29;
                    else
                        if (x < T31)
                            return 30;
                        else
                            return 31;
    else
        if (x < T48)
            if (x < T40)
                if (x < T36)
                    if (x < T34)
                        if (x < T33)
                            return 32;
                        else
                            return 33;
                    else
                        if (x < T35)
                            return 34;
                        else
                            return 35;
                else
                    if (x < T38)
                        if (x < T37)
                            return 36;
                        else
                            return 37;
                    else
                        if (x < T39)
                            return 38;
                        else
                            return 39;
            else
                if (x < T44)
                    if (x < T42)
                        if (x < T41)
                            return 40;
                        else
                            return 41;
                    else
                        if (x < T43)
                            return 42;
                        else
                            return 43;
                else
                    if (x < T46)
                        if (x < T45)
                            return 44;
                        else
                            return 45;
                    else
                        if (x < T47)
                            return 46;
                        else
                            return 47;
        else
            if (x < T56)
                if (x < T52)
                    if (x < T50)
                        if (x < T49)
                            return 48;
                        else
                            return 49;
                    else
                        if (x < T51)
                            return 50;
                        else
                            return 51;
                else
                    if (x < T54)
                        if (x < T53)
                            return 52;
                        else
                            return 53;
                    else
                        if (x < T55)
                            return 54;
                        else
                            return 55;
            else
                if (x < T60)
                    if (x < T58)
                        if (x < T57)
                            return 56;
                        else
                            return 57;
                    else
                        if (x < T59)
                            return 58;
                        else
                            return 59;
                else
                    if (x < T62)
                        if (x < T61)
                            return 60;
                        else
                            return 61;
                    else
                        if (x < T63)
                            return 62;
                        else
                            return 63;
}


#endif
