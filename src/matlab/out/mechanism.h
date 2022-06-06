#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 30
/* Species Indexes
0  H
1  O
2  OH
3  HO2
4  H2
5  H2O
6  H2O2
7  O2
8  CH2
9  CH2*
10  CH3
11  CH4
12  HCO
13  CH2O
14  CH3O
15  CO
16  CO2
17  C2H2
18  C2H3
19  C2H4
20  C2H5
21  C2H6
22  CH2CHO
23  aC3H5
24  C3H6
25  nC3H7
26  C2H3CHO
27  C4H81
28  NC12H26
29  C6H12
30  N2
*/

//Number of species
#define NSP 31
//Number of variables. NN = NSP + 1 (temperature)
#define NN 32
//Number of forward reactions
#define FWD_RATES 193
//Number of reversible reactions
#define REV_RATES 177
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 23

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

