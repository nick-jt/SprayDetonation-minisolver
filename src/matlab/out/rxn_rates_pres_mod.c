#include <math.h>
#include "header.h"
#include "rates.h"

void get_rxn_pres_mod (const double T, const double pres, const double * __restrict__ C, double * __restrict__ pres_mod) {
  // third body variable declaration
  double thd;

  // pressure dependence variable declarations
  double k0;
  double kinf;
  double Pr;

  // troe variable declarations
  double logFcent;
  double A;
  double B;

  double logT = log(T);
  double m = pres / (8.31446210e+03 * T);

  // reaction 4;
  pres_mod[0] = m - 1.0 * C[4] - 1.0 * C[5] - 1.0 * C[16];

  // reaction 8;
  pres_mod[1] = m + 1.0 * C[4] + 5.3 * C[5] + 0.75 * C[15] + 2.6 * C[16];

  // reaction 9;
  pres_mod[2] = m + 1.0 * C[4] + 11.0 * C[5] + 0.75 * C[15] + 2.6 * C[16];

  // reaction 10;
  pres_mod[3] = m + 1.4 * C[4] + 14.4 * C[5] + 0.75 * C[15] + 2.6 * C[16];

  // reaction 11;
  thd = m - 0.15000000000000002 * C[7] + 10.89 * C[5] + 0.09000000000000008 * C[15] + 1.1800000000000002 * C[16];
  k0 = exp(3.1778590439387948e+01 - 1.4 * logT);
  kinf = exp(2.2355638720663009e+01 + 0.44 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.00000000e-01 * exp(-T / 1.00000000e-30) + 5.00000000e-01 * exp(-T / 1.00000000e+30), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[4] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 13;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 0.75 * C[15] + 2.6 * C[16];
  k0 = exp(2.6026570745005486e+01 - 0.584 * logT - (-1.1538824622220602e+03 / T));
  kinf = exp(2.5432796038258747e+01 - 0.37 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.65400000e-01 * exp(-T / 9.40000000e+01) + 7.34600000e-01 * exp(-T / 1.75600000e+03) + exp(-5.18200000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[5] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 29;
  thd = m + 1.0 * C[4] + 11.0 * C[5] + 0.75 * C[15] + 2.6 * C[16];
  k0 = exp(4.1606096243564160e+01 - 2.79 * logT - (2.1089931963247509e+03 / T));
  kinf = exp(1.6427049858685642e+01 - (1.1996754426242439e+03 / T));
  Pr = k0 * thd / kinf;
  pres_mod[6] =  Pr / (1.0 + Pr);

  // reaction 38;
  pres_mod[7] = m + 1.0 * C[4] - 1.0 * C[5] + 0.75 * C[15] + 2.6 * C[16];

  // reaction 41;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(4.9977627770478051e+01 - 3.42 * logT - (4.2446570295870377e+04 / T));
  kinf = exp(1.0668955394675699e+01 + 1.5 * logT - (4.0056277362789355e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.80000000e-02 * exp(-T / 1.97000000e+02) + 9.32000000e-01 * exp(-T / 1.54000000e+03) + exp(-1.03000000e+04 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[8] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 42;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(4.1746636266343160e+01 - 2.57 * logT - (7.1708787992430689e+02 / T));
  kinf = exp(2.0809443533187462e+01 + 0.48 * logT - (-1.3083708686338232e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.17600000e-01 * exp(-T / 2.71000000e+02) + 7.82400000e-01 * exp(-T / 2.75500000e+03) + exp(-6.57000000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[9] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 43;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(4.9517437762680643e+01 - 3.14 * logT - (6.1896006477677020e+02 / T));
  kinf = exp(3.0849896940796750e+01 - 0.8 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(3.20000000e-01 * exp(-T / 7.80000000e+01) + 6.80000000e-01 * exp(-T / 1.99500000e+03) + exp(-5.59000000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[10] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 62;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(5.6050499592221364e+01 - 4.8 * logT - (2.7979007806169448e+03 / T));
  kinf = exp(2.0107079697522593e+01 + 0.454 * logT - (1.3083708686338232e+03 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.42000000e-01 * exp(-T / 9.40000000e+01) + 7.58000000e-01 * exp(-T / 1.55500000e+03) + exp(-4.20000000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[11] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 68;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(6.3076845661346454e+01 - 4.76 * logT - (1.2278557382563572e+03 / T));
  kinf = exp(3.0172623109393093e+01 - 0.63 * logT - (1.9273309334105934e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.17000000e-01 * exp(-T / 7.40000000e+01) + 7.83000000e-01 * exp(-T / 2.94100000e+03) + exp(-6.96400000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[12] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 81;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(1.0188472363832375e+02 - 9.67 * logT - (3.1300256934239924e+03 / T));
  kinf = exp(3.0685022297606515e+01 - 0.97 * logT - (3.1199613021268090e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(4.67500000e-01 * exp(-T / 1.51000000e+02) + 5.32500000e-01 * exp(-T / 1.03800000e+03) + exp(-4.97000000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[13] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 94;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21] + 2.0 * C[17] + 2.0 * C[19];
  k0 = exp(5.6204000710479832e+01 - 3.4 * logT - (1.8014616300915008e+04 / T));
  kinf = exp(1.9771347927429105e+01 + 1.62 * logT - (1.8643379082815231e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(-9.81600000e-01 * exp(-T / 5.38370000e+03) + 1.98160000e+00 * exp(-T / 4.29320000e+00) + exp(7.95000000e-02 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[14] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 99;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21] + 2.0 * C[17] + 2.0 * C[19];
  k0 = exp(5.5598514468478307e+01 - 3.86 * logT - (1.6706889553324204e+03 / T));
  kinf = exp(2.2528270532924488e+01 + 0.27 * logT - (1.4090147816056557e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.18000000e-01 * exp(-T / 2.07500000e+02) + 7.82000000e-01 * exp(-T / 2.66300000e+03) + exp(-6.09500000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[15] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 111;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21] + 2.0 * C[17] + 2.0 * C[19];
  k0 = exp(1.2118603866293091e+02 - 11.94 * logT - (4.9163545047610478e+03 / T));
  kinf = 25000000000.0;
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(8.25000000e-01 * exp(-T / 1.34060000e+03) + 1.75000000e-01 * exp(-T / 6.00000000e+04) + exp(-1.01398000e+04 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[16] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 117;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(7.6691864936273376e+01 - 6.642 * logT - (2.9030736696725098e+03 / T));
  kinf = exp(1.4128129115706086e+01 + 1.463 * logT - (6.8186251038416549e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(1.56900000e+00 * exp(-T / 2.99000000e+02) + -5.69000000e-01 * exp(-T / 9.14700000e+03) + exp(1.52400000e+02 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[17] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 129;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(8.1278612893528006e+01 - 7.08 * logT - (3.3640227910835029e+03 / T));
  kinf = exp(3.3886771157681913e+01 - 0.99 * logT - (7.9508691247747720e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(1.57800000e-01 * exp(-T / 1.25000000e+02) + 8.42200000e-01 * exp(-T / 2.21900000e+03) + exp(-6.88200000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[18] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 137;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(1.1556750958063344e+02 - 11.79 * logT - (4.5211761804771477e+03 / T));
  kinf = 15000000000.0;
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(8.02000000e-01 * exp(-T / 2.27790000e+03) + 1.98000000e-01 * exp(-T / 6.00000000e+04) + exp(-5.72320000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[19] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 144;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(1.2462477396391213e+02 - 12.0 * logT - (3.0031137191665116e+03 / T));
  kinf = 200000000000.0;
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(9.80000000e-01 * exp(-T / 1.09660000e+03) + 2.00000000e-02 * exp(-T / 1.09660000e+03) + exp(-6.85950000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[20] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 151;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(1.2570313239567574e+02 - 12.81 * logT - (3.1451222803697674e+03 / T));
  kinf = exp(2.5328436022934504e+01 - 0.32 * logT - (-1.3199449186255839e+02 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(8.96000000e-01 * exp(-T / 1.60600000e+03) + 1.04000000e-01 * exp(-T / 6.00000000e+04) + exp(-6.11840000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[21] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 152;
  thd = m + 1.0 * C[4] + 5.0 * C[5] + 1.0 * C[11] + 0.5 * C[15] + 1.0 * C[16] + 2.0 * C[21];
  k0 = exp(7.5516903160921473e+01 - 6.66 * logT - (3.5225369540141392e+03 / T));
  kinf = exp(2.3311029872174121e+01 - (1.6408480351362718e+03 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(0.00000000e+00 * exp(-T / 1.00000000e+03) + 1.00000000e+00 * exp(-T / 1.31000000e+03) + exp(-4.80970000e+04 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[22] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

} // end get_rxn_pres_mod

