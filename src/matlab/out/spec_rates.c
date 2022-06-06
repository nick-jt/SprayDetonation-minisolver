#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 0
  sp_rates[0] = -(fwd_rates[0] - rev_rates[0]);
  //sp 1
  sp_rates[1] = (fwd_rates[0] - rev_rates[0]);
  //sp 2
  sp_rates[2] = (fwd_rates[0] - rev_rates[0]);
  //sp 7
  sp_rates[7] = -(fwd_rates[0] - rev_rates[0]);

  //rxn 1
  //sp 0
  sp_rates[0] += (fwd_rates[1] - rev_rates[1]);
  //sp 1
  sp_rates[1] -= (fwd_rates[1] - rev_rates[1]);
  //sp 2
  sp_rates[2] += (fwd_rates[1] - rev_rates[1]);
  //sp 4
  sp_rates[4] = -(fwd_rates[1] - rev_rates[1]);

  //rxn 2
  //sp 0
  sp_rates[0] += (fwd_rates[2] - rev_rates[2]);
  //sp 2
  sp_rates[2] -= (fwd_rates[2] - rev_rates[2]);
  //sp 4
  sp_rates[4] -= (fwd_rates[2] - rev_rates[2]);
  //sp 5
  sp_rates[5] = (fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 1
  sp_rates[1] += (fwd_rates[3] - rev_rates[3]);
  //sp 2
  sp_rates[2] -= 2.0 * (fwd_rates[3] - rev_rates[3]);
  //sp 5
  sp_rates[5] += (fwd_rates[3] - rev_rates[3]);

  //rxn 4
  //sp 0
  sp_rates[0] -= 2.0 * (fwd_rates[4] - rev_rates[4]) * pres_mod[0];
  //sp 4
  sp_rates[4] += (fwd_rates[4] - rev_rates[4]) * pres_mod[0];

  //rxn 5
  //sp 0
  sp_rates[0] -= 2.0 * (fwd_rates[5] - rev_rates[5]);
  //sp 4
  sp_rates[4] += (fwd_rates[5] - rev_rates[5]);

  //rxn 6
  //sp 0
  sp_rates[0] -= 2.0 * (fwd_rates[6] - rev_rates[6]);
  //sp 4
  sp_rates[4] += (fwd_rates[6] - rev_rates[6]);

  //rxn 7
  //sp 0
  sp_rates[0] -= 2.0 * (fwd_rates[7] - rev_rates[7]);
  //sp 4
  sp_rates[4] += (fwd_rates[7] - rev_rates[7]);

  //rxn 8
  //sp 0
  sp_rates[0] -= (fwd_rates[8] - rev_rates[8]) * pres_mod[1];
  //sp 2
  sp_rates[2] -= (fwd_rates[8] - rev_rates[8]) * pres_mod[1];
  //sp 5
  sp_rates[5] += (fwd_rates[8] - rev_rates[8]) * pres_mod[1];

  //rxn 9
  //sp 0
  sp_rates[0] -= (fwd_rates[9] - rev_rates[9]) * pres_mod[2];
  //sp 1
  sp_rates[1] -= (fwd_rates[9] - rev_rates[9]) * pres_mod[2];
  //sp 2
  sp_rates[2] += (fwd_rates[9] - rev_rates[9]) * pres_mod[2];

  //rxn 10
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[10] - rev_rates[10]) * pres_mod[3];
  //sp 7
  sp_rates[7] += (fwd_rates[10] - rev_rates[10]) * pres_mod[3];

  //rxn 11
  //sp 0
  sp_rates[0] -= (fwd_rates[11] - rev_rates[11]) * pres_mod[4];
  //sp 3
  sp_rates[3] = (fwd_rates[11] - rev_rates[11]) * pres_mod[4];
  //sp 7
  sp_rates[7] -= (fwd_rates[11] - rev_rates[11]) * pres_mod[4];

  //rxn 12
  //sp 0
  sp_rates[0] += (fwd_rates[12] - rev_rates[12]);
  //sp 3
  sp_rates[3] += (fwd_rates[12] - rev_rates[12]);
  //sp 4
  sp_rates[4] -= (fwd_rates[12] - rev_rates[12]);
  //sp 7
  sp_rates[7] -= (fwd_rates[12] - rev_rates[12]);

  //rxn 13
  //sp 2
  sp_rates[2] -= 2.0 * (fwd_rates[13] - rev_rates[13]) * pres_mod[5];
  //sp 6
  sp_rates[6] = (fwd_rates[13] - rev_rates[13]) * pres_mod[5];

  //rxn 14
  //sp 0
  sp_rates[0] -= (fwd_rates[14] - rev_rates[14]);
  //sp 1
  sp_rates[1] += (fwd_rates[14] - rev_rates[14]);
  //sp 3
  sp_rates[3] -= (fwd_rates[14] - rev_rates[14]);
  //sp 5
  sp_rates[5] += (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 0
  sp_rates[0] -= (fwd_rates[15] - rev_rates[15]);
  //sp 2
  sp_rates[2] += 2.0 * (fwd_rates[15] - rev_rates[15]);
  //sp 3
  sp_rates[3] -= (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 1
  sp_rates[1] -= (fwd_rates[16] - rev_rates[16]);
  //sp 2
  sp_rates[2] += (fwd_rates[16] - rev_rates[16]);
  //sp 3
  sp_rates[3] -= (fwd_rates[16] - rev_rates[16]);
  //sp 7
  sp_rates[7] += (fwd_rates[16] - rev_rates[16]);

  //rxn 17
  //sp 3
  sp_rates[3] -= 2.0 * (fwd_rates[17] - rev_rates[17]);
  //sp 6
  sp_rates[6] += (fwd_rates[17] - rev_rates[17]);
  //sp 7
  sp_rates[7] += (fwd_rates[17] - rev_rates[17]);

  //rxn 18
  //sp 3
  sp_rates[3] -= 2.0 * (fwd_rates[18] - rev_rates[18]);
  //sp 6
  sp_rates[6] += (fwd_rates[18] - rev_rates[18]);
  //sp 7
  sp_rates[7] += (fwd_rates[18] - rev_rates[18]);

  //rxn 19
  //sp 2
  sp_rates[2] -= (fwd_rates[19] - rev_rates[19]);
  //sp 3
  sp_rates[3] -= (fwd_rates[19] - rev_rates[19]);
  //sp 5
  sp_rates[5] += (fwd_rates[19] - rev_rates[19]);
  //sp 7
  sp_rates[7] += (fwd_rates[19] - rev_rates[19]);

  //rxn 20
  //sp 2
  sp_rates[2] -= (fwd_rates[20] - rev_rates[20]);
  //sp 3
  sp_rates[3] -= (fwd_rates[20] - rev_rates[20]);
  //sp 5
  sp_rates[5] += (fwd_rates[20] - rev_rates[20]);
  //sp 7
  sp_rates[7] += (fwd_rates[20] - rev_rates[20]);

  //rxn 21
  //sp 2
  sp_rates[2] -= (fwd_rates[21] - rev_rates[21]);
  //sp 3
  sp_rates[3] -= (fwd_rates[21] - rev_rates[21]);
  //sp 5
  sp_rates[5] += (fwd_rates[21] - rev_rates[21]);
  //sp 7
  sp_rates[7] += (fwd_rates[21] - rev_rates[21]);

  //rxn 22
  //sp 2
  sp_rates[2] -= (fwd_rates[22] - rev_rates[22]);
  //sp 3
  sp_rates[3] -= (fwd_rates[22] - rev_rates[22]);
  //sp 5
  sp_rates[5] += (fwd_rates[22] - rev_rates[22]);
  //sp 7
  sp_rates[7] += (fwd_rates[22] - rev_rates[22]);

  //rxn 23
  //sp 2
  sp_rates[2] -= (fwd_rates[23] - rev_rates[23]);
  //sp 3
  sp_rates[3] -= (fwd_rates[23] - rev_rates[23]);
  //sp 5
  sp_rates[5] += (fwd_rates[23] - rev_rates[23]);
  //sp 7
  sp_rates[7] += (fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 0
  sp_rates[0] -= (fwd_rates[24] - rev_rates[24]);
  //sp 3
  sp_rates[3] += (fwd_rates[24] - rev_rates[24]);
  //sp 4
  sp_rates[4] += (fwd_rates[24] - rev_rates[24]);
  //sp 6
  sp_rates[6] -= (fwd_rates[24] - rev_rates[24]);

  //rxn 25
  //sp 0
  sp_rates[0] -= (fwd_rates[25] - rev_rates[25]);
  //sp 2
  sp_rates[2] += (fwd_rates[25] - rev_rates[25]);
  //sp 5
  sp_rates[5] += (fwd_rates[25] - rev_rates[25]);
  //sp 6
  sp_rates[6] -= (fwd_rates[25] - rev_rates[25]);

  //rxn 26
  //sp 1
  sp_rates[1] -= (fwd_rates[26] - rev_rates[26]);
  //sp 2
  sp_rates[2] += (fwd_rates[26] - rev_rates[26]);
  //sp 3
  sp_rates[3] += (fwd_rates[26] - rev_rates[26]);
  //sp 6
  sp_rates[6] -= (fwd_rates[26] - rev_rates[26]);

  //rxn 27
  //sp 2
  sp_rates[2] -= (fwd_rates[27] - rev_rates[27]);
  //sp 3
  sp_rates[3] += (fwd_rates[27] - rev_rates[27]);
  //sp 5
  sp_rates[5] += (fwd_rates[27] - rev_rates[27]);
  //sp 6
  sp_rates[6] -= (fwd_rates[27] - rev_rates[27]);

  //rxn 28
  //sp 2
  sp_rates[2] -= (fwd_rates[28] - rev_rates[28]);
  //sp 3
  sp_rates[3] += (fwd_rates[28] - rev_rates[28]);
  //sp 5
  sp_rates[5] += (fwd_rates[28] - rev_rates[28]);
  //sp 6
  sp_rates[6] -= (fwd_rates[28] - rev_rates[28]);

  //rxn 29
  //sp 16
  sp_rates[16] = (fwd_rates[29] - rev_rates[29]) * pres_mod[6];
  //sp 1
  sp_rates[1] -= (fwd_rates[29] - rev_rates[29]) * pres_mod[6];
  //sp 15
  sp_rates[15] = -(fwd_rates[29] - rev_rates[29]) * pres_mod[6];

  //rxn 30
  //sp 16
  sp_rates[16] += (fwd_rates[30] - rev_rates[30]);
  //sp 0
  sp_rates[0] += (fwd_rates[30] - rev_rates[30]);
  //sp 2
  sp_rates[2] -= (fwd_rates[30] - rev_rates[30]);
  //sp 15
  sp_rates[15] -= (fwd_rates[30] - rev_rates[30]);

  //rxn 31
  //sp 16
  sp_rates[16] += (fwd_rates[31] - rev_rates[31]);
  //sp 0
  sp_rates[0] += (fwd_rates[31] - rev_rates[31]);
  //sp 2
  sp_rates[2] -= (fwd_rates[31] - rev_rates[31]);
  //sp 15
  sp_rates[15] -= (fwd_rates[31] - rev_rates[31]);

  //rxn 32
  //sp 16
  sp_rates[16] += (fwd_rates[32] - rev_rates[32]);
  //sp 1
  sp_rates[1] += (fwd_rates[32] - rev_rates[32]);
  //sp 7
  sp_rates[7] -= (fwd_rates[32] - rev_rates[32]);
  //sp 15
  sp_rates[15] -= (fwd_rates[32] - rev_rates[32]);

  //rxn 33
  //sp 16
  sp_rates[16] += (fwd_rates[33] - rev_rates[33]);
  //sp 2
  sp_rates[2] += (fwd_rates[33] - rev_rates[33]);
  //sp 3
  sp_rates[3] -= (fwd_rates[33] - rev_rates[33]);
  //sp 15
  sp_rates[15] -= (fwd_rates[33] - rev_rates[33]);

  //rxn 34
  //sp 0
  sp_rates[0] -= (fwd_rates[34] - rev_rates[34]);
  //sp 12
  sp_rates[12] = -(fwd_rates[34] - rev_rates[34]);
  //sp 4
  sp_rates[4] += (fwd_rates[34] - rev_rates[34]);
  //sp 15
  sp_rates[15] += (fwd_rates[34] - rev_rates[34]);

  //rxn 35
  //sp 1
  sp_rates[1] -= (fwd_rates[35] - rev_rates[35]);
  //sp 2
  sp_rates[2] += (fwd_rates[35] - rev_rates[35]);
  //sp 12
  sp_rates[12] -= (fwd_rates[35] - rev_rates[35]);
  //sp 15
  sp_rates[15] += (fwd_rates[35] - rev_rates[35]);

  //rxn 36
  //sp 16
  sp_rates[16] += (fwd_rates[36] - rev_rates[36]);
  //sp 1
  sp_rates[1] -= (fwd_rates[36] - rev_rates[36]);
  //sp 12
  sp_rates[12] -= (fwd_rates[36] - rev_rates[36]);
  //sp 0
  sp_rates[0] += (fwd_rates[36] - rev_rates[36]);

  //rxn 37
  //sp 2
  sp_rates[2] -= (fwd_rates[37] - rev_rates[37]);
  //sp 12
  sp_rates[12] -= (fwd_rates[37] - rev_rates[37]);
  //sp 5
  sp_rates[5] += (fwd_rates[37] - rev_rates[37]);
  //sp 15
  sp_rates[15] += (fwd_rates[37] - rev_rates[37]);

  //rxn 38
  //sp 0
  sp_rates[0] += (fwd_rates[38] - rev_rates[38]) * pres_mod[7];
  //sp 12
  sp_rates[12] -= (fwd_rates[38] - rev_rates[38]) * pres_mod[7];
  //sp 15
  sp_rates[15] += (fwd_rates[38] - rev_rates[38]) * pres_mod[7];

  //rxn 39
  //sp 0
  sp_rates[0] += (fwd_rates[39] - rev_rates[39]);
  //sp 12
  sp_rates[12] -= (fwd_rates[39] - rev_rates[39]);
  //sp 15
  sp_rates[15] += (fwd_rates[39] - rev_rates[39]);

  //rxn 40
  //sp 3
  sp_rates[3] += (fwd_rates[40] - rev_rates[40]);
  //sp 12
  sp_rates[12] -= (fwd_rates[40] - rev_rates[40]);
  //sp 15
  sp_rates[15] += (fwd_rates[40] - rev_rates[40]);
  //sp 7
  sp_rates[7] -= (fwd_rates[40] - rev_rates[40]);

  //rxn 41
  //sp 4
  sp_rates[4] -= (fwd_rates[41] - rev_rates[41]) * pres_mod[8];
  //sp 13
  sp_rates[13] = (fwd_rates[41] - rev_rates[41]) * pres_mod[8];
  //sp 15
  sp_rates[15] -= (fwd_rates[41] - rev_rates[41]) * pres_mod[8];

  //rxn 42
  //sp 0
  sp_rates[0] -= (fwd_rates[42] - rev_rates[42]) * pres_mod[9];
  //sp 12
  sp_rates[12] -= (fwd_rates[42] - rev_rates[42]) * pres_mod[9];
  //sp 13
  sp_rates[13] += (fwd_rates[42] - rev_rates[42]) * pres_mod[9];

  //rxn 43
  //sp 8
  sp_rates[8] = -(fwd_rates[43] - rev_rates[43]) * pres_mod[10];
  //sp 0
  sp_rates[0] -= (fwd_rates[43] - rev_rates[43]) * pres_mod[10];
  //sp 10
  sp_rates[10] = (fwd_rates[43] - rev_rates[43]) * pres_mod[10];

  //rxn 44
  //sp 8
  sp_rates[8] -= (fwd_rates[44] - rev_rates[44]);
  //sp 1
  sp_rates[1] -= (fwd_rates[44] - rev_rates[44]);
  //sp 12
  sp_rates[12] += (fwd_rates[44] - rev_rates[44]);
  //sp 0
  sp_rates[0] += (fwd_rates[44] - rev_rates[44]);

  //rxn 45
  //sp 8
  sp_rates[8] -= (fwd_rates[45] - rev_rates[45]);
  //sp 0
  sp_rates[0] += (fwd_rates[45] - rev_rates[45]);
  //sp 2
  sp_rates[2] -= (fwd_rates[45] - rev_rates[45]);
  //sp 13
  sp_rates[13] += (fwd_rates[45] - rev_rates[45]);

  //rxn 46
  //sp 8
  sp_rates[8] -= (fwd_rates[46] - rev_rates[46]);
  //sp 0
  sp_rates[0] += (fwd_rates[46] - rev_rates[46]);
  //sp 10
  sp_rates[10] += (fwd_rates[46] - rev_rates[46]);
  //sp 4
  sp_rates[4] -= (fwd_rates[46] - rev_rates[46]);

  //rxn 47
  //sp 8
  sp_rates[8] -= (fwd_rates[47] - rev_rates[47]);
  //sp 2
  sp_rates[2] += (fwd_rates[47] - rev_rates[47]);
  //sp 12
  sp_rates[12] += (fwd_rates[47] - rev_rates[47]);
  //sp 7
  sp_rates[7] -= (fwd_rates[47] - rev_rates[47]);

  //rxn 48
  //sp 8
  sp_rates[8] -= (fwd_rates[48] - rev_rates[48]);
  //sp 16
  sp_rates[16] += (fwd_rates[48] - rev_rates[48]);
  //sp 0
  sp_rates[0] += 2.0 * (fwd_rates[48] - rev_rates[48]);
  //sp 7
  sp_rates[7] -= (fwd_rates[48] - rev_rates[48]);

  //rxn 49
  //sp 8
  sp_rates[8] -= (fwd_rates[49] - rev_rates[49]);
  //sp 2
  sp_rates[2] += (fwd_rates[49] - rev_rates[49]);
  //sp 3
  sp_rates[3] -= (fwd_rates[49] - rev_rates[49]);
  //sp 13
  sp_rates[13] += (fwd_rates[49] - rev_rates[49]);

  //rxn 50
  //sp 8
  sp_rates[8] -= 2.0 * (fwd_rates[50] - rev_rates[50]);
  //sp 17
  sp_rates[17] = (fwd_rates[50] - rev_rates[50]);
  //sp 4
  sp_rates[4] += (fwd_rates[50] - rev_rates[50]);

  //rxn 51
  //sp 8
  sp_rates[8] += (fwd_rates[51] - rev_rates[51]);
  //sp 9
  sp_rates[9] = -(fwd_rates[51] - rev_rates[51]);

  //rxn 52
  //sp 9
  sp_rates[9] -= (fwd_rates[52] - rev_rates[52]);
  //sp 4
  sp_rates[4] += (fwd_rates[52] - rev_rates[52]);
  //sp 1
  sp_rates[1] -= (fwd_rates[52] - rev_rates[52]);
  //sp 15
  sp_rates[15] += (fwd_rates[52] - rev_rates[52]);

  //rxn 53
  //sp 0
  sp_rates[0] += (fwd_rates[53] - rev_rates[53]);
  //sp 9
  sp_rates[9] -= (fwd_rates[53] - rev_rates[53]);
  //sp 12
  sp_rates[12] += (fwd_rates[53] - rev_rates[53]);
  //sp 1
  sp_rates[1] -= (fwd_rates[53] - rev_rates[53]);

  //rxn 54
  //sp 0
  sp_rates[0] += (fwd_rates[54] - rev_rates[54]);
  //sp 9
  sp_rates[9] -= (fwd_rates[54] - rev_rates[54]);
  //sp 2
  sp_rates[2] -= (fwd_rates[54] - rev_rates[54]);
  //sp 13
  sp_rates[13] += (fwd_rates[54] - rev_rates[54]);

  //rxn 55
  //sp 0
  sp_rates[0] += (fwd_rates[55] - rev_rates[55]);
  //sp 9
  sp_rates[9] -= (fwd_rates[55] - rev_rates[55]);
  //sp 10
  sp_rates[10] += (fwd_rates[55] - rev_rates[55]);
  //sp 4
  sp_rates[4] -= (fwd_rates[55] - rev_rates[55]);

  //rxn 56
  //sp 0
  sp_rates[0] += (fwd_rates[56] - rev_rates[56]);
  //sp 2
  sp_rates[2] += (fwd_rates[56] - rev_rates[56]);
  //sp 7
  sp_rates[7] -= (fwd_rates[56] - rev_rates[56]);
  //sp 9
  sp_rates[9] -= (fwd_rates[56] - rev_rates[56]);
  //sp 15
  sp_rates[15] += (fwd_rates[56] - rev_rates[56]);

  //rxn 57
  //sp 9
  sp_rates[9] -= (fwd_rates[57] - rev_rates[57]);
  //sp 15
  sp_rates[15] += (fwd_rates[57] - rev_rates[57]);
  //sp 5
  sp_rates[5] += (fwd_rates[57] - rev_rates[57]);
  //sp 7
  sp_rates[7] -= (fwd_rates[57] - rev_rates[57]);

  //rxn 58
  //sp 8
  sp_rates[8] += (fwd_rates[58] - rev_rates[58]);
  //sp 9
  sp_rates[9] -= (fwd_rates[58] - rev_rates[58]);

  //rxn 59
  //sp 8
  sp_rates[8] += (fwd_rates[59] - rev_rates[59]);
  //sp 9
  sp_rates[9] -= (fwd_rates[59] - rev_rates[59]);

  //rxn 60
  //sp 9
  sp_rates[9] -= (fwd_rates[60] - rev_rates[60]);
  //sp 8
  sp_rates[8] += (fwd_rates[60] - rev_rates[60]);

  //rxn 61
  //sp 16
  sp_rates[16] -= (fwd_rates[61] - rev_rates[61]);
  //sp 9
  sp_rates[9] -= (fwd_rates[61] - rev_rates[61]);
  //sp 13
  sp_rates[13] += (fwd_rates[61] - rev_rates[61]);
  //sp 15
  sp_rates[15] += (fwd_rates[61] - rev_rates[61]);

  //rxn 62
  //sp 0
  sp_rates[0] -= (fwd_rates[62] - rev_rates[62]) * pres_mod[11];
  //sp 13
  sp_rates[13] -= (fwd_rates[62] - rev_rates[62]) * pres_mod[11];
  //sp 14
  sp_rates[14] = (fwd_rates[62] - rev_rates[62]) * pres_mod[11];

  //rxn 63
  //sp 0
  sp_rates[0] -= (fwd_rates[63] - rev_rates[63]);
  //sp 4
  sp_rates[4] += (fwd_rates[63] - rev_rates[63]);
  //sp 12
  sp_rates[12] += (fwd_rates[63] - rev_rates[63]);
  //sp 13
  sp_rates[13] -= (fwd_rates[63] - rev_rates[63]);

  //rxn 64
  //sp 1
  sp_rates[1] -= (fwd_rates[64] - rev_rates[64]);
  //sp 2
  sp_rates[2] += (fwd_rates[64] - rev_rates[64]);
  //sp 12
  sp_rates[12] += (fwd_rates[64] - rev_rates[64]);
  //sp 13
  sp_rates[13] -= (fwd_rates[64] - rev_rates[64]);

  //rxn 65
  //sp 2
  sp_rates[2] -= (fwd_rates[65] - rev_rates[65]);
  //sp 5
  sp_rates[5] += (fwd_rates[65] - rev_rates[65]);
  //sp 12
  sp_rates[12] += (fwd_rates[65] - rev_rates[65]);
  //sp 13
  sp_rates[13] -= (fwd_rates[65] - rev_rates[65]);

  //rxn 66
  //sp 3
  sp_rates[3] += (fwd_rates[66] - rev_rates[66]);
  //sp 12
  sp_rates[12] += (fwd_rates[66] - rev_rates[66]);
  //sp 13
  sp_rates[13] -= (fwd_rates[66] - rev_rates[66]);
  //sp 7
  sp_rates[7] -= (fwd_rates[66] - rev_rates[66]);

  //rxn 67
  //sp 3
  sp_rates[3] -= (fwd_rates[67] - rev_rates[67]);
  //sp 12
  sp_rates[12] += (fwd_rates[67] - rev_rates[67]);
  //sp 13
  sp_rates[13] -= (fwd_rates[67] - rev_rates[67]);
  //sp 6
  sp_rates[6] += (fwd_rates[67] - rev_rates[67]);

  //rxn 68
  //sp 0
  sp_rates[0] -= (fwd_rates[68] - rev_rates[68]) * pres_mod[12];
  //sp 10
  sp_rates[10] -= (fwd_rates[68] - rev_rates[68]) * pres_mod[12];
  //sp 11
  sp_rates[11] = (fwd_rates[68] - rev_rates[68]) * pres_mod[12];

  //rxn 69
  //sp 0
  sp_rates[0] += (fwd_rates[69] - rev_rates[69]);
  //sp 1
  sp_rates[1] -= (fwd_rates[69] - rev_rates[69]);
  //sp 10
  sp_rates[10] -= (fwd_rates[69] - rev_rates[69]);
  //sp 13
  sp_rates[13] += (fwd_rates[69] - rev_rates[69]);

  //rxn 70
  //sp 8
  sp_rates[8] += (fwd_rates[70] - rev_rates[70]);
  //sp 10
  sp_rates[10] -= (fwd_rates[70] - rev_rates[70]);
  //sp 2
  sp_rates[2] -= (fwd_rates[70] - rev_rates[70]);
  //sp 5
  sp_rates[5] += (fwd_rates[70] - rev_rates[70]);

  //rxn 71
  //sp 9
  sp_rates[9] += (fwd_rates[71] - rev_rates[71]);
  //sp 10
  sp_rates[10] -= (fwd_rates[71] - rev_rates[71]);
  //sp 2
  sp_rates[2] -= (fwd_rates[71] - rev_rates[71]);
  //sp 5
  sp_rates[5] += (fwd_rates[71] - rev_rates[71]);

  //rxn 72
  //sp 1
  sp_rates[1] += (fwd_rates[72] - rev_rates[72]);
  //sp 10
  sp_rates[10] -= (fwd_rates[72] - rev_rates[72]);
  //sp 14
  sp_rates[14] += (fwd_rates[72] - rev_rates[72]);
  //sp 7
  sp_rates[7] -= (fwd_rates[72] - rev_rates[72]);

  //rxn 73
  //sp 10
  sp_rates[10] -= (fwd_rates[73] - rev_rates[73]);
  //sp 2
  sp_rates[2] += (fwd_rates[73] - rev_rates[73]);
  //sp 13
  sp_rates[13] += (fwd_rates[73] - rev_rates[73]);
  //sp 7
  sp_rates[7] -= (fwd_rates[73] - rev_rates[73]);

  //rxn 74
  //sp 11
  sp_rates[11] += (fwd_rates[74] - rev_rates[74]);
  //sp 10
  sp_rates[10] -= (fwd_rates[74] - rev_rates[74]);
  //sp 3
  sp_rates[3] -= (fwd_rates[74] - rev_rates[74]);
  //sp 7
  sp_rates[7] += (fwd_rates[74] - rev_rates[74]);

  //rxn 75
  //sp 2
  sp_rates[2] += (fwd_rates[75] - rev_rates[75]);
  //sp 10
  sp_rates[10] -= (fwd_rates[75] - rev_rates[75]);
  //sp 3
  sp_rates[3] -= (fwd_rates[75] - rev_rates[75]);
  //sp 14
  sp_rates[14] += (fwd_rates[75] - rev_rates[75]);

  //rxn 76
  //sp 3
  sp_rates[3] += (fwd_rates[76] - rev_rates[76]);
  //sp 10
  sp_rates[10] -= (fwd_rates[76] - rev_rates[76]);
  //sp 11
  sp_rates[11] += (fwd_rates[76] - rev_rates[76]);
  //sp 6
  sp_rates[6] -= (fwd_rates[76] - rev_rates[76]);

  //rxn 77
  //sp 10
  sp_rates[10] -= (fwd_rates[77] - rev_rates[77]);
  //sp 11
  sp_rates[11] += (fwd_rates[77] - rev_rates[77]);
  //sp 12
  sp_rates[12] -= (fwd_rates[77] - rev_rates[77]);
  //sp 15
  sp_rates[15] += (fwd_rates[77] - rev_rates[77]);

  //rxn 78
  //sp 10
  sp_rates[10] -= (fwd_rates[78] - rev_rates[78]);
  //sp 11
  sp_rates[11] += (fwd_rates[78] - rev_rates[78]);
  //sp 12
  sp_rates[12] += (fwd_rates[78] - rev_rates[78]);
  //sp 13
  sp_rates[13] -= (fwd_rates[78] - rev_rates[78]);

  //rxn 79
  //sp 8
  sp_rates[8] -= (fwd_rates[79] - rev_rates[79]);
  //sp 0
  sp_rates[0] += (fwd_rates[79] - rev_rates[79]);
  //sp 10
  sp_rates[10] -= (fwd_rates[79] - rev_rates[79]);
  //sp 19
  sp_rates[19] = (fwd_rates[79] - rev_rates[79]);

  //rxn 80
  //sp 0
  sp_rates[0] += (fwd_rates[80] - rev_rates[80]);
  //sp 9
  sp_rates[9] -= (fwd_rates[80] - rev_rates[80]);
  //sp 10
  sp_rates[10] -= (fwd_rates[80] - rev_rates[80]);
  //sp 19
  sp_rates[19] += (fwd_rates[80] - rev_rates[80]);

  //rxn 81
  //sp 10
  sp_rates[10] -= 2.0 * (fwd_rates[81] - rev_rates[81]) * pres_mod[13];
  //sp 21
  sp_rates[21] = (fwd_rates[81] - rev_rates[81]) * pres_mod[13];

  //rxn 82
  //sp 0
  sp_rates[0] += (fwd_rates[82] - rev_rates[82]);
  //sp 10
  sp_rates[10] -= 2.0 * (fwd_rates[82] - rev_rates[82]);
  //sp 20
  sp_rates[20] = (fwd_rates[82] - rev_rates[82]);

  //rxn 83
  //sp 0
  sp_rates[0] -= (fwd_rates[83] - rev_rates[83]);
  //sp 4
  sp_rates[4] += (fwd_rates[83] - rev_rates[83]);
  //sp 13
  sp_rates[13] += (fwd_rates[83] - rev_rates[83]);
  //sp 14
  sp_rates[14] -= (fwd_rates[83] - rev_rates[83]);

  //rxn 84
  //sp 0
  sp_rates[0] -= (fwd_rates[84] - rev_rates[84]);
  //sp 10
  sp_rates[10] += (fwd_rates[84] - rev_rates[84]);
  //sp 2
  sp_rates[2] += (fwd_rates[84] - rev_rates[84]);
  //sp 14
  sp_rates[14] -= (fwd_rates[84] - rev_rates[84]);

  //rxn 85
  //sp 0
  sp_rates[0] -= (fwd_rates[85] - rev_rates[85]);
  //sp 9
  sp_rates[9] += (fwd_rates[85] - rev_rates[85]);
  //sp 5
  sp_rates[5] += (fwd_rates[85] - rev_rates[85]);
  //sp 14
  sp_rates[14] -= (fwd_rates[85] - rev_rates[85]);

  //rxn 86
  //sp 1
  sp_rates[1] -= (fwd_rates[86] - rev_rates[86]);
  //sp 2
  sp_rates[2] += (fwd_rates[86] - rev_rates[86]);
  //sp 13
  sp_rates[13] += (fwd_rates[86] - rev_rates[86]);
  //sp 14
  sp_rates[14] -= (fwd_rates[86] - rev_rates[86]);

  //rxn 87
  //sp 2
  sp_rates[2] -= (fwd_rates[87] - rev_rates[87]);
  //sp 5
  sp_rates[5] += (fwd_rates[87] - rev_rates[87]);
  //sp 13
  sp_rates[13] += (fwd_rates[87] - rev_rates[87]);
  //sp 14
  sp_rates[14] -= (fwd_rates[87] - rev_rates[87]);

  //rxn 88
  //sp 3
  sp_rates[3] += (fwd_rates[88] - rev_rates[88]);
  //sp 13
  sp_rates[13] += (fwd_rates[88] - rev_rates[88]);
  //sp 14
  sp_rates[14] -= (fwd_rates[88] - rev_rates[88]);
  //sp 7
  sp_rates[7] -= (fwd_rates[88] - rev_rates[88]);

  //rxn 89
  //sp 0
  sp_rates[0] -= (fwd_rates[89] - rev_rates[89]);
  //sp 10
  sp_rates[10] += (fwd_rates[89] - rev_rates[89]);
  //sp 11
  sp_rates[11] -= (fwd_rates[89] - rev_rates[89]);
  //sp 4
  sp_rates[4] += (fwd_rates[89] - rev_rates[89]);

  //rxn 90
  //sp 2
  sp_rates[2] += (fwd_rates[90] - rev_rates[90]);
  //sp 1
  sp_rates[1] -= (fwd_rates[90] - rev_rates[90]);
  //sp 10
  sp_rates[10] += (fwd_rates[90] - rev_rates[90]);
  //sp 11
  sp_rates[11] -= (fwd_rates[90] - rev_rates[90]);

  //rxn 91
  //sp 10
  sp_rates[10] += (fwd_rates[91] - rev_rates[91]);
  //sp 2
  sp_rates[2] -= (fwd_rates[91] - rev_rates[91]);
  //sp 11
  sp_rates[11] -= (fwd_rates[91] - rev_rates[91]);
  //sp 5
  sp_rates[5] += (fwd_rates[91] - rev_rates[91]);

  //rxn 92
  //sp 8
  sp_rates[8] -= (fwd_rates[92] - rev_rates[92]);
  //sp 10
  sp_rates[10] += 2.0 * (fwd_rates[92] - rev_rates[92]);
  //sp 11
  sp_rates[11] -= (fwd_rates[92] - rev_rates[92]);

  //rxn 93
  //sp 9
  sp_rates[9] -= (fwd_rates[93] - rev_rates[93]);
  //sp 10
  sp_rates[10] += 2.0 * (fwd_rates[93] - rev_rates[93]);
  //sp 11
  sp_rates[11] -= (fwd_rates[93] - rev_rates[93]);

  //rxn 94
  //sp 0
  sp_rates[0] += (fwd_rates[94] - rev_rates[94]) * pres_mod[14];
  //sp 17
  sp_rates[17] += (fwd_rates[94] - rev_rates[94]) * pres_mod[14];
  //sp 18
  sp_rates[18] = -(fwd_rates[94] - rev_rates[94]) * pres_mod[14];

  //rxn 95
  //sp 8
  sp_rates[8] += (fwd_rates[95] - rev_rates[95]);
  //sp 17
  sp_rates[17] -= (fwd_rates[95] - rev_rates[95]);
  //sp 1
  sp_rates[1] -= (fwd_rates[95] - rev_rates[95]);
  //sp 15
  sp_rates[15] += (fwd_rates[95] - rev_rates[95]);

  //rxn 96
  //sp 17
  sp_rates[17] -= (fwd_rates[96] - rev_rates[96]);
  //sp 2
  sp_rates[2] -= (fwd_rates[96] - rev_rates[96]);
  //sp 10
  sp_rates[10] += (fwd_rates[96] - rev_rates[96]);
  //sp 15
  sp_rates[15] += (fwd_rates[96] - rev_rates[96]);

  //rxn 97
  //sp 17
  sp_rates[17] -= (fwd_rates[97] - rev_rates[97]);
  //sp 18
  sp_rates[18] += (fwd_rates[97] - rev_rates[97]);
  //sp 12
  sp_rates[12] -= (fwd_rates[97] - rev_rates[97]);
  //sp 15
  sp_rates[15] += (fwd_rates[97] - rev_rates[97]);

  //rxn 98
  //sp 17
  sp_rates[17] -= (fwd_rates[98] - rev_rates[98]);
  //sp 10
  sp_rates[10] -= (fwd_rates[98] - rev_rates[98]);
  //sp 23
  sp_rates[23] = (fwd_rates[98] - rev_rates[98]);

  //rxn 99
  //sp 0
  sp_rates[0] -= (fwd_rates[99] - rev_rates[99]) * pres_mod[15];
  //sp 18
  sp_rates[18] -= (fwd_rates[99] - rev_rates[99]) * pres_mod[15];
  //sp 19
  sp_rates[19] += (fwd_rates[99] - rev_rates[99]) * pres_mod[15];

  //rxn 100
  //sp 0
  sp_rates[0] -= (fwd_rates[100] - rev_rates[100]);
  //sp 17
  sp_rates[17] += (fwd_rates[100] - rev_rates[100]);
  //sp 18
  sp_rates[18] -= (fwd_rates[100] - rev_rates[100]);
  //sp 4
  sp_rates[4] += (fwd_rates[100] - rev_rates[100]);

  //rxn 101
  //sp 1
  sp_rates[1] -= (fwd_rates[101] - rev_rates[101]);
  //sp 18
  sp_rates[18] -= (fwd_rates[101] - rev_rates[101]);
  //sp 10
  sp_rates[10] += (fwd_rates[101] - rev_rates[101]);
  //sp 15
  sp_rates[15] += (fwd_rates[101] - rev_rates[101]);

  //rxn 102
  //sp 17
  sp_rates[17] += (fwd_rates[102] - rev_rates[102]);
  //sp 18
  sp_rates[18] -= (fwd_rates[102] - rev_rates[102]);
  //sp 2
  sp_rates[2] -= (fwd_rates[102] - rev_rates[102]);
  //sp 5
  sp_rates[5] += (fwd_rates[102] - rev_rates[102]);

  //rxn 103
  //sp 17
  sp_rates[17] += (fwd_rates[103] - rev_rates[103]);
  //sp 18
  sp_rates[18] -= (fwd_rates[103] - rev_rates[103]);
  //sp 3
  sp_rates[3] += (fwd_rates[103] - rev_rates[103]);
  //sp 7
  sp_rates[7] -= (fwd_rates[103] - rev_rates[103]);

  //rxn 104
  //sp 1
  sp_rates[1] += (fwd_rates[104] - rev_rates[104]);
  //sp 18
  sp_rates[18] -= (fwd_rates[104] - rev_rates[104]);
  //sp 22
  sp_rates[22] = (fwd_rates[104] - rev_rates[104]);
  //sp 7
  sp_rates[7] -= (fwd_rates[104] - rev_rates[104]);

  //rxn 105
  //sp 18
  sp_rates[18] -= (fwd_rates[105] - rev_rates[105]);
  //sp 12
  sp_rates[12] += (fwd_rates[105] - rev_rates[105]);
  //sp 13
  sp_rates[13] += (fwd_rates[105] - rev_rates[105]);
  //sp 7
  sp_rates[7] -= (fwd_rates[105] - rev_rates[105]);

  //rxn 106
  //sp 2
  sp_rates[2] += (fwd_rates[106] - rev_rates[106]);
  //sp 18
  sp_rates[18] -= (fwd_rates[106] - rev_rates[106]);
  //sp 3
  sp_rates[3] -= (fwd_rates[106] - rev_rates[106]);
  //sp 22
  sp_rates[22] += (fwd_rates[106] - rev_rates[106]);

  //rxn 107
  //sp 3
  sp_rates[3] += (fwd_rates[107] - rev_rates[107]);
  //sp 18
  sp_rates[18] -= (fwd_rates[107] - rev_rates[107]);
  //sp 19
  sp_rates[19] += (fwd_rates[107] - rev_rates[107]);
  //sp 6
  sp_rates[6] -= (fwd_rates[107] - rev_rates[107]);

  //rxn 108
  //sp 18
  sp_rates[18] -= (fwd_rates[108] - rev_rates[108]);
  //sp 19
  sp_rates[19] += (fwd_rates[108] - rev_rates[108]);
  //sp 12
  sp_rates[12] -= (fwd_rates[108] - rev_rates[108]);
  //sp 15
  sp_rates[15] += (fwd_rates[108] - rev_rates[108]);

  //rxn 109
  //sp 18
  sp_rates[18] -= (fwd_rates[109] - rev_rates[109]);
  //sp 26
  sp_rates[26] = (fwd_rates[109] - rev_rates[109]);
  //sp 12
  sp_rates[12] -= (fwd_rates[109] - rev_rates[109]);

  //rxn 110
  //sp 11
  sp_rates[11] += (fwd_rates[110] - rev_rates[110]);
  //sp 17
  sp_rates[17] += (fwd_rates[110] - rev_rates[110]);
  //sp 18
  sp_rates[18] -= (fwd_rates[110] - rev_rates[110]);
  //sp 10
  sp_rates[10] -= (fwd_rates[110] - rev_rates[110]);

  //rxn 111
  //sp 24
  sp_rates[24] = (fwd_rates[111] - rev_rates[111]) * pres_mod[16];
  //sp 18
  sp_rates[18] -= (fwd_rates[111] - rev_rates[111]) * pres_mod[16];
  //sp 10
  sp_rates[10] -= (fwd_rates[111] - rev_rates[111]) * pres_mod[16];

  //rxn 112
  //sp 0
  sp_rates[0] += (fwd_rates[112] - rev_rates[112]);
  //sp 18
  sp_rates[18] -= (fwd_rates[112] - rev_rates[112]);
  //sp 10
  sp_rates[10] -= (fwd_rates[112] - rev_rates[112]);
  //sp 23
  sp_rates[23] += (fwd_rates[112] - rev_rates[112]);

  //rxn 113
  //sp 17
  sp_rates[17] += (fwd_rates[113] - rev_rates[113]);
  //sp 18
  sp_rates[18] -= 2.0 * (fwd_rates[113] - rev_rates[113]);
  //sp 19
  sp_rates[19] += (fwd_rates[113] - rev_rates[113]);

  //rxn 114
  //sp 10
  sp_rates[10] += (fwd_rates[114] - rev_rates[114]);
  //sp 22
  sp_rates[22] -= (fwd_rates[114] - rev_rates[114]);
  //sp 15
  sp_rates[15] += (fwd_rates[114] - rev_rates[114]);

  //rxn 115
  //sp 0
  sp_rates[0] -= (fwd_rates[115] - rev_rates[115]);
  //sp 10
  sp_rates[10] += (fwd_rates[115] - rev_rates[115]);
  //sp 12
  sp_rates[12] += (fwd_rates[115] - rev_rates[115]);
  //sp 22
  sp_rates[22] -= (fwd_rates[115] - rev_rates[115]);

  //rxn 116
  //sp 2
  sp_rates[2] += (fwd_rates[116] - rev_rates[116]);
  //sp 7
  sp_rates[7] -= (fwd_rates[116] - rev_rates[116]);
  //sp 13
  sp_rates[13] += (fwd_rates[116] - rev_rates[116]);
  //sp 15
  sp_rates[15] += (fwd_rates[116] - rev_rates[116]);
  //sp 22
  sp_rates[22] -= (fwd_rates[116] - rev_rates[116]);

  //rxn 117
  //sp 0
  sp_rates[0] -= (fwd_rates[117] - rev_rates[117]) * pres_mod[17];
  //sp 19
  sp_rates[19] -= (fwd_rates[117] - rev_rates[117]) * pres_mod[17];
  //sp 20
  sp_rates[20] += (fwd_rates[117] - rev_rates[117]) * pres_mod[17];

  //rxn 118
  //sp 0
  sp_rates[0] -= (fwd_rates[118] - rev_rates[118]);
  //sp 18
  sp_rates[18] += (fwd_rates[118] - rev_rates[118]);
  //sp 19
  sp_rates[19] -= (fwd_rates[118] - rev_rates[118]);
  //sp 4
  sp_rates[4] += (fwd_rates[118] - rev_rates[118]);

  //rxn 119
  //sp 2
  sp_rates[2] += (fwd_rates[119] - rev_rates[119]);
  //sp 1
  sp_rates[1] -= (fwd_rates[119] - rev_rates[119]);
  //sp 18
  sp_rates[18] += (fwd_rates[119] - rev_rates[119]);
  //sp 19
  sp_rates[19] -= (fwd_rates[119] - rev_rates[119]);

  //rxn 120
  //sp 1
  sp_rates[1] -= (fwd_rates[120] - rev_rates[120]);
  //sp 10
  sp_rates[10] += (fwd_rates[120] - rev_rates[120]);
  //sp 19
  sp_rates[19] -= (fwd_rates[120] - rev_rates[120]);
  //sp 12
  sp_rates[12] += (fwd_rates[120] - rev_rates[120]);

  //rxn 121
  //sp 8
  sp_rates[8] += (fwd_rates[121] - rev_rates[121]);
  //sp 1
  sp_rates[1] -= (fwd_rates[121] - rev_rates[121]);
  //sp 19
  sp_rates[19] -= (fwd_rates[121] - rev_rates[121]);
  //sp 13
  sp_rates[13] += (fwd_rates[121] - rev_rates[121]);

  //rxn 122
  //sp 18
  sp_rates[18] += (fwd_rates[122] - rev_rates[122]);
  //sp 2
  sp_rates[2] -= (fwd_rates[122] - rev_rates[122]);
  //sp 19
  sp_rates[19] -= (fwd_rates[122] - rev_rates[122]);
  //sp 5
  sp_rates[5] += (fwd_rates[122] - rev_rates[122]);

  //rxn 123
  //sp 19
  sp_rates[19] -= (fwd_rates[123] - rev_rates[123]);
  //sp 12
  sp_rates[12] -= (fwd_rates[123] - rev_rates[123]);
  //sp 20
  sp_rates[20] += (fwd_rates[123] - rev_rates[123]);
  //sp 15
  sp_rates[15] += (fwd_rates[123] - rev_rates[123]);

  //rxn 124
  //sp 8
  sp_rates[8] -= (fwd_rates[124] - rev_rates[124]);
  //sp 0
  sp_rates[0] += (fwd_rates[124] - rev_rates[124]);
  //sp 19
  sp_rates[19] -= (fwd_rates[124] - rev_rates[124]);
  //sp 23
  sp_rates[23] += (fwd_rates[124] - rev_rates[124]);

  //rxn 125
  //sp 0
  sp_rates[0] += (fwd_rates[125] - rev_rates[125]);
  //sp 9
  sp_rates[9] -= (fwd_rates[125] - rev_rates[125]);
  //sp 19
  sp_rates[19] -= (fwd_rates[125] - rev_rates[125]);
  //sp 23
  sp_rates[23] += (fwd_rates[125] - rev_rates[125]);

  //rxn 126
  //sp 18
  sp_rates[18] += (fwd_rates[126] - rev_rates[126]);
  //sp 11
  sp_rates[11] += (fwd_rates[126] - rev_rates[126]);
  //sp 10
  sp_rates[10] -= (fwd_rates[126] - rev_rates[126]);
  //sp 19
  sp_rates[19] -= (fwd_rates[126] - rev_rates[126]);

  //rxn 127
  //sp 25
  sp_rates[25] = (fwd_rates[127] - rev_rates[127]);
  //sp 10
  sp_rates[10] -= (fwd_rates[127] - rev_rates[127]);
  //sp 19
  sp_rates[19] -= (fwd_rates[127] - rev_rates[127]);

  //rxn 128
  //sp 3
  sp_rates[3] += (fwd_rates[128] - rev_rates[128]);
  //sp 18
  sp_rates[18] += (fwd_rates[128] - rev_rates[128]);
  //sp 19
  sp_rates[19] -= (fwd_rates[128] - rev_rates[128]);
  //sp 7
  sp_rates[7] -= (fwd_rates[128] - rev_rates[128]);

  //rxn 129
  //sp 0
  sp_rates[0] -= (fwd_rates[129] - rev_rates[129]) * pres_mod[18];
  //sp 20
  sp_rates[20] -= (fwd_rates[129] - rev_rates[129]) * pres_mod[18];
  //sp 21
  sp_rates[21] += (fwd_rates[129] - rev_rates[129]) * pres_mod[18];

  //rxn 130
  //sp 0
  sp_rates[0] -= (fwd_rates[130] - rev_rates[130]);
  //sp 19
  sp_rates[19] += (fwd_rates[130] - rev_rates[130]);
  //sp 20
  sp_rates[20] -= (fwd_rates[130] - rev_rates[130]);
  //sp 4
  sp_rates[4] += (fwd_rates[130] - rev_rates[130]);

  //rxn 131
  //sp 1
  sp_rates[1] -= (fwd_rates[131] - rev_rates[131]);
  //sp 10
  sp_rates[10] += (fwd_rates[131] - rev_rates[131]);
  //sp 20
  sp_rates[20] -= (fwd_rates[131] - rev_rates[131]);
  //sp 13
  sp_rates[13] += (fwd_rates[131] - rev_rates[131]);

  //rxn 132
  //sp 3
  sp_rates[3] += (fwd_rates[132] - rev_rates[132]);
  //sp 19
  sp_rates[19] += (fwd_rates[132] - rev_rates[132]);
  //sp 20
  sp_rates[20] -= (fwd_rates[132] - rev_rates[132]);
  //sp 7
  sp_rates[7] -= (fwd_rates[132] - rev_rates[132]);

  //rxn 133
  //sp 3
  sp_rates[3] -= (fwd_rates[133] - rev_rates[133]);
  //sp 20
  sp_rates[20] -= (fwd_rates[133] - rev_rates[133]);
  //sp 21
  sp_rates[21] += (fwd_rates[133] - rev_rates[133]);
  //sp 7
  sp_rates[7] += (fwd_rates[133] - rev_rates[133]);

  //rxn 134
  //sp 19
  sp_rates[19] += (fwd_rates[134] - rev_rates[134]);
  //sp 3
  sp_rates[3] -= (fwd_rates[134] - rev_rates[134]);
  //sp 20
  sp_rates[20] -= (fwd_rates[134] - rev_rates[134]);
  //sp 6
  sp_rates[6] += (fwd_rates[134] - rev_rates[134]);

  //rxn 135
  //sp 2
  sp_rates[2] += (fwd_rates[135] - rev_rates[135]);
  //sp 3
  sp_rates[3] -= (fwd_rates[135] - rev_rates[135]);
  //sp 10
  sp_rates[10] += (fwd_rates[135] - rev_rates[135]);
  //sp 13
  sp_rates[13] += (fwd_rates[135] - rev_rates[135]);
  //sp 20
  sp_rates[20] -= (fwd_rates[135] - rev_rates[135]);

  //rxn 136
  //sp 3
  sp_rates[3] += (fwd_rates[136] - rev_rates[136]);
  //sp 20
  sp_rates[20] -= (fwd_rates[136] - rev_rates[136]);
  //sp 21
  sp_rates[21] += (fwd_rates[136] - rev_rates[136]);
  //sp 6
  sp_rates[6] -= (fwd_rates[136] - rev_rates[136]);

  //rxn 137
  //sp 18
  sp_rates[18] -= (fwd_rates[137] - rev_rates[137]) * pres_mod[19];
  //sp 27
  sp_rates[27] = (fwd_rates[137] - rev_rates[137]) * pres_mod[19];
  //sp 20
  sp_rates[20] -= (fwd_rates[137] - rev_rates[137]) * pres_mod[19];

  //rxn 138
  //sp 18
  sp_rates[18] -= (fwd_rates[138] - rev_rates[138]);
  //sp 10
  sp_rates[10] += (fwd_rates[138] - rev_rates[138]);
  //sp 20
  sp_rates[20] -= (fwd_rates[138] - rev_rates[138]);
  //sp 23
  sp_rates[23] += (fwd_rates[138] - rev_rates[138]);

  //rxn 139
  //sp 0
  sp_rates[0] -= (fwd_rates[139] - rev_rates[139]);
  //sp 4
  sp_rates[4] += (fwd_rates[139] - rev_rates[139]);
  //sp 20
  sp_rates[20] += (fwd_rates[139] - rev_rates[139]);
  //sp 21
  sp_rates[21] -= (fwd_rates[139] - rev_rates[139]);

  //rxn 140
  //sp 1
  sp_rates[1] -= (fwd_rates[140] - rev_rates[140]);
  //sp 2
  sp_rates[2] += (fwd_rates[140] - rev_rates[140]);
  //sp 20
  sp_rates[20] += (fwd_rates[140] - rev_rates[140]);
  //sp 21
  sp_rates[21] -= (fwd_rates[140] - rev_rates[140]);

  //rxn 141
  //sp 2
  sp_rates[2] -= (fwd_rates[141] - rev_rates[141]);
  //sp 5
  sp_rates[5] += (fwd_rates[141] - rev_rates[141]);
  //sp 20
  sp_rates[20] += (fwd_rates[141] - rev_rates[141]);
  //sp 21
  sp_rates[21] -= (fwd_rates[141] - rev_rates[141]);

  //rxn 142
  //sp 9
  sp_rates[9] -= (fwd_rates[142] - rev_rates[142]);
  //sp 10
  sp_rates[10] += (fwd_rates[142] - rev_rates[142]);
  //sp 20
  sp_rates[20] += (fwd_rates[142] - rev_rates[142]);
  //sp 21
  sp_rates[21] -= (fwd_rates[142] - rev_rates[142]);

  //rxn 143
  //sp 10
  sp_rates[10] -= (fwd_rates[143] - rev_rates[143]);
  //sp 11
  sp_rates[11] += (fwd_rates[143] - rev_rates[143]);
  //sp 20
  sp_rates[20] += (fwd_rates[143] - rev_rates[143]);
  //sp 21
  sp_rates[21] -= (fwd_rates[143] - rev_rates[143]);

  //rxn 144
  //sp 0
  sp_rates[0] -= (fwd_rates[144] - rev_rates[144]) * pres_mod[20];
  //sp 24
  sp_rates[24] += (fwd_rates[144] - rev_rates[144]) * pres_mod[20];
  //sp 23
  sp_rates[23] -= (fwd_rates[144] - rev_rates[144]) * pres_mod[20];

  //rxn 145
  //sp 0
  sp_rates[0] += (fwd_rates[145] - rev_rates[145]);
  //sp 1
  sp_rates[1] -= (fwd_rates[145] - rev_rates[145]);
  //sp 26
  sp_rates[26] += (fwd_rates[145] - rev_rates[145]);
  //sp 23
  sp_rates[23] -= (fwd_rates[145] - rev_rates[145]);

  //rxn 146
  //sp 0
  sp_rates[0] += 2.0 * (fwd_rates[146] - rev_rates[146]);
  //sp 2
  sp_rates[2] -= (fwd_rates[146] - rev_rates[146]);
  //sp 26
  sp_rates[26] += (fwd_rates[146] - rev_rates[146]);
  //sp 23
  sp_rates[23] -= (fwd_rates[146] - rev_rates[146]);

  //rxn 147
  //sp 26
  sp_rates[26] += (fwd_rates[147] - rev_rates[147]);
  //sp 2
  sp_rates[2] += (fwd_rates[147] - rev_rates[147]);
  //sp 7
  sp_rates[7] -= (fwd_rates[147] - rev_rates[147]);
  //sp 23
  sp_rates[23] -= (fwd_rates[147] - rev_rates[147]);

  //rxn 148
  //sp 24
  sp_rates[24] += (fwd_rates[148] - rev_rates[148]);
  //sp 3
  sp_rates[3] -= (fwd_rates[148] - rev_rates[148]);
  //sp 7
  sp_rates[7] += (fwd_rates[148] - rev_rates[148]);
  //sp 23
  sp_rates[23] -= (fwd_rates[148] - rev_rates[148]);

  //rxn 149
  //sp 2
  sp_rates[2] += (fwd_rates[149] - rev_rates[149]);
  //sp 3
  sp_rates[3] -= (fwd_rates[149] - rev_rates[149]);
  //sp 13
  sp_rates[13] += (fwd_rates[149] - rev_rates[149]);
  //sp 18
  sp_rates[18] += (fwd_rates[149] - rev_rates[149]);
  //sp 23
  sp_rates[23] -= (fwd_rates[149] - rev_rates[149]);

  //rxn 150
  //sp 24
  sp_rates[24] += (fwd_rates[150] - rev_rates[150]);
  //sp 12
  sp_rates[12] -= (fwd_rates[150] - rev_rates[150]);
  //sp 15
  sp_rates[15] += (fwd_rates[150] - rev_rates[150]);
  //sp 23
  sp_rates[23] -= (fwd_rates[150] - rev_rates[150]);

  //rxn 151
  //sp 10
  sp_rates[10] -= (fwd_rates[151] - rev_rates[151]) * pres_mod[21];
  //sp 27
  sp_rates[27] += (fwd_rates[151] - rev_rates[151]) * pres_mod[21];
  //sp 23
  sp_rates[23] -= (fwd_rates[151] - rev_rates[151]) * pres_mod[21];

  //rxn 152
  //sp 24
  sp_rates[24] -= (fwd_rates[152] - rev_rates[152]) * pres_mod[22];
  //sp 0
  sp_rates[0] -= (fwd_rates[152] - rev_rates[152]) * pres_mod[22];
  //sp 25
  sp_rates[25] += (fwd_rates[152] - rev_rates[152]) * pres_mod[22];

  //rxn 153
  //sp 24
  sp_rates[24] -= (fwd_rates[153] - rev_rates[153]);
  //sp 0
  sp_rates[0] -= (fwd_rates[153] - rev_rates[153]);
  //sp 10
  sp_rates[10] += (fwd_rates[153] - rev_rates[153]);
  //sp 19
  sp_rates[19] += (fwd_rates[153] - rev_rates[153]);

  //rxn 154
  //sp 24
  sp_rates[24] -= (fwd_rates[154] - rev_rates[154]);
  //sp 0
  sp_rates[0] -= (fwd_rates[154] - rev_rates[154]);
  //sp 4
  sp_rates[4] += (fwd_rates[154] - rev_rates[154]);
  //sp 23
  sp_rates[23] += (fwd_rates[154] - rev_rates[154]);

  //rxn 155
  //sp 24
  sp_rates[24] -= (fwd_rates[155] - rev_rates[155]);
  //sp 1
  sp_rates[1] -= (fwd_rates[155] - rev_rates[155]);
  //sp 26
  sp_rates[26] += (fwd_rates[155] - rev_rates[155]);
  //sp 0
  sp_rates[0] += 2.0 * (fwd_rates[155] - rev_rates[155]);

  //rxn 156
  //sp 24
  sp_rates[24] -= (fwd_rates[156] - rev_rates[156]);
  //sp 1
  sp_rates[1] -= (fwd_rates[156] - rev_rates[156]);
  //sp 20
  sp_rates[20] += (fwd_rates[156] - rev_rates[156]);
  //sp 12
  sp_rates[12] += (fwd_rates[156] - rev_rates[156]);

  //rxn 157
  //sp 24
  sp_rates[24] -= (fwd_rates[157] - rev_rates[157]);
  //sp 1
  sp_rates[1] -= (fwd_rates[157] - rev_rates[157]);
  //sp 2
  sp_rates[2] += (fwd_rates[157] - rev_rates[157]);
  //sp 23
  sp_rates[23] += (fwd_rates[157] - rev_rates[157]);

  //rxn 158
  //sp 24
  sp_rates[24] -= (fwd_rates[158] - rev_rates[158]);
  //sp 2
  sp_rates[2] -= (fwd_rates[158] - rev_rates[158]);
  //sp 5
  sp_rates[5] += (fwd_rates[158] - rev_rates[158]);
  //sp 23
  sp_rates[23] += (fwd_rates[158] - rev_rates[158]);

  //rxn 159
  //sp 24
  sp_rates[24] -= (fwd_rates[159] - rev_rates[159]);
  //sp 3
  sp_rates[3] -= (fwd_rates[159] - rev_rates[159]);
  //sp 6
  sp_rates[6] += (fwd_rates[159] - rev_rates[159]);
  //sp 23
  sp_rates[23] += (fwd_rates[159] - rev_rates[159]);

  //rxn 160
  //sp 24
  sp_rates[24] -= (fwd_rates[160] - rev_rates[160]);
  //sp 10
  sp_rates[10] -= (fwd_rates[160] - rev_rates[160]);
  //sp 11
  sp_rates[11] += (fwd_rates[160] - rev_rates[160]);
  //sp 23
  sp_rates[23] += (fwd_rates[160] - rev_rates[160]);

  //rxn 161
  //sp 0
  sp_rates[0] -= (fwd_rates[161] - rev_rates[161]);
  //sp 26
  sp_rates[26] -= (fwd_rates[161] - rev_rates[161]);
  //sp 19
  sp_rates[19] += (fwd_rates[161] - rev_rates[161]);
  //sp 12
  sp_rates[12] += (fwd_rates[161] - rev_rates[161]);

  //rxn 162
  //sp 1
  sp_rates[1] -= (fwd_rates[162] - rev_rates[162]);
  //sp 2
  sp_rates[2] += (fwd_rates[162] - rev_rates[162]);
  //sp 15
  sp_rates[15] += (fwd_rates[162] - rev_rates[162]);
  //sp 18
  sp_rates[18] += (fwd_rates[162] - rev_rates[162]);
  //sp 26
  sp_rates[26] -= (fwd_rates[162] - rev_rates[162]);

  //rxn 163
  //sp 2
  sp_rates[2] -= (fwd_rates[163] - rev_rates[163]);
  //sp 5
  sp_rates[5] += (fwd_rates[163] - rev_rates[163]);
  //sp 15
  sp_rates[15] += (fwd_rates[163] - rev_rates[163]);
  //sp 18
  sp_rates[18] += (fwd_rates[163] - rev_rates[163]);
  //sp 26
  sp_rates[26] -= (fwd_rates[163] - rev_rates[163]);

  //rxn 164
  //sp 0
  sp_rates[0] -= (fwd_rates[164] - rev_rates[164]);
  //sp 25
  sp_rates[25] -= (fwd_rates[164] - rev_rates[164]);
  //sp 10
  sp_rates[10] += (fwd_rates[164] - rev_rates[164]);
  //sp 20
  sp_rates[20] += (fwd_rates[164] - rev_rates[164]);

  //rxn 165
  //sp 0
  sp_rates[0] -= (fwd_rates[165] - rev_rates[165]);
  //sp 25
  sp_rates[25] -= (fwd_rates[165] - rev_rates[165]);
  //sp 4
  sp_rates[4] += (fwd_rates[165] - rev_rates[165]);
  //sp 24
  sp_rates[24] += (fwd_rates[165] - rev_rates[165]);

  //rxn 166
  //sp 25
  sp_rates[25] -= (fwd_rates[166] - rev_rates[166]);
  //sp 20
  sp_rates[20] += (fwd_rates[166] - rev_rates[166]);
  //sp 13
  sp_rates[13] += (fwd_rates[166] - rev_rates[166]);
  //sp 1
  sp_rates[1] -= (fwd_rates[166] - rev_rates[166]);

  //rxn 167
  //sp 24
  sp_rates[24] += (fwd_rates[167] - rev_rates[167]);
  //sp 25
  sp_rates[25] -= (fwd_rates[167] - rev_rates[167]);
  //sp 2
  sp_rates[2] -= (fwd_rates[167] - rev_rates[167]);
  //sp 5
  sp_rates[5] += (fwd_rates[167] - rev_rates[167]);

  //rxn 168
  //sp 24
  sp_rates[24] += (fwd_rates[168] - rev_rates[168]);
  //sp 25
  sp_rates[25] -= (fwd_rates[168] - rev_rates[168]);
  //sp 3
  sp_rates[3] += (fwd_rates[168] - rev_rates[168]);
  //sp 7
  sp_rates[7] -= (fwd_rates[168] - rev_rates[168]);

  //rxn 169
  //sp 2
  sp_rates[2] += (fwd_rates[169] - rev_rates[169]);
  //sp 3
  sp_rates[3] -= (fwd_rates[169] - rev_rates[169]);
  //sp 13
  sp_rates[13] += (fwd_rates[169] - rev_rates[169]);
  //sp 20
  sp_rates[20] += (fwd_rates[169] - rev_rates[169]);
  //sp 25
  sp_rates[25] -= (fwd_rates[169] - rev_rates[169]);

  //rxn 170
  //sp 24
  sp_rates[24] += (fwd_rates[170] - rev_rates[170]);
  //sp 25
  sp_rates[25] -= (fwd_rates[170] - rev_rates[170]);
  //sp 10
  sp_rates[10] -= (fwd_rates[170] - rev_rates[170]);
  //sp 11
  sp_rates[11] += (fwd_rates[170] - rev_rates[170]);

  //rxn 171
  //sp 0
  sp_rates[0] -= (fwd_rates[171] - rev_rates[171]);
  //sp 19
  sp_rates[19] += (fwd_rates[171] - rev_rates[171]);
  //sp 27
  sp_rates[27] -= (fwd_rates[171] - rev_rates[171]);
  //sp 20
  sp_rates[20] += (fwd_rates[171] - rev_rates[171]);

  //rxn 172
  //sp 0
  sp_rates[0] -= (fwd_rates[172] - rev_rates[172]);
  //sp 24
  sp_rates[24] += (fwd_rates[172] - rev_rates[172]);
  //sp 10
  sp_rates[10] += (fwd_rates[172] - rev_rates[172]);
  //sp 27
  sp_rates[27] -= (fwd_rates[172] - rev_rates[172]);

  //rxn 173
  //sp 1
  sp_rates[1] -= (fwd_rates[173] - rev_rates[173]);
  //sp 27
  sp_rates[27] -= (fwd_rates[173] - rev_rates[173]);
  //sp 12
  sp_rates[12] += (fwd_rates[173] - rev_rates[173]);
  //sp 25
  sp_rates[25] += (fwd_rates[173] - rev_rates[173]);

  //rxn 174
  //sp 10
  sp_rates[10] += fwd_rates[174];
  //sp 19
  sp_rates[19] += fwd_rates[174];
  //sp 25
  sp_rates[25] += fwd_rates[174];
  //sp 28
  sp_rates[28] = -fwd_rates[174];
  //sp 29
  sp_rates[29] = fwd_rates[174];

  //rxn 175
  //sp 25
  sp_rates[25] += 2.0 * fwd_rates[175];
  //sp 28
  sp_rates[28] -= fwd_rates[175];
  //sp 29
  sp_rates[29] += fwd_rates[175];

  //rxn 176
  //sp 19
  sp_rates[19] += fwd_rates[176];
  //sp 25
  sp_rates[25] += 2.0 * fwd_rates[176];
  //sp 27
  sp_rates[27] += fwd_rates[176];
  //sp 28
  sp_rates[28] -= fwd_rates[176];

  //rxn 177
  //sp 24
  sp_rates[24] += 2.0 * fwd_rates[177];
  //sp 25
  sp_rates[25] += 2.0 * fwd_rates[177];
  //sp 28
  sp_rates[28] -= fwd_rates[177];

  //rxn 178
  //sp 0
  sp_rates[0] -= fwd_rates[178];
  //sp 4
  sp_rates[4] += fwd_rates[178];
  //sp 24
  sp_rates[24] += fwd_rates[178];
  //sp 25
  sp_rates[25] += fwd_rates[178];
  //sp 28
  sp_rates[28] -= fwd_rates[178];
  //sp 29
  sp_rates[29] += fwd_rates[178];

  //rxn 179
  //sp 0
  sp_rates[0] -= fwd_rates[179];
  //sp 4
  sp_rates[4] += fwd_rates[179];
  //sp 24
  sp_rates[24] += fwd_rates[179];
  //sp 25
  sp_rates[25] += fwd_rates[179];
  //sp 28
  sp_rates[28] -= fwd_rates[179];
  //sp 29
  sp_rates[29] += fwd_rates[179];

  //rxn 180
  //sp 1
  sp_rates[1] -= fwd_rates[180];
  //sp 2
  sp_rates[2] += fwd_rates[180];
  //sp 24
  sp_rates[24] += fwd_rates[180];
  //sp 25
  sp_rates[25] += fwd_rates[180];
  //sp 28
  sp_rates[28] -= fwd_rates[180];
  //sp 29
  sp_rates[29] += fwd_rates[180];

  //rxn 181
  //sp 1
  sp_rates[1] -= fwd_rates[181];
  //sp 2
  sp_rates[2] += fwd_rates[181];
  //sp 24
  sp_rates[24] += fwd_rates[181];
  //sp 25
  sp_rates[25] += fwd_rates[181];
  //sp 28
  sp_rates[28] -= fwd_rates[181];
  //sp 29
  sp_rates[29] += fwd_rates[181];

  //rxn 182
  //sp 2
  sp_rates[2] -= fwd_rates[182];
  //sp 5
  sp_rates[5] += fwd_rates[182];
  //sp 24
  sp_rates[24] += fwd_rates[182];
  //sp 25
  sp_rates[25] += fwd_rates[182];
  //sp 28
  sp_rates[28] -= fwd_rates[182];
  //sp 29
  sp_rates[29] += fwd_rates[182];

  //rxn 183
  //sp 2
  sp_rates[2] -= fwd_rates[183];
  //sp 5
  sp_rates[5] += fwd_rates[183];
  //sp 24
  sp_rates[24] += fwd_rates[183];
  //sp 25
  sp_rates[25] += fwd_rates[183];
  //sp 28
  sp_rates[28] -= fwd_rates[183];
  //sp 29
  sp_rates[29] += fwd_rates[183];

  //rxn 184
  //sp 10
  sp_rates[10] -= fwd_rates[184];
  //sp 11
  sp_rates[11] += fwd_rates[184];
  //sp 24
  sp_rates[24] += fwd_rates[184];
  //sp 25
  sp_rates[25] += fwd_rates[184];
  //sp 28
  sp_rates[28] -= fwd_rates[184];
  //sp 29
  sp_rates[29] += fwd_rates[184];

  //rxn 185
  //sp 10
  sp_rates[10] -= fwd_rates[185];
  //sp 11
  sp_rates[11] += fwd_rates[185];
  //sp 24
  sp_rates[24] += fwd_rates[185];
  //sp 25
  sp_rates[25] += fwd_rates[185];
  //sp 28
  sp_rates[28] -= fwd_rates[185];
  //sp 29
  sp_rates[29] += fwd_rates[185];

  //rxn 186
  //sp 3
  sp_rates[3] += fwd_rates[186];
  //sp 7
  sp_rates[7] -= fwd_rates[186];
  //sp 24
  sp_rates[24] += fwd_rates[186];
  //sp 25
  sp_rates[25] += fwd_rates[186];
  //sp 28
  sp_rates[28] -= fwd_rates[186];
  //sp 29
  sp_rates[29] += fwd_rates[186];

  //rxn 187
  //sp 3
  sp_rates[3] += fwd_rates[187];
  //sp 7
  sp_rates[7] -= fwd_rates[187];
  //sp 24
  sp_rates[24] += fwd_rates[187];
  //sp 25
  sp_rates[25] += fwd_rates[187];
  //sp 28
  sp_rates[28] -= fwd_rates[187];
  //sp 29
  sp_rates[29] += fwd_rates[187];

  //rxn 188
  //sp 3
  sp_rates[3] -= fwd_rates[188];
  //sp 6
  sp_rates[6] += fwd_rates[188];
  //sp 24
  sp_rates[24] += fwd_rates[188];
  //sp 25
  sp_rates[25] += fwd_rates[188];
  //sp 28
  sp_rates[28] -= fwd_rates[188];
  //sp 29
  sp_rates[29] += fwd_rates[188];

  //rxn 189
  //sp 3
  sp_rates[3] -= fwd_rates[189];
  //sp 6
  sp_rates[6] += fwd_rates[189];
  //sp 24
  sp_rates[24] += fwd_rates[189];
  //sp 25
  sp_rates[25] += fwd_rates[189];
  //sp 28
  sp_rates[28] -= fwd_rates[189];
  //sp 29
  sp_rates[29] += fwd_rates[189];

  //rxn 190
  //sp 25
  sp_rates[25] += (fwd_rates[190] - rev_rates[174]);
  //sp 29
  sp_rates[29] -= (fwd_rates[190] - rev_rates[174]);
  //sp 23
  sp_rates[23] += (fwd_rates[190] - rev_rates[174]);

  //rxn 191
  //sp 0
  sp_rates[0] -= (fwd_rates[191] - rev_rates[175]);
  //sp 24
  sp_rates[24] += (fwd_rates[191] - rev_rates[175]);
  //sp 29
  sp_rates[29] -= (fwd_rates[191] - rev_rates[175]);
  //sp 25
  sp_rates[25] += (fwd_rates[191] - rev_rates[175]);

  //rxn 192
  //sp 1
  sp_rates[1] -= (fwd_rates[192] - rev_rates[176]);
  //sp 12
  sp_rates[12] += (fwd_rates[192] - rev_rates[176]);
  //sp 19
  sp_rates[19] += (fwd_rates[192] - rev_rates[176]);
  //sp 25
  sp_rates[25] += (fwd_rates[192] - rev_rates[176]);
  //sp 29
  sp_rates[29] -= (fwd_rates[192] - rev_rates[176]);

  //sp 30
  (*dy_N) = 0.0;
} // end eval_spec_rates

