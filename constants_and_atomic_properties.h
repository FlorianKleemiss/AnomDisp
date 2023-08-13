#pragma once
#include "convenience.h"

const double prefactor = pow(el_charge,2)/(2*pow(PI,2)*el_mass);
const double r_0 = pow(el_charge,2)/el_charge/pow(speed_of_light,2);
const double r_0_for_nu = pow(el_charge,2)/el_charge;
const double constant_factor = 4*PI2*el_mass/h/h; //p_0 according to Hönl's Waller paper (1933) Page 646

double q(const double& nu){
  return 2*PI*nu/speed_of_light;
}

cdouble n_prime_from_z(const double& z, const double& n_0){
  return n_0/std::sqrt(cdouble(z)-1.0);
}

cdouble exp_from_z(const double& z, const double& n_0){
  const cdouble z_1 = z-1;
  cdouble temp = 0.0;
  const cdouble sqrt_var = std::sqrt(abs(z_1));
  if(z<=1)
    temp -= 2.0/3.0*z_1*hypergeometric(0.75,1.0,1.75,z_1*z_1);
  temp += std::atan(sqrt_var)/sqrt_var;
  temp = std::exp(-2.0*n_0*temp);
  return temp;
}

double get_Zeff_1s(const int& Z){
  //Get the shielding constant according to the data in https://research.unl.pt/ws/portalfiles/portal/13073849/ADNDT_2016_7_Revision_2.pdf
  const double Zeff[] = {    1.0,                                                                                                                                                                                                                  1.618,
           2.618, 3.614,                                                                                                                                                                         4.603, 5.587, 6.570, 7.553, 8.534, 9.515,
          10.498,11.484,                                                                                                                                                                        12.471,13.459,14.448,15.437,16.427,17.417,
          18.408,19.400,                                                                                                  20.394,21.388,22.382,23.378,24.372,25.368,26.364,27.360,28.357,29.353,30.349,31.345,32.341,33.337,34.333,35.329,
          36.325,37.322,                                                                                                  38.318,39.315,40.311,41.308,42.308,43.302,44.299,45.296,46.294,47.291,48.288,49.285,50.282,51.280,52.277,53.274,
          54.271,55.269,56.266,57.264,58.262,59.259,60.257,61.255,62.252,63.250,64.248,65.246,66.244,67.241,68.239,69.237,70.234,71.232,72.230,73.227,74.225,75.222,76.220,77.217,78.214,79.211,80.208,81.204,82.200,83.197,84.192,85.186,
          86.182,87.176,88.171,89.165,90.158,91.149,92.143,93.134,94.123,95.115,96.105,97.093,98.081,99.067,100.051,101.034,102.017,102.997,103.976,104.941,105.914,106.884,107.851,108.815,109.775,110.731,111.684,112.631,113.574,114.511,115.442,116.367
  };
  if (Z<=118) return Zeff[Z-1];
  else exit(-1);
}

double get_Zeff_2s(const int& Z){
  //Get the shielding constant according to the data in https://research.unl.pt/ws/portalfiles/portal/13073849/ADNDT_2016_7_Revision_2.pdf
  const double Zeff[] = {   0.0,                                                                                                                                                                                                                      0.0,
            1.550, 2.266,                                                                                                                                                                         3.036, 3.776, 4.505, 5.257, 5.997, 6.731,
            7.708, 8.699,                                                                                                                                                                         9.686,10.669,11.648,12.625,13.599,14.572,
           15.550,16.530,                                                                                                  17.510,18.490,19.468,20.446,21.421,22.397,23.373,24.348,25.324,26.298,27.273,28.248,29.224,30.201,31.179,32.158,
           33.139,34.121,                                                                                                  35.105,36.089,37.074,38.060,39.048,40.036,41.025,42.014,43.005,43.997,44.989,45.982,46.976,47.971,48.967,49.963,
           50.960,51.958,52.957,53.960,54.966,55.971,56.976,57.982,58.989,59.994,61.004,62.013,63.023,64.033,65.044,66.055,67.065,68.076,69.087,70.099,71.111,72.123,73.136,74.149,75.163,76.177,77.192,78.207,79.222,80.238,81.254,82.270,
           83.287,84.304,85.322,86.340,87.359,88.376,89.397,90.416,91.435,92.455,93.476,94.495,95.515,96.534,97.552,98.571,99.589,100.607,101.623,102.631,103.644,104.656,105.666,106.674,107.680,108.684,109.684,110.682,111.676,112.666,113.651,114.631
  };
  if (Z<=118) return Zeff[Z-1];
  else exit(-1);
}

double get_Zeff_2p_1_2(const int& Z){
  //Get the shielding constant according to the data in https://research.unl.pt/ws/portalfiles/portal/13073849/ADNDT_2016_7_Revision_2.pdf
  const double Zeff[] = {    0.0,                                                                                                                                                                                                                     0.0,
              0.0,   0.0,                                                                                                                                                                         2.268, 2.916, 3.547, 4.022, 4.598, 5.185,
            6.268, 7.306,                                                                                                                                                                         8.343, 9.362,10.355,11.356,12.357,13.340,
           14.330,15.320,                                                                                                  16.303,17.286,18.270,19.263,20.224,21.218,22.194,23.172,24.161,25.137,26.114,27.093,28.071,29.051,30.033,31.014,
           31.998,32.983,                                                                                                  33.968,34.955,35.943,36.932,37.922,38.912,39.903,40.897,41.890,42.884,43.879,44.875,45.872,46.869,47.868,48.867,
           49.867,50.868,51.870,52.876,53.886,54.894,55.903,56.913,57.924,58.933,59.947,60.960,61.973,62.988,64.003,65.019,66.034,67.049,68.065,69.082,70.099,71.116,72.135,73.154,74.173,75.193,76.214,77.236,78.258,79.281,80.304,81.329,
           82.354,83.380,84.406,85.434,86.463,87.492,88.523,89.554,90.586,91.618,92.652,93.687,94.722,95.757,96.794,97.831,98.868,99.906,100.945,101.982,103.021,104.061,105.101,106.141,107.181,108.220,109.260,110.299,111.337,112.374,113.410,114.444
  };
  if (Z<=118) return Zeff[Z-1];
  else exit(-1);
}

double get_Zeff_2p_3_2(const int& Z){
  //Get the shielding constant according to the data in https://research.unl.pt/ws/portalfiles/portal/13073849/ADNDT_2016_7_Revision_2.pdf
  const double Zeff[] = {    0.0,                                                                                                                                                                                                                     0.0,
              0.0,   0.0,                                                                                                                                                                         2.268, 2.916, 3.546, 4.082, 4.615, 5.175,
            6.256, 7.292,                                                                                                                                                                         8.313, 9.323,10.332,11.321,12.306,13.304,
           14.289,15.275,                                                                                                  16.262,17.244,18.222,19.193,20.165,21.140,22.112,23.081,24.045,25.011,25.976,26.942,27.909,28.876,29.844,30.813,
           31.782,32.752,                                                                                                  33.723,34.695,35.666,36.637,37.610,38.582,39.555,40.528,41.501,42.474,43.447,44.421,45.395,46.370,47.345,48.319,
           49.295,50.270,51.245,52.225,53.209,54.189,55.170,56.151,57.132,58.111,59.095,60.077,61.058,62.040,63.021,64.003,64.983,65.962,66.941,67.919,68.897,69.875,70.853,71.830,72.807,73.784,74.761,75.737,76.713,77.689,78.665,79.640,
           80.615,81.590,82.565,83.539,84.515,85.489,86.464,87.439,88.413,89.386,90.361,91.334,92.307,93.279,94.251,95.223,96.194,97.164,98.134,99.103,100.072,101.040,102.007,102.974,103.940,104.905,105.870,106.834,107.797,108.759,109.721,110.682
  };
  if (Z<=118) return Zeff[Z-1];
  else exit(-1);
}

double get_Zeff_3s(const int& Z){
  //Get the shielding constant according to the data in https://research.unl.pt/ws/portalfiles/portal/13073849/ADNDT_2016_7_Revision_2.pdf
  const double Zeff[] = {    0.0,                                                                                                                                                                                                                     0.0,
              0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
            3.212, 4.156,                                                                                                                                                                         5.201, 6.127, 6.999, 7.863, 8.699, 9.519,
           10.601,11.677,                                                                                                  12.549,13.377,14.189,14.875,15.776,16.578,17.370,18.157,18.835,19.719,20.623,21.549,22.489,23.439,24.397,25.359,
           26.331,27.308,                                                                                                  28.285,29.264,30.240,31.220,32.203,33.181,34.163,35.143,36.127,37.113,38.100,39.088,40.079,41.071,42.064,43.060,
           44.058,45.058,46.059,47.065,48.075,49.081,50.089,51.097,52.106,53.112,54.126,55.137,56.150,57.163,58.178,59.193,60.207,61.222,62.238,63.256,64.275,65.296,66.318,67.341,68.367,69.395,70.425,71.456,72.489,73.525,74.562,75.601,
           76.642,77.685,78.730,79.777,80.827,81.876,82.930,83.985,85.041,86.100,87.161,88.222,89.286,90.352,91.418,92.487,93.558,94.629,95.702,96.770,97.844,98.920,99.997,101.074,102.151,103.230,104.308,105.386,106.464,107.541,108.617,109.691
  };
  if (Z<=118) return Zeff[Z-1];
  else exit(-1);
}

double get_Zeff_3p_1_2(const int& Z){
  //Get the shielding constant according to the data in https://research.unl.pt/ws/portalfiles/portal/13073849/ADNDT_2016_7_Revision_2.pdf
  const double Zeff[] = {    0.0,                                                                                                                                                                                                                     0.0,
              0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
              0.0,   0.0,                                                                                                                                                                         3.642, 4.543, 5.386, 6.041, 6.786, 7.541,
            8.726, 9.837,                                                                                                  10.673,11.471,12.251,12.864,13.766,14.537,15.297,16.053,16.667,17.544,18.468,19.420,20.370,21.342,22.324,23.298,
           24.283,25.271,                                                                                                  26.261,27.250,28.238,29.226,30.214,31.201,32.189,33.176,34.164,35.152,36.144,37.138,38.131,39.129,40.128,41.128,
           42.130,43.134,44.140,45.148,46.159,47.170,48.182,49.195,50.214,51.224,52.238,53.253,54.272,55.292,56.313,57.337,58.356,59.377,60.398,61.421,62.446,63.473,64.501,65.531,66.563,67.598,68.635,69.674,70.715,71.759,72.805,73.853,
           74.904,75.957,77.012,78.069,79.130,80.193,81.258,82.327,83.398,84.470,85.546,86.624,87.706,88.789,89.876,90.965,92.057,93.151,94.248,95.346,96.448,97.554,98.662,99.772,100.886,102.003,103.122,104.245,105.370,106.497,107.627,108.759
  };
  if (Z<=118) return Zeff[Z-1];
  else exit(-1);
}

double get_Zeff_3p_3_2(const int& Z){
  //Get the shielding constant according to the data in https://research.unl.pt/ws/portalfiles/portal/13073849/ADNDT_2016_7_Revision_2.pdf
  const double Zeff[] = {    0.0,                                                                                                                                                                                                                     0.0,
              0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
              0.0,   0.0,                                                                                                                                                                           0.0, 4.541, 5.376, 6.083, 6.781, 7.502,
            8.685, 9.791,                                                                                                  10.633,11.415,12.178,12.774,13.666,14.422,15.163,15.898,16.495,17.356,18.257,19.184,20.136,21.083,22.037,23.012,
           23.977,24.946,                                                                                                  25.914,26.882,27.848,28.812,29.776,30.738,31.698,32.656,33.614,34.573,35.531,36.491,37.454,38.416,39.379,40.344,
           41.309,42.275,43.241,44.211,45.186,46.153,47.120,48.085,49.047,50.008,50.978,51.942,52.904,53.865,54.825,55.785,56.742,57.699,58.656,59.613,60.569,61.527,62.484,63.440,64.398,65.356,66.314,67.273,68.233,69.193,70.154,71.115,
           72.076,73.038,74.001,74.963,75.927,76.891,77.855,78.819,79.783,80.747,81.712,82.676,83.641,84.605,85.570,86.535,87.499,88.463,89.527,90.391,91.355,92.319,93.283,94.246,95.209,96.172,97.135,98.098,99.060,100.023,100.984,101.946
  };
  if (Z<=118) return Zeff[Z-1];
  else exit(-1);
}

double get_Zeff_3d_3_2(const int& Z){
  //Get the shielding constant according to the data in https://research.unl.pt/ws/portalfiles/portal/13073849/ADNDT_2016_7_Revision_2.pdf
  const double Zeff[] = {    0.0,                                                                                                                                                                                                                     0.0,
              0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
              0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
              0.0,   0.0,                                                                                                   6.216, 7.146, 7.904, 7.631, 9.246, 9.670,10.240,10.829,10.584,11.978,13.272,14.501,15.694,16.822,17.929,19.018,
           20.096,21.166,                                                                                                  22.234,23.290,24.335,25.361,26.384,27.406,28.416,29.429,30.430,31.428,32.424,33.419,34.413,35.405,36.396,37.387,
           38.377,39.368,40.358,41.345,42.322,43.308,44.294,45.280,46.275,47.253,48.220,49.196,50.176,51.156,52.138,53.121,54.092,55.062,56.032,57.001,57.969,58.939,59.908,60.876,61.844,62.813,63.782,64.752,65.721,66.691,67.661,68.632,
           69.602,70.573,71.545,72.516,73.487,74.459,75.431,76.403,77.379,78.351,79.322,80.295,81.268,82.242,83.215,84.189,85.161,86.133,87.105,88.077,89.048,90.020,90.991,91.961,92.932,93.902,94.872,95.842,96.811,97.781,98.749,99.718
  };
  if (Z<=118) return Zeff[Z-1];
  else exit(-1);
}

double get_Zeff_3d_5_2(const int& Z){
  //Get the shielding constant according to the data in https://research.unl.pt/ws/portalfiles/portal/13073849/ADNDT_2016_7_Revision_2.pdf
  const double Zeff[] = {    0.0,                                                                                                                                                                                                                     0.0,
              0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
              0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
              0.0,   0.0,                                                                                                     0.0, 7.139, 7.845, 7.608, 9.225, 9.779,10.300,10.805,10.479,11.886,13.228,14.457,15.613,16.751,17.850,18.924,
           19.998,21.062,                                                                                                  22.117,23.160,24.199,25.232,26.247,27.251,28.253,29.271,30.264,31.254,32.242,33.228,34.212,35.195,36.176,37.157,
           38.137,39.116,40.094,41.078,42.068,43.044,44.016,44.985,45.944,46.907,47.886,48.851,49.811,50.768,51.723,52.674,53.626,54.578,55.528,56.479,57.430,58.379,59.329,60.279,61.228,62.177,63.125,64.074,65.022,65.970,66.918,67.866,
           68.814,69.762,70.709,71.657,72.606,73.554,74.502,75.450,76.395,77.341,78.288,79.233,80.178,81.122,82.066,83.009,83.952,84.893,85.835,86.776,87.716,88.656,89.595,90.534,91.472,92.409,93.346,94.282,95.217,96.151,97.085,98.018
  };
  if (Z<=118) return Zeff[Z-1];
  else exit(-1);
}

double get_ionization_energy_1s(const int& Z){
  //#table from X-ray data booklet, in eV
  const double ionization_energies[] = {  13.6,                                                                                                                                                                                                                    24.6,
                           54.7, 111.5,                                                                                                                                                                         188.0, 284.2, 409.9, 543.1, 696.7, 870.2,
                         1070.8,1303.0,                                                                                                                                                                        1559.6,1839.0,2145.5,2472.0,2822.4,3205.9,
                         3608.4,4038.5,                                                                                                    4492,  4966,  5465,  5989,  6539,  7112,  7709,  8333,  8979,  9659, 10367, 11103, 11867, 12658, 13474, 14326,
                          15200, 16105,                                                                                                   17038, 17998, 18986, 20000, 21044, 22117, 23220, 24350, 25514, 26711, 27940, 29200, 30491, 31814, 33169, 34561,
                          35985, 37441, 38925, 40443, 41991, 43569, 45184, 46834, 48519, 50239, 51996, 53789, 55618, 57486, 59390, 61332, 63314, 65351, 67416, 69525, 71676, 73871, 76111, 78395, 80725, 83102, 85530, 88005, 90524, 93105, 95730, 98404,
                         101137,103922,106755,109651,112601,115606
  };
  if (Z<=92)
    return ionization_energies[Z-1];
  else{
    //#extrapoalte from fourth order regression of tabulated data
    const double a1 =  4.00729846e-04;
    const double a2 = -2.27683735e-02;
    const double a3 =  1.30249774e+01;
    const double a4 = -6.43729353e+01;
    const double c  =  1.99592659e+02;
    return a1 * pow(Z,4) + a2 * pow(Z,3) + a3 * pow(Z,2) + a4 * Z + c;
  }
  /*import np as np
  from scipy.optimize import curve_fit
  def func(x,a1,a2,a3,a4,c):
    return a1*x*x*x*x + a2*x*x*x + a3*x*x + a4*x + c
  x = range(2,93)
  result,resultcov = curve_fit(func,x,ionization_energies[1:])
  print (result)
  print(resultcov)
  coef = np.polyfit(x,ionization_energies[1:],3)
  funct = np.polyval(coef,Z)
  print(coef)
  import matplotlib.pyplot as plt
  plt.scatter(x,ionization_energies[1:],marker="*")
  plt.plot(x,func(x, *result), 'r-')
  plt.show()
  */
}

double get_ionization_energy_2s(const int& Z){
  //#table from X-ray data booklet, in eV
  const double ionization_energies[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,  37.3,  41.6, 0.000,  48.5,
                           63.5,  88.7,                                                                                                                                                                         117.8, 149.7, 189.0, 230.9, 270.0, 326.3,
                          378.6, 438.4,                                                                                                   498.0, 560.9, 626.7, 696.0, 769.1, 844.6, 925.1,1008.6,1096.7,1196.2,1299.0,1414.6,1527.0,1652.0,  1782,  1921,
                           2065,  2216,                                                                                                    2373,  2532,  2698,  2866,  3043,  3224,  3412,  3604,  3806,  4018,  4238,  4465,  4698,  4939,  5188,  5453,
                           5714,  5989,  6266,  6549,  6835,  7126,  7428,  7737,  8052,  8376,  8708,  9046,  9394,  9751, 10116, 10486, 10870, 11271, 11682, 12100, 12527, 12968, 13419, 13880, 14353, 14839, 15347, 15861, 16388, 16939, 17493, 18049,
                          18639, 19237, 19840, 20472, 21105, 21757
  };
  //# Coefficients from fit to data
  const double a1 = 1.44070614e+01;
  const double a2 = 1.49753416e-01;
  const double a3 = 1.15499698e-02;
  const double a4 = 1.61081136e+00;
  const double a5 =-2.27553738e+01;
  const double b1 = 4.99999128e+01;
  const double c  = 1.27894569e+02;
  if (Z<=92){
    double result = ionization_energies[Z-1];
    if (result == 0.0)
      if (Z>=3)
        result = a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*pow(Z,2) + a5*Z + c;
    return result;
  }
  else
    //#extrapoalte from fourth order regression of tabulated data
    return a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*pow(Z,2) + a5*Z + c;
    /*import np as np
    from scipy.optimize import curve_fit
    def func(x,a1,a2,a3,a4,a5,b1,c):
      return a1*np.exp2(a2*(x-b1)) + a3*x*x*x + a4*x*x + a5*x + c
    x = range(10,93)
    x2 = range(1,103)
    y = ionization_energies[9:]
    result,resultcov = curve_fit(func,x,y,p0=[8.36040351e+02,5.75329301e-02,2.15148012e+00,0,0,2.51302235e+01,0],bounds=([0,0,-np.inf,-np.inf,-np.inf,-50,-500],[2000,1000,1000,2000,2000,50,500]))
    print (result)
    print(resultcov)
    #coef = np.polyfit(x,ionization_energies[9:],3)
    #funct = np.polyval(coef,Z)
    #print(coef)
    import matplotlib.pyplot as plt
    plt.scatter(x,y,marker="*")
    plt.plot(x2,func(x2, *result), 'r-')
    plt.show() */
}

double get_ionization_energy_2p1_2(const int& Z){
  //#table from X-ray data booklet, in eV
  const double ionization_energies[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0, 0.000,  21.7,
                          30.65, 49.78,                                                                                                                                                                         72.95, 99.82,   136, 163.6,   202, 250.6,
                          297.3, 349.7,                                                                                                   403.6, 460.2, 519.8, 583.8, 649.9, 719.9, 793.2, 870.0, 952.3,1044.9,1143.2,1248.1,1359.1,1474.3,  1596,1730.9,
                           1864,  2007,                                                                                                    2156,  2307,  2456,  2625,  2793,  2967,  3146,  3330,  3524,  3727,  3938,  4156,  4380,  4612,  4852,  5107,
                           5359,  5624,  5891,  6164,  6440,  6722,  7013,  7312,  7617,  7930,  8252,  8581,  8918,  9264,  9617,  9978, 10349, 10739, 11136, 11544, 11959, 12385, 12824, 13273, 13734, 14209, 14698, 15200, 15711, 16244, 16785, 17337,
                          17907, 18484, 19083, 19693, 20314, 20948
  };
  //# Coefficients from fit to data
  const double a1 = 1.30802855e+01;
  const double a2 = 1.51798476e-01;
  const double a3 = 1.14195659e-02;
  const double a4 = 1.58049679e+00;
  const double a5 =-2.74880264e+01;
  const double b1 = 4.99597550e+01;
  const double c  = 1.48993379e+02;
  if (Z<=92){
    double result = ionization_energies[Z-1];
    if (result == 0.0)
      if (Z>=5)
        result = a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
    return result;
  }
  else
    //#extrapoalte from fourth order regression of tabulated data
    return a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
    /*
    import np as np
    from scipy.optimize import curve_fit
    def func(x,n1,n2,n3,n4,n5,m1,o):
      return n1*np.exp2(n2*(x-m1)) + n3*x*x*x + n4*x*x + n5*x + o
    x = range(10,93)
    x2 = range(1,103)
    y = ionization_energies[9:]
    result,resultcov = curve_fit(func,x,y,p0=[a1,a2,a3,a4,a5,b1,c],bounds=([0,0,-np.inf,-np.inf,-np.inf,-50,-500],[2000,1000,1000,2000,2000,50,500]))
    print (result)
    print(resultcov)
    #coef = np.polyfit(x,ionization_energies[9:],3)
    #funct = np.polyval(coef,Z)
    #print(coef)
    import matplotlib.pyplot as plt
    plt.scatter(x,y,marker="*")
    plt.plot(x2,func(x2, *result), 'r-')
    plt.show()*/
}

double get_ionization_energy_2p3_2(const int& Z){
  //#table from X-ray data booklet, in eV
  const double ionization_energies[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0, 0.000,  21.6,
                          30.81, 49.50,                                                                                                                                                                         72.55, 99.42,   135, 163.6,   200, 248.4,
                          294.6, 346.2,                                                                                                   398.7, 453.8, 512.1, 574.1, 638.7, 706.8, 778.1, 852.7, 932.7,1021.8,1116.4,1217.0,1323.6,1433.9,  1550,1678.4,
                           1804,  1940,                                                                                                    2080,  2223,  2371,  2520,  2677,  2838,  3004,  3173,  3351,  3538,  3730,  3929,  4132,  4341,  4557,  4786,
                           5012,  5247,  5483,  5723,  5964,  6208,  6459,  6716,  6977,  7243,  7514,  7790,  8071,  8358,  8648,  8944,  9244,  9561,  9881, 10207, 10535, 10871, 11215, 11564, 11919, 12284, 12658, 13035, 13419, 13814, 14214, 14619,
                          15031, 15444, 15871, 16300, 16733, 17166
  };
  //# Coefficients from fit to data
  const double a1 = 6.26079141e-01;
  const double a2 = 2.48380219e-01;
  const double a3 = 3.01836854e-03;
  const double a4 = 2.15937528e+00;
  const double a5 =-4.18326000e+01;
  const double b1 = 6.02739274e+01;
  const double c  = 2.60817688e+02;
  if (Z<=92){
    double result = ionization_energies[Z-1];
    if (result == 0.0)
      if (Z>=7)
        result = a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
    return result;
  }
  else
    //#extrapoalte from fourth order regression of tabulated data
    return a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
    /*from scipy.optimize import curve_fit
    def func(x,n1,n2,n3,n4,n5,m1,o):
      return n1*np.exp2(n2*(x-m1)) + n3*x*x*x + n4*x*x + n5*x + o
    x = range(10,93)
    x2 = range(7,103)
    y = ionization_energies[9:]
    result,resultcov = curve_fit(func,x,y,p0=[a1,a2,a3,a4,a5,b1,c],bounds=([0,0,-100,-100,-1000,-100,0],[2000,100,100,1000,1000,100,500]),maxfev=100000000)
    print (result)
    print(resultcov)
    #coef = np.polyfit(x,ionization_energies[9:],3)
    #funct = np.polyval(coef,Z)
    #print(coef)
    import matplotlib.pyplot as plt
    plt.scatter(x,y,marker="*")
    plt.plot(x2,func(x2, *result), 'r-')
    plt.show()*/
}

double get_ionization_energy_3s(const int& Z){
  //#table from X-ray data booklet, in eV
  const double ionization_energies[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0, 0.000,   0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,  29.3,
                           34.8,  44.3,                                                                                                    51.1,  58.7,  66.3,  74.1,  82.3,  91.3, 101.0, 110.8, 122.5, 139.8, 159.5, 180.1, 204.7, 229.6,   257, 292.8,
                          326.7, 358.7,                                                                                                   392.0, 430.3, 466.6, 506.3,   544, 586.1, 628.1, 671.6, 719.0, 772.0, 827.2, 884.7,   946,  1006,  1072,1148.7,
                           1211,  1293,  1362,  1436,  1511,  1575,   0.0,  1723,  1800,  1881,  1968,  2047,  2128,  2207,  2307,  2398,  2491,  2601,  2708,  2820,  2932,  3049,  3174,  3296,  3425,  3562,  3704,  3851,  3999,  4149,  4317,  4482,
                           4652,  4822,  5002,  5182,  5367,  5548
  };
  //# Coefficients from fit to data
  const double a1 = 1.43433111e+02;
  const double a2 = 5.93352595e-02;
  const double a3 = 3.63930366e-03;
  const double a4 = 1.96661938e-01;
  const double a5 =-8.55996542e+00;
  const double b1 = 2.98052563e+01;
  const double c  = 6.63407582e-05;
  if (Z<=92){
    double result = ionization_energies[Z-1];
    if (result == 0.0)
      if (Z>=18)
        result = a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
    return result;
  }
  else
    //#extrapoalte from fourth order regression of tabulated data
    return a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
  /*from scipy.optimize import curve_fit
  def func(x,n1,n2,n3,n4,n5,m1,o):
    return n1*np.exp2(n2*(x-m1)) + n3*x*x*x + n4*x*x + n5*x + o
  x = range(19,93)
  x2 = range(7,103)
  y = ionization_energies[18:]
  result,resultcov = curve_fit(func,x,y,p0=[a1,a2,a3,a4,a5,b1,c],bounds=([0,0,-100,-100,-1000,-100,0],[2000,100,100,1000,1000,100,500]),maxfev=100000000)
  print (result)
  print(resultcov)
  #coef = np.polyfit(x,ionization_energies[9:],3)
  #funct = np.polyval(coef,Z)
  #print(coef)
  import matplotlib.pyplot as plt
  plt.scatter(x,y,marker="*")
  plt.plot(x2,func(x2, *result), 'r-')
  plt.show()*/
}

double get_ionization_energy_3p_1_2(const int& Z){
  //#table from X-ray data booklet, in eV
  const double ionization_energies[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0, 0.000,   0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,  15.9,
                           18.3,  25.4,                                                                                                    28.3,  32.6,  37.2,  42.2,  47.2,  52.7,  58.9,  68.0,  77.3,  91.4, 103.5, 124.9, 146.9, 166.5,   189, 222.2,
                          248.7, 280.3,                                                                                                   310.6, 343.5, 376.1, 411.6, 447.6, 483.5, 521.3, 559.9, 603.8, 652.6, 703.2, 756.5, 812.7, 870.8,   931,1002.1,
                           1071,  1137,  1209,  1274,  1337,  1403,  1471,  1541,  1614,  1688,  1768,  1842,  1923,  2006,  2090,  2173,  2264,  2365,  2469,  2575,  2682,  2792,  2909,  3027,  3148,  3279,  3416,  3554,  3696,  3854,  4008,  4159,
                           4327,  4490,  4656,  4830,  5001,  5182
  };
  //# Coefficients from fit to data
  const double a1 = 1.70056228e+02;
  const double a2 = 5.84139911e-02;
  const double a3 =-1.06353074e-02;
  const double a4 = 1.50935984e+00;
  const double a5 =-6.65267026e+01;
  const double b1 = 2.65575820e+00;
  const double c  = 4.99863585e+02;
  if (Z<=92){
    double result = ionization_energies[Z-1];
    if (result == 0.0)
      if (Z>=18)
        result = a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
    return result;
  }
  else
    //#extrapoalte from fourth order regression of tabulated data
    return a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
  /*from scipy.optimize import curve_fit
  def func(x,n1,n2,n3,n4,n5,m1,o):
    return n1*np.exp2(n2*(x-m1)) + n3*x*x*x + n4*x*x + n5*x + o
  x = range(19,93)
  x2 = range(7,103)
  y = ionization_energies[18:]
  result,resultcov = curve_fit(func,x,y,p0=[a1,a2,a3,a4,a5,b1,c],bounds=([0,0,-100,-100,-1000,-100,0],[2000,100,100,1000,1000,100,500]),maxfev=100000000)
  print (result)
  print(resultcov)
  #coef = np.polyfit(x,ionization_energies[9:],3)
  #funct = np.polyval(coef,Z)
  #print(coef)
  import matplotlib.pyplot as plt
  plt.scatter(x,y,marker="*")
  plt.plot(x2,func(x2, *result), 'r-')
  plt.show()*/
}

double get_ionization_energy_3p_3_2(const int& Z){
  //#table from X-ray data booklet, in eV
  const double ionization_energies[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0, 0.000,   0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,  15.7,
                           18.3,  25.4,                                                                                                    28.3,  32.6,  37.2,  42.2,  47.2,  52.7,  59.9,  66.2,  75.1,  88.6, 100.0, 120.8, 141.2, 160.7,   182, 214.4,
                          239.1, 270.0,                                                                                                   298.8, 329.8, 360.6, 394.0, 417.7, 461.4, 496.5, 532.3, 573.0, 618.4, 665.3, 714.6, 766.4, 820.0,   875, 940.6,
                           1003,  1063,  1128,  1187,  1242,  1297,  1357,  1420,  1481,  1544,  1611,  1676,  1741,  1812,  1885,  1950,  2024,  2108,  2194,  2281,  2367,  2457,  2551,  2645,  2743,  2847,  2957,  3066,  3177,  3302,  3426,  3538,
                           3663,  3792,  3909,  4046,  4174,  4303
  };
  //# Coefficients from fit to data
  const double a1 = 2.09817052e+02;
  const double a2 = 5.25431495e-02;
  const double a3 =-1.34329853e-02;
  const double a4 = 1.70347853e+00;
  const double a5 =-7.71441584e+01;
  const double b1 =-4.18496392e+00;
  const double c  = 4.99999983e+02;
  if (Z<=92){
    double result = ionization_energies[Z-1];
    if (result == 0.0)
      if (Z>=18)
        result = a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
    return result;
  }
  else
    //#extrapoalte from fourth order regression of tabulated data
    return a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
  /*from scipy.optimize import curve_fit
  def func(x,n1,n2,n3,n4,n5,m1,o):
    return n1*np.exp2(n2*(x-m1)) + n3*x*x*x + n4*x*x + n5*x + o
  x = range(19,93)
  x2 = range(7,103)
  y = ionization_energies[18:]
  result,resultcov = curve_fit(func,x,y,p0=[a1,a2,a3,a4,a5,b1,c],bounds=([0,0,-100,-100,-1000,-100,0],[2000,100,100,1000,1000,100,500]),maxfev=100000000)
  print (result)
  print(resultcov)
  #coef = np.polyfit(x,ionization_energies[9:],3)
  #funct = np.polyval(coef,Z)
  #print(coef)
  import matplotlib.pyplot as plt
  plt.scatter(x,y,marker="*")
  plt.plot(x2,func(x2, *result), 'r-')
  plt.show()*/
}

double get_ionization_energy_3d_3_2(const int& Z){
  //#table from X-ray data booklet, in eV
  const double ionization_energies[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0, 0.000,   0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                            0.0,   0.0,                                                                                                     0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  10.2,  18.7,  29.8,  41.7,  55.5,    70,  95.0,
                          113.0, 136.0,                                                                                                   157.7, 181.1, 205.0, 231.1, 257.6, 284.2, 311.9, 340.5, 374.0, 411.9, 451.4, 493.2, 537.5, 583.4, 630.8, 689.0,
                          740.5, 795.7,   853, 902.4, 948.3,1003.3,  1052,1110.9,1158.6,1221.9,1276.9,  1333,  1392,  1453,  1515,  1576,  1639,  1716,  1793,  1872,  1949,  2031,  2116,  2202,  2291,  2385,  2485,  2586,  2688,  2798,  2909,  3022,
                           3136,  3248,  3370,  3491,  3611,  3728
  };
  //# Coefficients from fit to data
  const double a1 = 2.19640503e+02;
  const double a2 = 5.02997104e-02;
  const double a3 =-1.30856531e-02;
  const double a4 = 1.59217604e+00;
  const double a5 =-8.09352821e+01;
  const double b1 =-8.92936626e+00;
  const double c  = 4.99999918e+02;
  if (Z<=92){
    double result = ionization_energies[Z-1];
    if (result == 0.0)
      if (Z>=18)
        result = a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
    return result;
  }
  else
    //#extrapoalte from fourth order regression of tabulated data
    return a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
  /*from scipy.optimize import curve_fit
  def func(x,n1,n2,n3,n4,n5,m1,o):
    return n1*np.exp2(n2*(x-m1)) + n3*x*x*x + n4*x*x + n5*x + o
  x = range(19,93)
  x2 = range(7,103)
  y = ionization_energies[18:]
  result,resultcov = curve_fit(func,x,y,p0=[a1,a2,a3,a4,a5,b1,c],bounds=([0,0,-100,-100,-1000,-100,0],[2000,100,100,1000,1000,100,500]),maxfev=100000000)
  print (result)
  print(resultcov)
  #coef = np.polyfit(x,ionization_energies[9:],3)
  #funct = np.polyval(coef,Z)
  #print(coef)
  import matplotlib.pyplot as plt
  plt.scatter(x,y,marker="*")
  plt.plot(x2,func(x2, *result), 'r-')
  plt.show()*/
}
  
double get_ionization_energy_3d_5_2(const int& Z){
  //#table from X-ray data booklet, in eV
  const double ionization_energies[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0, 0.000,   0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,   0.0,   0.0,   0.0,   0.0,   0.0,
                            0.0,   0.0,                                                                                                     0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,  10.1,  18.7,  29.2,  41.7,  54.6,    69,  93.8,
                            112, 134.2,                                                                                                   155.8, 178.8, 202.3, 227.9, 253.9, 280.0, 307.2, 335.2, 368.3, 405.2, 443.9, 484.9, 528.2, 573.0, 619.3, 676.4,
                          726.6, 780.5,   836, 883.8, 928.8, 980.4,  1027,1083.4,1127.5,1189.6,1241.1,1292.6,  1351,  1409,  1468,  1528,  1589,  1662,  1735,  1809,  1883,  1960,  2040,  2122,  2206,  2295,  2389,  2484,  2580,  2683,  2787,  2892,
                           3000,  3105,  3219,  3332,  3442,  3552
  };
  //# Coefficients from fit to data
  const double a1 = 2.19703182e+02;
  const double a2 = 4.99554795e-02;
  const double a3 =-1.31912578e-02;
  const double a4 = 1.60333067e+00;
  const double a5 =-8.08169197e+01;
  const double b1 =-8.82796315e+00;
  const double c  = 4.99877841e+02;
  if (Z<=92){
    double result = ionization_energies[Z-1];
    if (result == 0.0)
      if (Z>=18)
        result = a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
    return result;
  }
  else
    //#extrapoalte from fourth order regression of tabulated data
    return a1*std::exp2(a2*(Z-b1)) + a3*pow(Z,3) + a4*Z*Z + a5*Z + c;
  /*from scipy.optimize import curve_fit
  def func(x,n1,n2,n3,n4,n5,m1,o):
    return n1*np.exp2(n2*(x-m1)) + n3*x*x*x + n4*x*x + n5*x + o
  x = range(19,93)
  x2 = range(7,103)
  y = ionization_energies[18:]
  result,resultcov = curve_fit(func,x,y,p0=[a1,a2,a3,a4,a5,b1,c],bounds=([0,0,-100,-100,-1000,-100,0],[2000,100,100,1000,1000,100,500]),maxfev=100000000)
  print (result)
  print(resultcov)
  #coef = np.polyfit(x,ionization_energies[9:],3)
  #funct = np.polyval(coef,Z)
  #print(coef)
  import matplotlib.pyplot as plt
  plt.scatter(x,y,marker="*")
  plt.plot(x2,func(x2, *result), 'r-')
  plt.show()*/
}

double get_line_width_K(const int& Z){
  //#returns line width in eV
  //#according to Kraus & Oliver J. Phys. Chem. Ref. Data, Vol.8, No.2, 1979
  const double Linewidths[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0, 0.000,  0.24,
                           0.30,  0.36,                                                                                                                                                                          0.42,  0.48,  0.53,  0.59,  0.64,  0.68,
                           0.74,  0.81,                                                                                                    0.86,  0.94,  1.01,  1.08,  1.16,  1.25,  1.33,  1.44,  1.55,  1.67,  1.82,  1.96,  2.14,  2.33,  2.52,  2.75,
                           2.99,  3.25,                                                                                                    3.52,  3.84,  4.14,  4.52,  4.91,  5.33,  5.77,  6.24,  6.75,  7.28,  7.91,  8.49,  9.16,  9.89,  10.6,  11.4,
                           12.3,  13.2,  14.1,  15.1,  16.2,  17.3,  18.5,  19.7,  21.0,  22.3,  23.8,  25.2,  26.8,  28.4,  30.1, 31.9,   33.7,  35.7,  37.7,  39.9,  42.1,  44.4,  46.8,  49.3,  52.0,  54.6,  57.4,  60.4,  63.4,  66.6,  69.8,  73.3,
                           76.8,  80.4,  84.1,  88.0,  91.9,  96.1, 100, 105, 109, 114, 119, 124, 129, 135, 140, 145, 150
};
  return Linewidths[Z-1];
}

double get_line_width_Lp_3_2(const int& Z){
  //#returns line width in eV
  //#according to Kraus & Oliver J. Phys. Chem. Ref. Data, Vol.8, No.2, 1979
  const double Linewidths[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0, 0.000,   0.1,
                           0.20,  0.41,                                                                                                                                                                          0.73,  1.03,  1.26,  1.49,  1.58,  1.63,
                           1.92,  2.07,                                                                                                    2.21,  2.34,  2.41,  2.54,  2.62,  2.76,  2.79,  2.89,  3.06,  3.28,  3.38,  3.53,  3.79,  3.94,  4.11,  4.28,
                           4.44,  4.67,                                                                                                    4.71,  4.78,  3.94,  4.25,  4.36,  4.58,  4.73,  4.93,  4.88,  4.87,  5.00,  2.97,  3.13,  3.32,  3.46,  3.64,
                           3.78,  3.92,  4.06,  4.21,  4.34,  4.52,  4.67,  4.80,  4.91,  5.05,  5.19,  5.25,  5.33,  5.43,  5.47, 5.53,   5.54,  5.63,  5.58,  5.61,  6.18,  7.25,  8.30,  9.39,  10.5,  11.3,  12.0,  12.2,  12.4,  12.6,  12.8,  13.1,
                           13.3,  13.4,  13.6,  13.7,  14.3,  14.0,    14, 13.5, 13.3, 13.6, 13.8, 14.0, 14.3, 14.4, 14.8, 15.1, 15.9
                        };
  return Linewidths[Z-1];
}

double get_line_width_Lp_1_2(const int& Z){
  //#returns line width in eV
  //#according to Kraus & Oliver J. Phys. Chem. Ref. Data, Vol.8, No.2, 1979
  const double Linewidths[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0,   0.0,   0.0,
                           0.00, 0.001,                                                                                                                                                                         0.004, 0.015, 0.032, 0.054, 0.083, 0.126,
                          0.152,  0.17,                                                                                                    0.19,  0.24,  0.26,  0.29,  0.34,  0.37,  0.43,  0.52,  0.62,  0.72,  0.83,  0.95,  1.03,  1.13,  1.21,  1.31,
                           1.43,  1.54,                                                                                                    1.65,  1.78,  1.87,  1.97,  2.08,  2.23,  2.35,  2.43,  2.57,  2.62,  2.72,  2.84,  3.00,  3.12,  3.25,  3.40,
                           3.51,  3.57,  3.68,  3.80,  3.89,  3.97,  4.06,  4.15,  4.23,  4.32,  4.43,  4.55,  4.66,  4.73,  4.79, 4.82,   4.92,  5.02,  5.15,  5.33,  5.48,  5.59,  5.69,  5.86,  6.00,  6.17,  6.32,  6.48,  6.67,  6.83,  7.01,  7.20,
                           7.47,  7.68,  7.95,  8.18,  8.75,  9.32,  9.91, 10.5, 10.9, 11.4, 11.8, 12.2, 12.7, 13.1, 13.6, 14.0, 14.4
  };
  return Linewidths[Z-1];
}
  
double get_line_width_Ls(const int& Z){
  //#returns line width in eV
  //#according to Kraus & Oliver J. Phys. Chem. Ref. Data, Vol.8, No.2, 1979
  const double Linewidths[] = {   0.0,                                                                                                                                                                                                                     0.0,
                            0.0,   0.0,                                                                                                                                                                           0.0,  0.00,   0.0,   0.0,   0.0,   0.0,
                           0.00, 0.001,                                                                                                                                                                         0.004, 0.014, 0.033, 0.054, 0.087, 0.128,
                          0.156,  0.17,                                                                                                    0.19,  0.22,  0.24,  0.27,  0.32,  0.36,  0.43,  0.48,  0.56,  0.65,  0.76,  0.82,  0.94,  1.00,  1.08,  1.17,
                           1.27,  1.39,                                                                                                    1.50,  1.57,  1.66,  1.78,  1.91,  2.00,  2.13,  2.25,  2.40,  2.50,  2.65,  2.75,  2.87,  2.95,  3.08,  3.13,
                           3.25,  3.32,  3.41,  3.48,  3.60,  3.65,  3.75,  3.86,  3.91,  4.01,  4.12,  4.17,  4.26,  4.35,  4.48,  4.60,  4.68,  4.80,  4.88,  4.98,  5.04,  5.16,  5.25,  5.31,  5.41,  5.50,  5.65,  5.81,  5.98,  6.13,  6.29,  6.41,
                           6.65,  6.82,  6.98,  7.13,  7.33,  7.43,  7.59,  7.82,  8.04,  8.26,  8.55,  8.75,  9.04,  9.33,  9.61,  9.90, 10.1
                         };
  return Linewidths[Z-1];
}

namespace kleemiss_peyerimhoff{
double one_minus_delta_edge(const int& Z, const int& n0, const int& l, const double& j){
  double Z_s_sq, e_ion;
  if (n0 == 1){
    if (l != 0) exit(-1);
    else{
      Z_s_sq = pow(get_Zeff_1s(Z),2);
      e_ion = get_ionization_energy_1s(Z);
      double delta_K = alpha_sq * Z_s_sq / 4 - e_ion/(Ryd_ener * Z_s_sq); // #21a in Hoenl
      return -delta_K;
    }
  }
  else if (n0 == 2){
    if (l == 0){
      Z_s_sq = pow(get_Zeff_2s(Z),2);
      e_ion = get_ionization_energy_2s(Z);
      double delta_l1 = alpha_sq * Z_s_sq * 0.3125 - 4*e_ion/(Ryd_ener * Z_s_sq); // #21a in Hoenl
      return -delta_l1;
    }
    else if (l == 1){
      if (j == 0.5){
        Z_s_sq = pow(get_Zeff_2p_1_2(Z),2);
        e_ion = get_ionization_energy_2p1_2(Z);
        double delta_l2 = alpha_sq * Z_s_sq * 0.3125 - 4*e_ion/(Ryd_ener * Z_s_sq); // #21a in Hoenl
        return -delta_l2;
      }
      else if (j == 1.5){
        Z_s_sq = pow(get_Zeff_2p_3_2(Z),2);
        e_ion = get_ionization_energy_2p3_2(Z);
        double delta_l2 = alpha_sq * Z_s_sq * 0.0625 - 4*e_ion/(Ryd_ener * Z_s_sq); // #21a in Hoenl
        return -delta_l2;
      }
      else exit(-1);
    }
    else exit(-1);
  }
  else if (n0 == 3){
    if (l == 0){
      Z_s_sq = pow(get_Zeff_3s(Z),2);
      e_ion = get_ionization_energy_3s(Z);
      double delta_m1 = alpha_sq * Z_s_sq * 0.25 - 9*e_ion/(Ryd_ener * Z_s_sq); // #21a in Hoenl
      return -delta_m1;
    }
    else if (l == 1){
      if (j == 0.5){
        Z_s_sq = pow(get_Zeff_3p_1_2(Z),2);
        e_ion = get_ionization_energy_3p_1_2(Z);
        double delta_m2 = alpha_sq * Z_s_sq * 0.25 - 9*e_ion/(Ryd_ener * Z_s_sq); // #21a in Hoenl
        return -delta_m2;
      }
      else if (j == 1.5){
        Z_s_sq = pow(get_Zeff_3p_3_2(Z),2);
        e_ion = get_ionization_energy_3p_3_2(Z);
        double delta_m3 = alpha_sq * Z_s_sq /12 - 9*e_ion/(Ryd_ener * Z_s_sq); // #21a in Hoenl
        return -delta_m3;
      }
      else exit(-1);
    }
    else if (l == 2){
      if (j == 1.5){
        Z_s_sq = pow(get_Zeff_3d_3_2(Z),2);
        e_ion = get_ionization_energy_3d_3_2(Z);
        double delta_m4 = alpha_sq * Z_s_sq /12 - 9*e_ion/(Ryd_ener * Z_s_sq); // #21a in Hoenl
        return -delta_m4;
      }
      else if (j == 2.5){
        Z_s_sq = pow(get_Zeff_3d_5_2(Z),2);
        e_ion = get_ionization_energy_3d_5_2(Z);
        double delta_m5 = alpha_sq * Z_s_sq /36 - 9*e_ion/(Ryd_ener * Z_s_sq); // #21a in Hoenl
        return -delta_m5;
      }
      else exit(-1);
    }
    else exit(-1);
  }
  else exit(-1);
  return -200;
}

double b(const int& n_0, const int& l_0 , const int& Z){
  double Z_eff = -20;
  if (n_0 == 1)
    Z_eff = get_Zeff_1s(Z);
  else if (n_0 == 2){
    if (l_0 == 0)
      Z_eff = get_Zeff_2s(Z);
    else if (l_0 == 1)
      Z_eff = get_Zeff_2p_1_2(Z);
    else if (l_0 == 2)
      Z_eff = get_Zeff_2p_3_2(Z);
  }
  else if (n_0 == 3){
    if (l_0 == 0)
      Z_eff = get_Zeff_3s(Z);
    else if (l_0 == 1)
      Z_eff = get_Zeff_3p_1_2(Z);
    else if (l_0 == 2)
      Z_eff = get_Zeff_3d_3_2(Z);
  }
  err_checkf(Z_eff != -20,"WRONG FUNCTION! " + std::to_string(Z) + " " + std::to_string(n_0) + " " +std::to_string(l_0) , std::cout);
  return Z_eff/(n_0*a0);
}

double N0_square(const double &b_){
  return pow(b_,3)/PI;
}

double N0(const double& b_){
  return std::sqrt(N0_square(b_));
}

double product_n_prime_from_z(const int& n_0, const double& z, const int& l){
  cdouble n_p = n_prime_from_z(z,n_0);
  double fact = 1.0;
  for (int nu=1; nu <= l; nu++)
    fact *= real(n_p * n_p) + nu * nu;
  double denom = 1.0 - std::exp(-2*PI*abs(n_p));
  return fact/denom;
}

//Introduces factors 2pi m_e /h^2 and m
double N_square_from_z(const int& l, const int& m, const double& b_, const int& n_0, const double& z){
  if (m > l)
    return 0;
  double result = (2*l+1)*ft_fun(l-m)/ft_fun(l+m) * constant_factor / PI * n_0 * b_ * product_n_prime_from_z(n_0,z,l);
  if (m >= 1)
    result *= 2;
  return result;
}

double N(const int& l, const int& m, const double& b_, const int& n_0, const double& z){
  return std::sqrt(N_square_from_z(l,m,b_,n_0,z));
}

double N_square(const int& l, const int& m, const double& b_, const int& n_0, const double& z){
  return N_square_from_z(l,m,b_,n_0,z);
}

double N_lm_from_z(const int& l, const int& m, const double& z, const double& b_, const int& n_0){
  return N(l,m,b_,n_0, z);
}

double N_lm_square_from_z(const int& l, const int& m, const double& z, const double& b_, const int& n_0){
  return N_square(l,m,b_,n_0,z);
}
}