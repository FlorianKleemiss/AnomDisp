#pragma once
#include "constants_and_atomic_properties.h"
#include "legendre_polynomials.h"
#include "matrix_coefficient.h"


double apply_angle_part_s_parallel(double a_l, double theta0, double alpha, int l){
  return a_l * alpha_coef(l,1,1,theta0,alpha);
}

double apply_angle_part_s_orthogonal(double a_l, double theta0, double alpha, int l){
  return a_l * beta_coef(l,1,1,theta0,alpha);
}

std::vector<double> apply_angle_part_p_parallel(std::vector<double> vals, int _k, double theta0, double alpha, int l){
  double ct0 =cos(theta0);
  double st0 =sin(theta0);
  double ca =cos(alpha);
  double sa =sin(alpha);
  if (_k == 0){
    vals[0] *= -st0 * alpha_coef(l,1,0,theta0,alpha);
    vals[1] *=  ct0 * alpha_coef(l,0,0,theta0,alpha) * ca;
    vals[2] *=  ct0 * alpha_coef(l,2,0,theta0,alpha) * ca;
    vals[3] *=  ct0 * alpha_bar_coef(l,2,0,theta0,alpha) * sa;
  }
  else if (_k == 1){
    vals[0] *= ct0 * alpha_coef(l,1,1,theta0,alpha);
    vals[1] *= st0 * alpha_coef(l,0,1,theta0,alpha) * ca;
    vals[2] *= st0 * alpha_coef(l,2,1,theta0,alpha) * ca;
    vals[3] *= st0 * alpha_bar_coef(l,2,1,theta0,alpha) * sa;
  }
  else if (_k == 2){
    vals[0] *= -st0 * alpha_coef(l,1,2,theta0,alpha);
    vals[1] *=  ct0 * alpha_coef(l,0,2,theta0,alpha) * ca - sa * beta_coef(l,0,2,theta0,alpha);
    vals[2] *=  ct0 * alpha_coef(l,2,2,theta0,alpha) * ca - sa * beta_coef(l,2,2,theta0,alpha);
    vals[3] *=  ct0 * alpha_bar_coef(l,2,2,theta0,alpha) * sa + ca * beta_bar_coef(l,2,2,theta0,alpha);
  }
  return vals;
}

std::vector<double> apply_angle_part_p_orthogonal(std::vector<double> vals,int _k,double theta0,double alpha,int l){
  double ct0 = cos(theta0);
  double st0 = sin(theta0);
  double ca = cos(alpha);
  double sa = sin(alpha);
  if (_k == 0){
    vals[0] *= 0;
    vals[1] *= -alpha_coef(l,0,0,theta0,alpha) * sa;
    vals[2] *= -alpha_coef(l,2,0,theta0,alpha) * sa;
    vals[3] *= alpha_bar_coef(l,2,0,theta0,alpha) * ca;
  }
  else if (_k == 1){
    vals[0] *= ct0 * beta_coef(l,1,1,theta0,alpha);
    vals[1] *= st0 * beta_coef(l,0,1,theta0,alpha) * ca;
    vals[2] *= st0 * beta_coef(l,2,1,theta0,alpha) * ca;
    vals[3] *= st0 * beta_bar_coef(l,2,1,theta0,alpha) * sa;
  }
  else if (_k == 2){
    vals[0] *= -st0 * beta_coef(l,1,2,theta0,alpha);
    vals[1] *=  ct0 * beta_coef(l,0,2,theta0,alpha) * ca + sa * alpha_coef(l,0,2,theta0,alpha);
    vals[2] *=  ct0 * beta_coef(l,2,2,theta0,alpha) * ca + sa * alpha_coef(l,2,2,theta0,alpha);
    vals[3] *=  ct0 * beta_bar_coef(l,2,2,theta0,alpha) * sa - ca * alpha_bar_coef(l,2,2,theta0,alpha);
  }
  return vals;
}

std::vector<double> apply_angle_part_d_parallel(std::vector<double> vals,int _k,double theta0, double alpha, int l){
  double ct0  = cos(theta0);
  double c2t0 = cos(2*theta0);
  double st0  = sin(theta0);
  double s2t0 = sin(2*theta0);
  double ca   = cos(alpha);
  double c2a  = cos(2*alpha);
  double sa   = sin(alpha);
  double s2a  = sin(2*alpha);
  if (_k == 0){
    vals[0] *=  c2t0 * alpha_coef(l,0,0,theta0,alpha) * ca;
    vals[1] *=  c2t0 * alpha_coef(l,2,0,theta0,alpha) * ca;
    vals[2] *=  c2t0 * sa * alpha_bar_coef(l,2,0,theta0,alpha);
    vals[3] *=  0.5*s2t0 * alpha_bar_coef(l,1,0,theta0,alpha) * s2a;
    vals[4] *=  0.5*s2t0 * alpha_bar_coef(l,3,0,theta0,alpha) * s2a;
    vals[5] *=  -std::sqrt(3.0)/2.0 * s2t0 * alpha_coef(l,1,0,theta0,alpha);
    vals[6] *=  0.5*s2t0*c2a *alpha_coef(l,1,0,theta0,alpha);
    vals[7] *=  0.5*s2t0*c2a *alpha_coef(l,3,0,theta0,alpha);
  }
  else if (_k == 1){
    vals[0] *=  st0*sa * beta_coef(l,0,1,theta0,alpha)      - 0.5*s2t0*ca*alpha_coef(l,0,1,theta0,alpha);
    vals[1] *=  st0*sa * beta_coef(l,2,1,theta0,alpha)      - 0.5*s2t0*ca*alpha_coef(l,2,1,theta0,alpha);
    vals[2] *= -st0*ca * beta_bar_coef(l,2,1,theta0,alpha)  - 0.5*s2t0*sa*alpha_bar_coef(l,2,1,theta0,alpha);
    vals[3] *=  ct0*c2a * beta_bar_coef(l,1,1,theta0,alpha) + (1-0.5*st0*st0)*s2a*alpha_bar_coef(l,1,1,theta0,alpha);
    vals[4] *=  ct0*c2a * beta_bar_coef(l,3,1,theta0,alpha) + (1-0.5*st0*st0)*s2a*alpha_bar_coef(l,3,1,theta0,alpha);
    vals[5] *=                                                std::sqrt(3.0)*0.5*st0*st0 * alpha_coef(l,1,1,theta0,alpha);
    vals[6] *= -ct0*s2a * beta_coef(l,1,1,theta0,alpha)     + (1-0.5*st0*st0)*c2a*alpha_coef(l,1,1,theta0,alpha);
    vals[7] *= -ct0*s2a * beta_coef(l,3,1,theta0,alpha)     + (1-0.5*st0*st0)*c2a*alpha_coef(l,3,1,theta0,alpha);
  }
  else if (_k == 2){
    vals[0] *=  c2t0*ca * alpha_coef(l,0,2,theta0,alpha)          - ct0*sa*beta_coef(l,0,2,theta0,alpha);
    vals[1] *=  c2t0*ca * alpha_coef(l,2,2,theta0,alpha)          - ct0*sa*beta_coef(l,2,2,theta0,alpha);
    vals[2] *=  c2t0*sa * alpha_bar_coef(l,2,2,theta0,alpha)      + ct0*ca*beta_bar_coef(l,2,2,theta0,alpha);
    vals[3] *=  0.5*s2t0*s2a * alpha_bar_coef(l,1,2,theta0,alpha) + st0*c2a*beta_bar_coef(l,1,2,theta0,alpha);
    vals[4] *=  0.5*s2t0*s2a * alpha_bar_coef(l,3,2,theta0,alpha) + st0*c2a*beta_bar_coef(l,3,2,theta0,alpha);
    vals[5] *=                                                    - std::sqrt(3.0)*0.5*s2t0 * alpha_coef(l,1,2,theta0,alpha);
    vals[6] *=  0.5*s2t0*c2a * alpha_coef(l,1,2,theta0,alpha)     - st0*s2a*beta_coef(l,1,2,theta0,alpha);
    vals[7] *=  0.5*s2t0*c2a * alpha_coef(l,3,2,theta0,alpha)     - st0*s2a*beta_coef(l,3,2,theta0,alpha);
  }
  else if (_k == 3){
    vals[0] *=  st0*sa * beta_coef(l,0,3,theta0,alpha)      - 0.5*s2t0*ca*alpha_coef(l,0,3,theta0,alpha);
    vals[1] *=  st0*sa * beta_coef(l,2,3,theta0,alpha)      - 0.5*s2t0*ca*alpha_coef(l,2,3,theta0,alpha);
    vals[2] *= -st0*ca * beta_bar_coef(l,2,3,theta0,alpha)  - 0.5*s2t0*sa*alpha_bar_coef(l,2,3,theta0,alpha);
    vals[3] *=  ct0*c2a * beta_bar_coef(l,1,3,theta0,alpha) + (1-0.5*st0*st0)*s2a*alpha_bar_coef(l,1,3,theta0,alpha);
    vals[4] *=  ct0*c2a * beta_bar_coef(l,3,3,theta0,alpha) + (1-0.5*st0*st0)*s2a*alpha_bar_coef(l,3,3,theta0,alpha);
    vals[5] *=                                                std::sqrt(3.0)*0.5*st0*st0 * alpha_coef(l,1,3,theta0,alpha);
    vals[6] *= -ct0*s2a * beta_coef(l,1,3,theta0,alpha)     + (1-0.5*st0*st0)*c2a*alpha_coef(l,1,3,theta0,alpha);
    vals[7] *= -ct0*s2a * beta_coef(l,3,3,theta0,alpha)     + (1-0.5*st0*st0)*c2a*alpha_coef(l,3,3,theta0,alpha);
  }
  else if (_k == 4){
    vals[0] *=  0.5*std::sqrt(3.0) *s2t0*ca * alpha_coef(l,0,1,theta0,alpha);
    vals[1] *=  0.5*std::sqrt(3.0) *s2t0*ca * alpha_coef(l,2,1,theta0,alpha);
    vals[2] *=  std::sqrt(3.0)*0.5*s2t0*sa * alpha_bar_coef(l,2,1,theta0,alpha);
    vals[3] *=  std::sqrt(3.0)*0.5*st0*st0*s2a * alpha_bar_coef(l,1,1,theta0,alpha);
    vals[4] *=  std::sqrt(3.0)*0.5*st0*st0*s2a * alpha_bar_coef(l,3,1,theta0,alpha);
    vals[5] *=  (1.5*ct0*ct0-0.5) * alpha_coef(l,1,1,theta0,alpha);
    vals[6] *=  std::sqrt(3.0)*0.5*st0*st0*c2a *alpha_coef(l,1,1,theta0,alpha);
    vals[7] *=  std::sqrt(3.0)*0.5*st0*st0*c2a *alpha_coef(l,3,1,theta0,alpha);
  }
  return vals;
}
std::vector<double> apply_angle_part_d_orthogonal(std::vector<double> vals,int _k, double theta0, double alpha, int l){
  double ct0  = cos(theta0);
  double c2t0 = cos(2*theta0);
  double st0  = sin(theta0);
  double s2t0 = sin(2*theta0);
  double ca   = cos(alpha);
  double c2a  = cos(2*alpha);
  double sa   = sin(alpha);
  double s2a  = sin(2*alpha);
  if (_k == 0){
    vals[0] *=  -ct0*sa * alpha_coef(l,0,0,theta0,alpha);
    vals[1] *=  -ct0*sa * alpha_coef(l,2,0,theta0,alpha);
    vals[2] *=  ct0 * ca * alpha_bar_coef(l,2,0,theta0,alpha);
    vals[3] *=  st0 * c2a * alpha_bar_coef(l,1,0,theta0,alpha);
    vals[4] *=  st0 * c2a * alpha_bar_coef(l,3,0,theta0,alpha);
    vals[5] *=  0;
    vals[6] *=  -st0*s2a * alpha_coef(l,1,0,theta0,alpha);
    vals[7] *=  -st0*s2a * alpha_coef(l,3,0,theta0,alpha);
  }
  else if (_k == 1){
    vals[0] *= 0.5*s2t0*ca * beta_coef(l,0,1,theta0,alpha)              + st0*sa*alpha_coef(l,0,1,theta0,alpha);
    vals[1] *= 0.5*s2t0*ca * beta_coef(l,2,1,theta0,alpha)              + st0*sa*alpha_coef(l,2,1,theta0,alpha);
    vals[2] *= 0.5*s2t0*sa * beta_bar_coef(l,2,1,theta0,alpha)          - st0*ca*alpha_bar_coef(l,2,1,theta0,alpha);
    vals[3] *= -(1-0.5*st0*st0)*s2a * beta_bar_coef(l,1,1,theta0,alpha) + ct0*c2a*alpha_bar_coef(l,1,1,theta0,alpha);
    vals[4] *= -(1-0.5*st0*st0)*s2a * beta_bar_coef(l,3,1,theta0,alpha) + ct0*c2a*alpha_bar_coef(l,3,1,theta0,alpha);
    vals[5] *= -std::sqrt(3.0)*0.5*st0*st0 * beta_coef(l,1,1,theta0,alpha);
    vals[6] *= -(1-0.5*st0*st0)*c2a * beta_coef(l,1,1,theta0,alpha)     - ct0*s2a*alpha_coef(l,1,1,theta0,alpha);
    vals[7] *= -(1-0.5*st0*st0)*c2a * beta_coef(l,3,1,theta0,alpha)     - ct0*s2a*alpha_coef(l,3,1,theta0,alpha);
  }
  else if (_k == 2){
    vals[0] *= c2t0*ca * beta_coef(l,0,2,theta0,alpha)           + ct0*sa*alpha_coef(l,0,2,theta0,alpha);
    vals[1] *= c2t0*ca * beta_coef(l,2,2,theta0,alpha)           + ct0*sa*alpha_coef(l,2,2,theta0,alpha); 
    vals[2] *= -ct0*ca * alpha_bar_coef(l,2,2,theta0,alpha)      + c2t0*sa*beta_bar_coef(l,2,2,theta0,alpha);
    vals[3] *= -st0*c2a * alpha_bar_coef(l,1,2,theta0,alpha) + 0.5*s2t0*s2a*beta_bar_coef(l,1,2,theta0,alpha);
    vals[4] *= -st0*c2a * alpha_bar_coef(l,3,2,theta0,alpha) + 0.5*s2t0*s2a*beta_bar_coef(l,3,2,theta0,alpha);
    vals[5] *= -std::sqrt(3.0)*0.5*s2t0 * beta_coef(l,1,2,theta0,alpha);
    vals[6] *= 0.5*s2t0*c2a * beta_coef(l,1,2,theta0,alpha)     + st0*s2a*alpha_coef(l,1,2,theta0,alpha);
    vals[7] *= 0.5*s2t0*c2a * beta_coef(l,3,2,theta0,alpha)     + st0*s2a*alpha_coef(l,3,2,theta0,alpha);
  }
  else if (_k == 3){
    vals[0] *= -0.5*s2t0*ca * beta_coef(l,0,3,theta0,alpha)      - st0*sa*alpha_coef(l,0,3,theta0,alpha);
    vals[1] *= -0.5*s2t0*ca * beta_coef(l,2,3,theta0,alpha)      - st0*sa*alpha_coef(l,2,3,theta0,alpha);
    vals[2] *= -0.5*s2t0*sa * beta_bar_coef(l,2,3,theta0,alpha)  + st0*ca*alpha_bar_coef(l,2,3,theta0,alpha);
    vals[3] *= (1-0.5*st0*st0)*s2a * beta_bar_coef(l,1,3,theta0,alpha) - ct0*c2a*alpha_bar_coef(l,1,3,theta0,alpha);
    vals[4] *= (1-0.5*st0*st0)*s2a * beta_bar_coef(l,3,3,theta0,alpha) - ct0*c2a*alpha_bar_coef(l,3,3,theta0,alpha);
    vals[5] *= std::sqrt(3.0)*0.5*st0*st0 * beta_coef(l,1,3,theta0,alpha);
    vals[6] *= (1-0.5*st0*st0)*c2a * beta_coef(l,1,3,theta0,alpha)     + ct0*s2a*alpha_coef(l,1,3,theta0,alpha);
    vals[7] *= (1-0.5*st0*st0)*c2a * beta_coef(l,3,3,theta0,alpha)     + ct0*s2a*alpha_coef(l,3,3,theta0,alpha);
  }
  else if (_k == 4){
    vals[0] *= 0.5*std::sqrt(3.0)*s2t0*ca * beta_coef(l,0,1,theta0,alpha);
    vals[1] *= 0.5*std::sqrt(3.0)*s2t0*ca * beta_coef(l,2,1,theta0,alpha);
    vals[2] *= 0.5*std::sqrt(3.0)*s2t0*sa * beta_bar_coef(l,2,1,theta0,alpha);
    vals[3] *= 0.5*std::sqrt(3.0)*st0*st0*s2a * beta_bar_coef(l,1,1,theta0,alpha);
    vals[4] *= 0.5*std::sqrt(3.0)*st0*st0*s2a * beta_bar_coef(l,3,1,theta0,alpha);
    vals[5] *= (1.5*ct0*ct0-0.5) * beta_coef(l,1,1,theta0,alpha);
    vals[6] *= std::sqrt(3.0)*0.5*st0*st0*c2a * beta_coef(l,1,1,theta0,alpha);
    vals[7] *= std::sqrt(3.0)*0.5*st0*st0*c2a * beta_coef(l,3,1,theta0,alpha);
  }
  return vals;
}

const double al00[] = {alpha_coef(0,0,0,0,0), 
                       alpha_coef(1,0,0,0,0), 
                       alpha_coef(2,0,0,0,0), 
                       alpha_coef(3,0,0,0,0), 
                       alpha_coef(4,0,0,0,0), 
                       alpha_coef(5,0,0,0,0), 
                       alpha_coef(6,0,0,0,0), 
                       alpha_coef(7,0,0,0,0), 
                       alpha_coef(8,0,0,0,0), 
                       alpha_coef(9,0,0,0,0),
                       alpha_coef(10,0,0,0,0)};
const double al11[] = {alpha_coef(0,1,1,0,0), 
                       alpha_coef(1,1,1,0,0), 
                       alpha_coef(2,1,1,0,0), 
                       alpha_coef(3,1,1,0,0), 
                       alpha_coef(4,1,1,0,0), 
                       alpha_coef(5,1,1,0,0), 
                       alpha_coef(6,1,1,0,0), 
                       alpha_coef(7,1,1,0,0), 
                       alpha_coef(8,1,1,0,0), 
                       alpha_coef(9,1,1,0,0),
                       alpha_coef(10,1,1,0,0)};
const double al22[] = {alpha_coef(0,2,2,0,0), 
                       alpha_coef(1,2,2,0,0), 
                       alpha_coef(2,2,2,0,0), 
                       alpha_coef(3,2,2,0,0), 
                       alpha_coef(4,2,2,0,0), 
                       alpha_coef(5,2,2,0,0), 
                       alpha_coef(6,2,2,0,0), 
                       alpha_coef(7,2,2,0,0), 
                       alpha_coef(8,2,2,0,0), 
                       alpha_coef(9,2,2,0,0),
                       alpha_coef(10,2,2,0,0)};
const double al33[] = {alpha_coef(0,3,3,0,0), 
                       alpha_coef(1,3,3,0,0), 
                       alpha_coef(2,3,3,0,0), 
                       alpha_coef(3,3,3,0,0), 
                       alpha_coef(4,3,3,0,0), 
                       alpha_coef(5,3,3,0,0), 
                       alpha_coef(6,3,3,0,0), 
                       alpha_coef(7,3,3,0,0), 
                       alpha_coef(8,3,3,0,0), 
                       alpha_coef(9,3,3,0,0),
                       alpha_coef(10,3,3,0,0)};
const double bbl22[] = {beta_bar_coef(0,2,2,0,0), 
                        beta_bar_coef(1,2,2,0,0), 
                        beta_bar_coef(2,2,2,0,0), 
                        beta_bar_coef(3,2,2,0,0), 
                        beta_bar_coef(4,2,2,0,0), 
                        beta_bar_coef(5,2,2,0,0), 
                        beta_bar_coef(6,2,2,0,0), 
                        beta_bar_coef(7,2,2,0,0), 
                        beta_bar_coef(8,2,2,0,0), 
                        beta_bar_coef(9,2,2,0,0),
                        beta_bar_coef(10,2,2,0,0)};
const double bbl11[] = {beta_bar_coef(0,1,1,0,0), 
                        beta_bar_coef(1,1,1,0,0), 
                        beta_bar_coef(2,1,1,0,0), 
                        beta_bar_coef(3,1,1,0,0), 
                        beta_bar_coef(4,1,1,0,0), 
                        beta_bar_coef(5,1,1,0,0), 
                        beta_bar_coef(6,1,1,0,0), 
                        beta_bar_coef(7,1,1,0,0), 
                        beta_bar_coef(8,1,1,0,0), 
                        beta_bar_coef(9,1,1,0,0),
                        beta_bar_coef(10,1,1,0,0)};
const double bbl33[] = {beta_bar_coef(0,3,3,0,0), 
                        beta_bar_coef(1,3,3,0,0), 
                        beta_bar_coef(2,3,3,0,0), 
                        beta_bar_coef(3,3,3,0,0), 
                        beta_bar_coef(4,3,3,0,0), 
                        beta_bar_coef(5,3,3,0,0), 
                        beta_bar_coef(6,3,3,0,0), 
                        beta_bar_coef(7,3,3,0,0), 
                        beta_bar_coef(8,3,3,0,0), 
                        beta_bar_coef(9,3,3,0,0),
                        beta_bar_coef(10,3,3,0,0)};

double calc_Intensity_s_orbital(double alpha_loc, 
                                double nu_in,
                                double t0, 
                                int l_max, 
                                int p_max, 
                                int n0,
                                int el_nr){
  double z_temp = 0.0;
  double de = 0.0;
  if (n0 == 1){
    z_temp = nu_in / (get_ionization_energy_1s(el_nr) / h);
    de = one_minus_delta_edge(el_nr,1,0,0);
  }
  else if (n0 == 2){
    z_temp = nu_in / (get_ionization_energy_2s(el_nr) / h);
    de = one_minus_delta_edge(el_nr,2,0,0);
  }
  else if (n0 == 3){
    z_temp = nu_in / (get_ionization_energy_3s(el_nr) / h);
    de = one_minus_delta_edge(el_nr,3,0,0);
  }
  if (z_temp < 1) return 0;
  z_temp *= de;

  double par = 0;
  double orth = 0;
  double temp;
  std::vector<int> fac;
  for (int p=0; p < int(p_max/2+1); p++)
    fac.push_back(p+1);
  for (int p=int(p_max/2+1); p>0; p++)
    fac.push_back(p+1);

  for (int l=0; l<=l_max; l+=1){
    temp = 0;
    for (int p=0; p<=p_max+1; p+=2){
      auto r = f_a_for_p(el_nr, l, 0, z_temp, nu_in, n0, p);
      for (int run=0; run < r.size(); run++)
        temp += real(r[run]);
    }
    par += alpha_coef(l,1,1,t0,alpha_loc) * temp;
    orth += beta_coef(l,1,1,t0,alpha_loc) * temp;
  }
  
  return (par*par + orth*orth);
}

double calc_Intensity_p_orbital(double alpha_loc,
                                double nu_in,
                                double t0,
                                int l_max,
                                int p_max,
                                int n0,
                                int subshell,
                                int el_nr){
  double z_temp = 0, de = 0;
  if (n0 == 1){
    printf("This shell doesn't have p orbitals");
    exit(-1);
  }
  else if (n0 == 2){
    if (subshell == 1){
      z_temp = nu_in / (get_ionization_energy_2p1_2(el_nr) / h);
      de = one_minus_delta_edge(el_nr,2,1,0.5);
    }
    else{
      z_temp = nu_in / (get_ionization_energy_2p3_2(el_nr) / h);
      de = one_minus_delta_edge(el_nr,2,1,1.5);
    }
  }
  else if (n0 == 3){
    if (subshell == 1){
      z_temp = nu_in / (get_ionization_energy_3p_1_2(el_nr) / h);
      de = one_minus_delta_edge(el_nr,3,1,0.5);
    }
    else{
      z_temp = nu_in / (get_ionization_energy_3p_3_2(el_nr) / h);
      de = one_minus_delta_edge(el_nr,3,1,1.5);
    }
  }
  if (z_temp < 1) return 0;
  z_temp *= de;

  double par = 0;
  double orth = 0;
  std::vector<int> fac;
  for (int p=0; p < int(p_max/2+1); p++)
    fac.push_back(p+1);
  for (int p=int(p_max/2+1); p>0; p++)
    fac.push_back(p+1);

  int ms[] = {0,1,2};
  for (int l=0; l<=l_max; l++){
    double mults[] = {al00[l],al11[l],al22[l]+bbl22[l]};
    double temp[] = {0,0,0};
    for (int p=0; p<=p_max; p+=2){
      for (int nummy=0; nummy < 3; nummy++){
        auto r = f_p_el_for_p(el_nr, l, ms[nummy], ms[nummy], z_temp, nu_in, n0, p);
        for (int run=0; run < r.size(); run++)
          temp[run] += real(r[run]) * mults[nummy];          
      }
      for (int nummy=0; nummy < 3; nummy++){
        par += apply_angle_part_s_parallel(temp[nummy], t0, alpha_loc, fac[nummy]);
        orth += apply_angle_part_s_orthogonal(temp[nummy], t0, alpha_loc, fac[nummy]);
      }
    }
  }
  
  //#Has division by 3**2 to avergae over all 3 p orbitals a.k.a. the p-wonder
  return (par*par + orth*orth) / 9.0;
}

double calc_Intensity_d_orbital(double alpha_loc, 
                                double nu_in, 
                                double t0, 
                                int l_max, 
                                int p_max, 
                                int n0, 
                                int subshell, 
                                int el_nr){
  double z_temp, de;
  if (n0 == 1){
    printf("This shell doesn't have d orbitals");
    exit(-1);
  }
  else if (n0 == 2){
    printf("This shell doesn't have d orbitals");
    exit(-1);
  }
  else if (n0 == 3){
    if (subshell == 1){
      z_temp = nu_in / (get_ionization_energy_3d_3_2(el_nr) / h);
      de = one_minus_delta_edge(el_nr,3,2,1.5);
    }
    else{
      z_temp = nu_in / (get_ionization_energy_3d_5_2(el_nr) / h);
      de = one_minus_delta_edge(el_nr,3,2,2.5);
    }
  }
  else{
    z_temp = 0;
    de = 0;
  }
  if (z_temp < 1) return 0;
  z_temp *= de;

  double par = 0;
  double orth = 0;

  std::vector<int> fac;
  for (int p=0; p < int(p_max/2+1); p++)
    fac.push_back(p+1);
  for (int p=int(p_max/2+1); p>0; p++)
    fac.push_back(p+1);

  int ms[] = {0,1,2,3,4};
  for (int l=0; l<=l_max; l++){
    double mults[] = {al00[l],bbl11[l]+al11[l],al22[l]+bbl22[l],bbl33[l]+al33[l],al11[l]};
    double temp[] = {0,0,0,0,0};
    for (int p=0; p<=p_max; p+=2){
      for (int nummy=0; nummy < 5; nummy++){
        auto r = f_d_el_for_p(el_nr, l, ms[nummy], ms[nummy], z_temp, nu_in, n0, p);
        for (int run=0; run < r.size(); run++)
          temp[run] += real(r[run]) * mults[nummy];          
      }
      for (int nummy=0; nummy < 5; nummy++){
        par += apply_angle_part_s_parallel(temp[nummy], t0, alpha_loc, fac[nummy]);
        orth += apply_angle_part_s_orthogonal(temp[nummy], t0, alpha_loc, fac[nummy]);
      }
    }
  }
  //#Has division by 5**2 to accomodate the averaging over 5 d orbtials  a.k.a. the d-wonder
  return (par*par + orth*orth) / 25.0;
}