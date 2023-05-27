#pragma once
#include "convenience.h"
#include "constants_and_atomic_properties.h"
#include "legendre_polynomials.h"
#include <vector>
#include <complex>
#include <functional>

M_type prepare_M(int p,int l,double z,int n_0){
  cdouble inprime(n_prime_from_z(z,n_0)*c_one);
  M_type M(p+2);
  for(int i=0; i<p+2; i++)
    M[i].resize(p+2,0.0);
  M[0][0] = 1.0;
  M[1][0] = inprime;
  M[1][1] = -2*(l+1);
  for (int j = 1; j <= p; j++){
    M[j+1][0] = -0.25*M[j][1] + inprime*M[j][0];
    for (int s = 1; s<j; s++)
      M[j+1][s] = -0.25*(s+1)*M[j][s+1] + inprime*M[j][s] + double(s-1-2*(l+j+1)) * M[j][s-1];
    M[j+1][j] = inprime * M[j][j] + double(j-1-2*(l+j+1)) * M[j][j-1];
    M[j+1][j+1] = double(j-2*(l+j+1)) * M[j][j];
  }
  return M;
}

cdouble g_from_M(M_type M, cdouble xhi, int j){
  cdouble sum = M[j][j];
  for (int i=j-1; i>-1; i--){
    sum *= xhi;
    sum += M[j][i];
  }
  return sum;
}

cdouble K_recursive_from_z(int p, int l, double b, double z, int n_0){
  if (l > p+1)
    return 0.0;
  cdouble zm1 = z-1;
  cdouble xhi = c_one/(2.0*std::sqrt(zm1));
  cdouble banane = g_from_M(prepare_M(p,l,z, n_0), xhi, p+1-l);
  cdouble part1 = -pow(2,p+3+l) * c_one * PI * pow(std::sqrt(zm1),p+2+l) / pow(-z,p+2) / pow(b*c_one,p+2-l);
  cdouble ex = exp_from_z(z,n_0);
  return part1 * banane * ex;
}

W_type make_matrix_W(int a,int c,int M){
  W_type Q(M+2);
  for(int i=0; i<M+2; i++)
    Q[i].resize(M+2,0.0);
  if (c == 0)
    Q[c][0] = 1;
  else if (c == 1)
    Q[c][0] = 1;
  else if (c == 2){
    Q[c][0] = 3;
    Q[c][2] = -3;
  }
  else if (c == 3){
    Q[c][0] = 15;
    Q[c][2] = -15;
  }
  else if (c == 4){
    Q[c][0] = 105;
    Q[c][2] = -210;
    Q[c][4] = 105;
  }
  for (int s=0; s<=M; s++)
    if (Q[c][s] != 0)
      Q[c+1][s+1] = (2*c+1) * Q[c][s];
  for (int k=c;k<M; k++){
    for (int s=0; s<=M+1; s++){
      if (Q[k][s] != 0){
        Q[k+2][s] = (-(k+1+c)/(k+2-c)) * Q[k][s];
      }
    }
    for (int s=0; s<=M; s++){
      if (Q[k+1][s] != 0){
        Q[k+2][s+1] += (2*k+3)/(k+2-c) * Q[k+1][s];
      }
    }
  }
  int h = int((a+c%2)/2);
  W_type W(Q);
  for (int run=0; run < h; run++)
    for (int k=c; k<M+2; k++)
      for (int s=0; s<M; s++)
        W[k][s+2] -= Q[k][s];
    Q = W;
  return W;
}

double value_from_W(W_type W, int b,int l,int M){
  double sum = 0;
  for (int s=0; s<M+1; s++){
    int x = b+s+1;
    if (x%2 != 0)
      sum += W[l][s] * 2/x;
  }
  return sum;
}

double J(int a,int b,int c,int l, W_type Matrix){
  int M = a+c+l;
  if (Matrix.size() == 0)
    Matrix = make_matrix_W(a,c, M);
  double result = value_from_W(Matrix,b,l,M);
  return result;
}

W_type W20 = make_matrix_W(2, 0, 20);
W_type W00 = make_matrix_W(0, 0, 20);
W_type W11 = make_matrix_W(1, 1, 20);
W_type W22 = make_matrix_W(2, 2, 20);
W_type W31 = make_matrix_W(3, 1, 20);
W_type W33 = make_matrix_W(3, 3, 20);

inline cdouble n1f(double nu,int p){
  return pow(std::complex<double> (0.0,-q(nu)),p) / double(ft[p]);
}

cdouble A_l_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J_ = J(1,p,1,l,W11);
  if (J_ == 0.0) return 0.0;
  cdouble K1 = K_recursive_from_z(p,l,b,z, n_0);
  if (K1 == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -b/2.0/pow(-2*b*std::sqrt(z-1),l+1);
  if (z<1) part1 *= std::sqrt(2.0);
  return part1 * n1 * J_ * K1;
}

cdouble C_l_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J_ = J(1,p,1,l,W11);
  if (J_ == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p, l, b, z, n_0);
  cdouble K2 = K_recursive_from_z(p+1, l, b, z, n_0);
  cdouble K2_mod = b/2*K2;
  cdouble Ktot = (K1-K2_mod);
  if (Ktot == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -b/pow(-2*b*std::sqrt(z-1),l+1);
  if (z<1) part1 *= std::sqrt(2.0);
  return part1 * n1 * J_ * Ktot;
}

cdouble E_l_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J_ = J(1,p,1,l,W11);
  if (J_ == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p, l, b, z, n_0);
  cdouble K2 = K_recursive_from_z(p+1, l, b, z, n_0);
  cdouble K3 = K_recursive_from_z(p+2, l, b, z, n_0);
  cdouble Ktot = (3.*K1-10.*b/3.*K2+2.*b*b/3.*K3);
  if (Ktot == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -b/2.0/pow(-2.0*b*std::sqrt(z-1),l+1);
  if (z<1) part1 *= std::sqrt(2.0);
  return part1 * n1 * J_ * Ktot;
}

cdouble B2_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J_ = J(2,p,2,l,W22);
  if (J_ == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p+1,l,b,z,n_0);
  if (K1 == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -b*b/(4.0*pow(-2.0*b*std::sqrt(z-1),l+1));
  if (z<1) part1 *= std::sqrt(2.0);
  return part1 * n1 * J_ * K1;
}

cdouble B1_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J_ = J(1,p+1,1,l,W11);
  if (J_ == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p+1,l,b,z,n_0);
  if (K1 == 0.) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -b*b/(2.0*pow(-2.0*b*std::sqrt(z-1),l+1));
  if (z<1) part1 *= std::sqrt(2.0);
  return part1 * n1 * J_ * K1;
}
  
cdouble B0_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J1 = J(2,p,0,l,W20);
  double J2 = J(0,p,0,l,W00);
  if (J1 == 0 && J2 == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p+1,l,b,z,n_0);
  cdouble K2 = K_recursive_from_z(p,l,b,z,n_0);
  cdouble tot_term = (-2. * J2 * K2 + b* J1 * K1);
  if (tot_term == 0.) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -b/(2.*pow(-2.*b*std::sqrt(z-1),l+1));
  if (z<1) part1 *= std::sqrt(2.);
  return part1 * n1 * tot_term;
}

cdouble D2_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J_ = J(2,p,2,l,W22);
  if (J_ == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p+1,l,b,z,n_0);
  cdouble K2 = K_recursive_from_z(p+2,l,b,z,n_0);
  cdouble Ktot = (3.0*K1 - b * K2);
  if (Ktot == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -sqrt(2./3.)*b*b/(4.0*pow(-2*b*std::sqrt(z-1),l+1));
  if (z<1) part1 *= sqrt(2.0);
  return part1 * n1 * J_ * Ktot;
}

cdouble D1_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J_ = J(1,p+1,1,l,W11);
  if (J_ == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p+1,l,b,z,n_0);
  cdouble K2 = K_recursive_from_z(p+2,l,b,z,n_0);
  cdouble Ktot = (3.0 * K1 - b * K2);
  if (Ktot == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -sqrt(2./3.)*b*b/(2.0*pow(-2.0*b*std::sqrt(z-1),l+1));
  if (z<1) part1 *= sqrt(2.0);
  return part1 * n1 * J_ * Ktot;
}
  
cdouble D0_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J1 = J(2,p,0,l,W20);
  double J2 = J(0,p,0,l,W00);
  if (J1 == 0 && J2 == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p,l,b,z,n_0);
  cdouble K2 = K_recursive_from_z(p+1,l,b,z,n_0);
  cdouble K3 = K_recursive_from_z(p+2,l,b,z,n_0);
  cdouble tot_term = (J2*(-4.0*K1 + 2*b*K2) + b * J1 * (3.0*K2-b*K3));
  if (tot_term == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -sqrt(2./3.)*b/(2*pow(-2*b*std::sqrt(z-1),l+1));
  if (z<1) part1 *= sqrt(2.0);
  return part1 * n1 * tot_term;
}

cdouble G0_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J1 = J(2,p+1,0,l,W20);
  double J2 = J(0,p+1,0,l,W00);
  if (J1 == 0 && J2 == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p+1,l,b,z,n_0);
  cdouble K2 = K_recursive_from_z(p+2,l,b,z,n_0);
  cdouble tot_term = (-2.0*J2*K1 + b*J1*K2);
  if (tot_term == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -sqrt(2./3.)*b*b/2./pow(-2*b*std::sqrt(z-1),l+1);
  if (z<1) part1 *= sqrt(2.0);
  return part1 * n1 * tot_term;
}

cdouble G1_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J1 = J(1,p,1,l,W11);
  double J2 = J(3,p,1,l,W31);
  if (J1 == 0 && J2 == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p+1,l,b,z,n_0);
  cdouble K2 = K_recursive_from_z(p+2,l,b,z,n_0);
  cdouble tot_term = (-J1*K1 + b/4.*J2*K2);
  if (tot_term == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -sqrt(2./3.)*b*b/2./pow(-2*b*std::sqrt(z-1),l+1);
  if (z<1) part1 *= sqrt(2.0);
  return part1 * n1 * tot_term;
}

cdouble G2_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J1 = J(2,p+1,2,l,W22);
  if (J1 == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p+2,l,b,z,n_0);
  if (K1 == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -sqrt(2./3.)*b*b*b/4./pow(-2*b*std::sqrt(z-1),l+1);
  if (z<1) part1 *= sqrt(2.0);
  return part1 * n1 * (J1*K1);
}

cdouble G3_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){
  double J1 = J(3,p,3,l,W33);
  if (J1 == 0) return 0.0;
  cdouble K1 = K_recursive_from_z(p+2,l,b,z,n_0);
  if (K1 == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -sqrt(2./3.)*b*b*b/8./pow(-2*b*std::sqrt(z-1),l+1);
  if (z<1) part1 *= sqrt(2.0);
  return part1 * n1 * (J1*K1);
}

cdouble G4_from_z_for_p(double b, double z, int n_0, int l, double nu, int p){ //'This is G_tilde'
  double J1 = J(1,p,1,l,W11);
  double J2 = J(1,p+2,1,l,W11);
  if (J1 == 0 && J2 == 0.0) return 0.0;
  cdouble K1 = K_recursive_from_z(p+1,l,b,z,n_0);
  cdouble K2 = K_recursive_from_z(p+2,l,b,z,n_0);
  cdouble tot_term = (1./3.*J1*(2.0*K1-b*K2) + b*J2*K2);
  if (tot_term == 0.0) return 0.0;
  cdouble n1 = n1f(nu,p);
  cdouble part1 = -sqrt(1./2.)*b*b/2./pow(-2*b*std::sqrt(z-1),l+1);
  if (z<1) part1 *= sqrt(2.0);
  return part1 * n1 * tot_term;
}

typedef std::function<cdouble(double, double, int, int, double, int)> matrix_func;

matrix_func function_selector(int n_0, int l_0, int g_m){
  matrix_func func;
  if (l_0 == 0){
    if (n_0 == 1)
      func = A_l_from_z_for_p;
    else if (n_0 == 2)
      func = C_l_from_z_for_p;
    else if (n_0 == 3)
      func = E_l_from_z_for_p;
    else
      exit(-1);
  }
  else if (l_0 == 1){
    if (n_0 == 2){
      if (g_m == 0)
        func = B0_from_z_for_p;
      else if (g_m == 1)
        func = B1_from_z_for_p;
      else if (g_m == 2)
        func = B2_from_z_for_p;
      else exit(-1);
    }
    else if (n_0 == 3){
      if (g_m == 0)
        func = D0_from_z_for_p;
      else if (g_m == 1)
        func = D1_from_z_for_p;
      else if (g_m == 2)
        func = D2_from_z_for_p;
      else exit(-1);
    }
    else
      exit(-1);
  }
  else if (l_0 == 2){
    if (n_0 == 3){
      if (g_m == 0)
        func = G0_from_z_for_p;
      else if (g_m == 1)
        func = G1_from_z_for_p;
      else if (g_m == 2)
        func = G2_from_z_for_p;
      else if (g_m == 3)
        func = G3_from_z_for_p;
      else if (g_m == 4)
        func = G4_from_z_for_p;
      else exit(-1);
    }
    else exit(-1);
  }
  else exit(-1);
  return func;
}

std::vector<cdouble> f_a_for_p(int Z, int l,int k,double z,double nu_in,int n_0,int p){
  std::vector<cdouble> result(2*p+1,0.);
  if (z <= 0) return result;
  double b_ = b(n_0,0,Z);
  //prefactor = N0_square(b_) * N_square(l,1,b_,n_0,z)
  cdouble prefactor = 1.0;
  matrix_func func = function_selector(n_0, 0,0);
  for (int j=0; j <= p; j++){
    cdouble matrix_value1 = func(b_,z,n_0,l,nu_in,j);
    if (matrix_value1 == 0.0){
      result[j] = 0.0;
      continue;
    }
    cdouble matrix_value2 = func(b_,z,n_0,l,nu_in,p-j);
    if (matrix_value2 == 0.0){
      result[j] = 0.0;
      continue;
    }
    cdouble postfactor = matrix_value1 * std::conj(matrix_value2);
    result[j] = prefactor * postfactor;
  }
  return result;
}

std::vector<cdouble> f_p_el_for_p(int Z,int l,int g_m,int g_k,double z,double nu_in,int n_0,int p){
  std::vector<cdouble> result(2*p+1,0.);
  if (z <= 0) return result;
  double b_ = b(n_0, 1, Z);
  //prefactor = N0_square(b_) * N_square(l,g_m,b_,n_0,z)
  cdouble prefactor = 1.0;
  matrix_func func = function_selector(n_0,1,g_m);
  matrix_func conjugate_function = function_selector(n_0,1,g_k);
  for (int j=0; j <= p; j++){
    cdouble matrix_value1 = func(b_,z,n_0,l,nu_in,j);
    if (matrix_value1 == 0.0){
      result[j] = 0.0;
      continue;
    }
    cdouble matrix_value2 = func(b_,z,n_0,l,nu_in,p-j);
    if (matrix_value2 == 0.0){
      result[j] = 0.0;
      continue;
    }
    cdouble postfactor = matrix_value1 * std::conj(matrix_value2);
    result[j] = prefactor * postfactor;
  }
  return result;
}

std::vector<cdouble> f_d_el_for_p(int Z,int l,int g_m,int g_k,double z,double nu_in,int n_0,int p){
  std::vector<cdouble> result(2*p+1,0.);
  if (z <= 0) return result;
  double b_ = b(n_0, 1, Z);
  cdouble prefactor = 1.0;
  //if g_m <= 3: prefactor = N0_square(b_) * N_square(l,g_m,b_,n_0,z)
  //else: prefactor = N0_square(b_) * N_square(l,1,b_,n_0,z)
  matrix_func func = function_selector(n_0,2,g_m);
  matrix_func conjugate_function = function_selector(n_0,2,g_k);
  for (int j=0; j <= p; j++){
    cdouble matrix_value1 = func(b_,z,n_0,l,nu_in,j);
    if (matrix_value1 == 0.0){
      result[j] = 0.0;
      continue;
    }
    cdouble matrix_value2 = func(b_,z,n_0,l,nu_in,p-j);
    if (matrix_value2 == 0.0){
      result[j] = 0.0;
      continue;
    }
    cdouble postfactor = matrix_value1 * std::conj(matrix_value2);
    result[j] = prefactor * postfactor;
  }
  return result;
}
