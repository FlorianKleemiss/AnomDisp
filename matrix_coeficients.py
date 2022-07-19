from distutils.errors import LibError
import math
import cmath
from ssl import SSL_ERROR_WANT_X509_LOOKUP
from turtle import speed
import numpy as np
import scipy.special
import legendre_plynomials
from constants_and_atomic_properties import *


def b(n_0, l_0, Z):
  factor = None
  Z_eff = None
  if n_0 == 1:
    factor = 1
    Z_eff = get_Zeff_1s(Z)
  elif n_0 == 2:
    factor = 2
    if l_0 == 0:
      Z_eff = get_Zeff_2s(Z)
    elif l_0 == 1:
      Z_eff = get_Zeff_2p_1_2(Z)
    elif l_0 == 2:
      Z_eff = get_Zeff_2p_3_2(Z)
  elif n_0 == 3:
    factor = 3
    if l_0 ==0:
      Z_eff = get_Zeff_3s(Z)
    elif l_0 == 1:
      Z_eff = get_Zeff_3p_1_2(Z)
    elif l_0 == 2:
      Z_eff = get_Zeff_3d_3_2(Z)
  return Z_eff/(factor*a0)

def k(E):
  return 2*math.pi / h * math.sqrt(2*el_mass*E)

def n_prime(E, Z, n,l):
  return b(n,l,Z)/k(E)

def N0(b_):
  return math.sqrt(pow(b_,3)/math.pi)

def product_n_prime(b_,k_,l):
  n_p = 2*b_/k_
  fact = 1.0
  for nu in range(1,l+1):
    #print(nu)
    fact *= n_p * n_p + nu*nu
  denom = 1-math.exp(-2*math.pi*n_p)
  return fact/denom
  
def N(l, m, b_, k_):
  result = 2*(2*l+1)*math.factorial(l-m)/math.factorial(l+m) * 4 * math.pi * el_mass/h/h * b_ * product_n_prime(b_,k_,l)
  if m >= 1:
    result *= 2
  return math.sqrt(result)
  
def q(nu):
  return 2*math.pi*nu/speed_of_light

def delta(a,b):
  if a == b:
      return 1
  else:
      return 0

######################### BEGIN OF MAKING K_p,l #######################################################
def C(n_prime, p, l, s):
  coef = scipy.special.binom(p+1-l,s)
  product1 = 1
  UL1 = s+1
  for t1 in range(1,s+1):
    print("t1: %d"%t1)
    product1 *= -complex(l+t1,n_prime)
  product2 = 1
  UL2 = p+1-l-s
  for t2 in range(1,p+2-l-s):
    print("t2: %d"%t2)
    product2 *= complex(-l-t2,n_prime)
  return coef * product1 * product2

def K(p, l, k, b):
  n_prime = 2*b/k
  exponential_part = cmath.exp(-2.0*n_prime * cmath.atan(1.0/(2.0*n_prime)))
  numerator = -exponential_part*complex(0.0, math.pi) * pow(2.0,p+3*l+5)
  denominator = pow(complex(0.0,k),p+2-l) * pow(complex(-2.0,n_prime),p+2) * pow(complex(2.0,n_prime),l+1)
  factor1 = numerator / denominator
  sum = 0
  factor2 = complex(-2.0,n_prime)/complex(2.0,n_prime)
  for s in range(p+2-l):
    C_res = C(n_prime, p, l, s)
    sum += C_res * pow(factor2,s)
  return factor1 * sum

def prepare_M(p,l,np4):
  np44 = 4*np4
  M = np.zeros((p+2,p+2),dtype=complex)
  M.fill(complex(0.0,0.0))
  M[0,0] = complex(1.0,0)
  M[1,0] = complex(0.0,np44)
  M[1,1] = -2*(l+1)
  for j in range(2,p+2):
    M[j,0] = -0.25*M[j-1,1] + complex(0,np44)*M[j-1,0]
    for s in range(1,j):
      M[j,s] = -0.25*(s+1)*M[j-1,s+1] + complex(0,np44)*M[j-1,s] + (s-1-2*(l+j))*M[j-1,s-1]
    M[j,j-1] = complex(0,np44) * M[j-1,j-1] + (j-2-2*(l+j))*M[j-1,j-2]
    M[j,j] = (j-1-2*(l+j))*M[j-1,j-1]
  return M

def g_from_M(M, zeta, j):
  sum = M[j,j]
  for i in range(j-1,-1,-1):
    sum *= zeta
    sum += M[j,i]
  return sum
  
def K_recursive(p, l, k, b):
  if l > p+1:
    return 0
  np4 = b/(k*2) #n'/4
  M = prepare_M(p,l,np4)
  part1 = -complex(0,2*math.pi) / pow(complex(0,2*k),p+2-l)
  g = g_from_M(M, complex(0,np4), p+1-l)
  part2 = g / pow(-pow(np4,2)-0.25,p+2)
  ex = math.exp(-8*np4 * math.atan(1/(8*np4)))
  return part1 * part2 * ex

################### END OF CALC K ############################################
################### BEGIN CALC J  ############################################

def make_matrix_W(a,c, M):
  Q = np.zeros((M+2,M+2))
  W = None
  #Q[0,0] = 1
  if c == 0:
    Q[c,0] = 1
  elif c == 1:
    Q[c,0] = 1
  elif c == 2:
    Q[c,0] = 3
    Q[c,2] = -3
  for s in range(0,M+1):
    if Q[c,s] != 0:
      Q[c+1,s+1] = (2*c+1) * Q[c,s]
  for k in range(c,M):
    for s in range(0,M+2):
      if Q[k,s] != 0:
        fact = (-(k+1+c)/(k+2-c))
        Q[k+2,s] = fact * Q[k,s]
    for s in range(0,M+1):
      if Q[k+1,s] != 0:
        fact = (2*k+3)/(k+2-c)
        Q[k+2,s+1] += fact * Q[k+1,s]
  h = int((a+c%2)/2)
  W = np.array(Q, copy=True)
  for run in range(0,h):
    for k in range(c,M+2):
      for s in range(0,M):
        W[k,s+2] -= Q[k,s]
    Q = np.array(W, copy=True)
  return W

def value_from_W(W,b,l,M):
  sum = 0
  for s in range(0,M+1):
    x = b+s+1
    if x%2 != 0:
      sum += W[l,s] * 2/x
  return sum

def J(a,b,c,l, Matrix = None):
  M = a+c+20
  if Matrix is None:
    Matrix = make_matrix_W(a,c, M)
  result = value_from_W(Matrix,b,l,M)
  return result

def J0_hat(p,l):
  if p == 0:
    return 1/3 * delta(l,0) - 1/15 * delta(l,2)
  elif p == 1:
    return 1/15 * delta(l,1) - 1/35 * delta(l,3)
  elif p == 2:
    return 1/15 * delta(l,0) + 1/105 * delta(l,2) - 4/315 * delta(l,4)
  elif p ==3:
    return 1/35 * delta(l,1) - 1/315 * delta(l,3) - 4/693 * delta(l,5)
  else:
    return 0.25 * J(2,p,0,l)
    
def J0_bar(p,l):
  if p == 0:
    return 2 * delta(l,0)
  elif p == 1:
    return 2/3 * delta(l,1)
  elif p == 2:
    return 2/3 * delta(l,0) + 4/15 * delta(l,2)
  elif p == 3:
    return 2/5 * delta(l,1) + 4/35 * delta(l,3)
  elif p == 4:
    return 2/5 * delta(l,0) + 8/35 * delta(l,2) + 16/315 * delta(l,4)
  elif p == 5:
    return 2/7 * delta(l,1) + 8/63 * delta(l,3) + 16/693 * delta(l,5)
  else:
    return J(0,p,0,l)

def J1(p,l):
  if p == 0:
    return 4/3 * delta(l,1)
  elif p == 1:
    return 4/5 * delta(l,2)
  elif p == 2:
    return 4/15 * delta(l,1) + 16/35 * delta(l,3)
  elif p == 3:
    return 12/35 * delta(l,2) + 16/63 * delta(l,4)
  elif p == 4:
    return 4/35 * delta(l,1) + 32/105 * delta(l,3) + 32/231 * delta(l,5)
  else:
    return J(1,p,1,l)
  
def J2(p,l):
  if p == 0:
    return 16/5 * delta(l,2)
  elif p == 1:
    return 16/7 * delta(l,3)
  elif p == 2:
    return 16/35 * delta(l,2) + 32/21 * delta(l,4)
  elif p == 3:
    return 16/21 * delta(l,3) + 32/33 * delta(l,5)
  else:
    return J(2,p,2,l)
  
################################# END OF CALC Js ###################################    
  
def polynom(k, eta):
  if k == 0: return 0
  elif k == 1: return math.sqrt(1-eta*eta)
  else: return (2*k-1)/(k-1) * eta * polynom(k-1,eta) - k/(k-1) * polynom(k-2,eta)

def C_1(b_, k_, l, nu, p_limit):
  #print("*"*30 + "Start C_1" + "*"*30)
  part1 = b_/pow(-2*k_,l+1)
  sum = 0
  for p in range(p_limit):
    q_l = q(nu)
    n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
    J = J1(p,l)
    K1 = K_recursive(p,l,k_,b_)
    K2 = K_recursive(p+1,l,k_,b_)
    K2_mod = b_/2*K2
    #print("p: {:d} n1: {:.5e}".format(p,n1))
    #print("  J: {:2f} K1: {:.5e} K2: {:.5e}".format(J,K1,K2))
    #print("  ->new_contribution:{:12.5e}+{:12.5e}j old value:{:.5e}<-".format((n1 * J * (K1-K2_mod)).real,
    #                                                                          (n1 * J * (K1-K2_mod)).imag,
    #                                                                          sum))
    sum += n1 * J * (K1-K2_mod)
  return part1 * sum
  
def B_2(b_, k_, l, nu, p_limit):
  #print("*"*30 + "Start B_2" + "*"*30)
  part1 = b_*b_/(4*pow(-2*k_,l+1))
  sum = 0
  for p in range(p_limit):
    n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
    J = J2(p,l)
    K1 = K_recursive(p+1,l,k_,b_)
    #print("p: {:d} n1: {:.5e}".format(p,n1))
    #print("  J: {:2f} K1: {:.5e}".format(J,K1))
    #print("  ->new_contribution:{:12.5e}+{:12.5e}j old value:{:.5e}<-".format((n1 * J * K1).real,
    #                                                                          (n1 * J * K1).imag,
    #                                                                          sum))
    sum += n1 * J * K1
  return part1 * sum

def B_1(b_, k_, l, nu, p_limit):
  #print("*"*30 + "Start B_1" + "*"*30)
  part1 = b_*b_/(2*pow(-2*k_,l+1))
  sum = 0
  for p in range(p_limit):
    n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
    J = J1(p+1,l)
    K1 = K_recursive(p+1,l,k_,b_)
    #print("p: {:d} n1: {:.5e}".format(p,n1))
    #print("  J: {:2f} K1: {:.5e}".format(J,K1))
    #print("  ->new_contribution:{:12.5e}+{:12.5e}j old value:{:.5e}<-".format((n1 * J * K1).real,
    #                                                                          (n1 * J * K1).imag,
    #                                                                          sum))
    sum += n1 * J * K1
  return part1 * sum
  
def B_0(b_, k_, l, nu, p_limit):
  #print("*"*30 + "Start B_0" + "*"*30)
  part1 = b_/pow(-2*k_,l+1)
  sum = 0
  for p in range(p_limit):
    n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
    J1 = J0_hat(p,l)
    K1 = K_recursive(p+1,l,k_,b_)
    J2 = J0_bar(p,l)
    K2 = K_recursive(p,l,k_,b_)
    #print("p: {:d} n1: {:.5e}".format(p,n1))
    #print("  J1: {:2f} K1: {:.5e} J2: {:2f} K2: {:.5e}".format(J1,K1,J2,K2))
    #print("  ->new_contribution:{:12.5e}+{:12.5e}j old value:{:.5e}<-".format((n1 * (2*b_*J1 * K1 - J2 * K2)).real,
    #                                                                          (n1 * (2*b_*J1 * K1 - J2 * K2)).imag,
    #                                                                          sum))
    sum += n1 * (2*b_*J1 * K1 - J2 * K2)
  return part1 * sum

def A(b_, k_, l, nu, p_limit):
  part1 = b_/pow(-2*k_,l+1)
  sum = 0
  for p in range(p_limit):
    n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
    J1 = J1(p,l)
    K1 = K_recursive(p+1,l,k_,b_)
    sum += n1 * J1 * K1

#Start of final matrixelement calculation
  
def A_a_1(Z, l, E, nu, p_limit):
  b_ = b(2, 0, Z)
  k_ = k(E)
  N_zero = N0(b_)
  Nlm = N(l, 0, b_, k_)
  C1 = C_1(b_, k_, l, nu, p_limit)
  #print("N_0: {:.2e}  Nlm: {:.2e}   C1:{:.2e}".format(N_zero,Nlm,C1))
  return N_zero * Nlm * C1

def A_b_1(Z, l, E, nu, p_limit):
  b_ = b(2, 1, Z)
  k_ = k(E)
  N_zero = N0(b_)
  Nlm = N(l, 1, b_, k_)
  B1 = B_1(b_, k_, l, nu, p_limit)
  #print("N_0: {:.2e}  Nlm: {:.2e}   B1:{:.2e}".format(N_zero,Nlm,B1))
  return N_zero * Nlm * B1

def A_c_0(Z, l, E, nu, p_limit):
  b_ = b(2, 2, Z)
  k_ = k(E)
  N_zero = N0(b_)
  Nlm = N(l, 0, b_, k_)
  B0 = B_0(b_, k_, l, nu, p_limit)
  #print("N_0: {:.2e}  Nlm: {:.2e}   B0:{:.2e}".format(N_zero,Nlm,B0))
  return N_zero * Nlm * B0

def A_c_2(Z, l, E, nu, p_limit):
  b_ = b(2, 2, Z)
  k_ = k(E)
  N_zero = N0(b_)
  Nlm = N(l, 2, b_, k_)
  B2 = B_2(b_, k_, l, nu, p_limit)
  #print("N_0: {:.2e}  Nlm: {:.2e}   B2:{:.2e}".format(N_zero,Nlm,B2))
  return N_zero * Nlm * B2

def A_d_2(Z, l, E, nu, p_limit):
  b_ = b(2, 2, Z)
  k_ = k(E)
  N_zero = N0(b_)
  Nlm = N(l, 2, b_, k_)
  B2 = B_2(b_, k_, l, nu, p_limit)
  #print("N_0: {:.2e}  Nlm: {:.2e}   B0:{:.2e}".format(N_zero,Nlm,B2))
  return N_zero*Nlm * B2

#Start of calculation of products of matrix elements

def A_product(z,l,E,nu,p_limit, theta0, alpha):
  b_ = b(1, l, z)
  k_ = k(E)
  alpha_l00 = legendre_plynomials.alpha_coef(l,0,0,theta0,alpha)
  beta_l00 = legendre_plynomials.beta_coef(l,0,0,theta0,alpha)
  N2N2A = pow(N0(b_),2) * pow(N(l, 0, b_, k_),2) * A(b_, k_, l, nu, p_limit)
  return (N2N2A*alpha_l00,N2N2A*beta_l00)

def A_a_1_product(z,l,E,nu,p_limit, theta0, alpha):
  alpha_l11 = legendre_plynomials.alpha_coef(l,1,1,theta0,alpha)
  beta_l11 = legendre_plynomials.beta_coef(l,1,1,theta0,alpha)
  NNC = pow(abs(A_a_1(z,l,E,nu,p_limit)),2)
  return (NNC*alpha_l11,NNC*beta_l11)

def A_b_1_product(z,l,E,nu,p_limit, theta0, alpha):
  alpha_l21 = legendre_plynomials.alpha_coef(l,1,2,theta0,alpha)
  alpha_l11 = legendre_plynomials.alpha_coef(l,1,1,theta0,alpha)
  alpha_l01 = legendre_plynomials.alpha_coef(l,1,0,theta0,alpha)
  beta_l11 = legendre_plynomials.beta_coef(l,1,1,theta0,alpha)
  beta_l21 = legendre_plynomials.beta_coef(l,1,2,theta0,alpha)
  b_ = b(2, 1, z)
  k_ = k(E)
  N2N2B1 = pow(N0(b_),2) * pow(N(l, 1, b_, k_),2) * B_1(b_, k_, l, nu, p_limit)
  Bl0star = B_0(b_,k_,l,nu,p_limit).conjugate()
  Bl1star = B_1(b_,k_,l,nu,p_limit).conjugate()
  Bl2star = B_2(b_,k_,l,nu,p_limit).conjugate()
  st0 = np.sin(theta0)
  ct0 = np.cos(theta0)
  parallel = N2N2B1*((alpha_l01*Bl0star+alpha_l21*Bl2star)*st0+alpha_l11*Bl1star*ct0)
  orthogonal = N2N2B1*(beta_l11*Bl1star*ct0+beta_l21*Bl2star*st0)
  return (parallel,orthogonal)

def A_c_0_product(z,l,E,nu,p_limit, theta0, alpha):
  alpha_l20 = legendre_plynomials.alpha_coef(l,0,2,theta0,alpha)
  alpha_l10 = legendre_plynomials.alpha_coef(l,0,1,theta0,alpha)
  alpha_l00 = legendre_plynomials.alpha_coef(l,0,0,theta0,alpha)
  beta_l10 = legendre_plynomials.beta_coef(l,0,1,theta0,alpha)
  beta_l20 = legendre_plynomials.beta_coef(l,0,2,theta0,alpha)
  b_ = b(2, 1, z)
  k_ = k(E)
  N2N2B0 = pow(N0(b_),2) * pow(N(l, 0, b_, k_),2) * B_0(b_, k_, l, nu, p_limit)
  Bl0star = B_0(b_,k_,l,nu,p_limit).conjugate()
  Bl1star = B_1(b_,k_,l,nu,p_limit).conjugate()
  Bl2star = B_2(b_,k_,l,nu,p_limit).conjugate()
  st0 = np.sin(theta0)
  ct0 = np.cos(theta0)
  sa = np.sin(alpha)
  ca = np.cos(alpha)
  parallel = N2N2B0*(-alpha_l00*Bl0star*ct0*ca+alpha_l10*Bl1star*st0*ca+Bl2star*(beta_l20*sa-alpha_l20*ct0*ca))
  orthogonal = N2N2B0*(Bl0star*alpha_l00*sa+Bl1star*beta_l10*st0*ca - Bl2star*(alpha_l20*sa+beta_l20*ct0*ca))
  return (parallel,orthogonal)

def A_c_2_product(z,l,E,nu,p_limit, theta0, alpha):
  alpha_l22 = legendre_plynomials.alpha_coef(l,2,2,theta0,alpha)
  alpha_l12 = legendre_plynomials.alpha_coef(l,2,1,theta0,alpha)
  alpha_l02 = legendre_plynomials.alpha_coef(l,2,0,theta0,alpha)
  beta_l12 =  legendre_plynomials.beta_coef( l,2,1,theta0,alpha)
  beta_l22 =  legendre_plynomials.beta_coef( l,2,2,theta0,alpha)
  b_ = b(2, 1, z)
  k_ = k(E)
  N2N2B2 = pow(N0(b_),2) * pow(N(l, 0, b_, k_),2) * B_2(b_, k_, l, nu, p_limit)
  Bl0star = B_0(b_,k_,l,nu,p_limit).conjugate()
  Bl1star = B_1(b_,k_,l,nu,p_limit).conjugate()
  Bl2star = B_2(b_,k_,l,nu,p_limit).conjugate()
  st0 = np.sin(theta0)
  ct0 = np.cos(theta0)
  sa = np.sin(alpha)
  ca = np.cos(alpha)
  parallel   = N2N2B2*(-alpha_l02*Bl0star*ct0*ca+alpha_l12*Bl1star*st0*ca+Bl2star*(beta_l22*sa-alpha_l22*ct0*ca))
  orthogonal = N2N2B2*(Bl0star*alpha_l02*sa+Bl1star*beta_l12*st0*ca - Bl2star*(alpha_l22*sa+beta_l22*ct0*ca))
  return (parallel,orthogonal)

def A_d_2_product(z,l,E,nu,p_limit, theta0, alpha):
  alpha_l22 = legendre_plynomials.alpha_bar_coef(l,2,2,theta0,alpha)
  alpha_l12 = legendre_plynomials.alpha_bar_coef(l,2,1,theta0,alpha)
  alpha_l02 = legendre_plynomials.alpha_bar_coef(l,2,0,theta0,alpha)
  beta_l12 =  legendre_plynomials.beta_bar_coef( l,2,1,theta0,alpha)
  beta_l22 =  legendre_plynomials.beta_bar_coef( l,2,2,theta0,alpha)
  b_ = b(2, 2, z)
  k_ = k(E)
  N2N2B2 = pow(N0(b_),2) * pow(N(l, 0, b_, k_),2) * B_2(b_, k_, l, nu, p_limit)
  Bl0star = B_0(b_,k_,l,nu,p_limit).conjugate()
  Bl1star = B_1(b_,k_,l,nu,p_limit).conjugate()
  Bl2star = B_2(b_,k_,l,nu,p_limit).conjugate()
  st0 = np.sin(theta0)
  ct0 = np.cos(theta0)
  sa = np.sin(alpha)
  ca = np.cos(alpha)
  parallel   = N2N2B2*(-alpha_l02*Bl0star*ct0*sa+alpha_l12*Bl1star*st0*sa+Bl2star*(beta_l22*ca-alpha_l22*ct0*sa))
  orthogonal = N2N2B2*(Bl0star*alpha_l02*ca+Bl1star*beta_l12*st0*sa - Bl2star*(alpha_l22*ca+beta_l22*ct0*sa))
  return (parallel,orthogonal)

################## END of matrix element calculation

def print_Js():
  p_limit = 6
  l_limit = 6
  ##Js = np.zeros((4,l_limit,p_limit))
  ##for l in range(l_limit):
  ##  for p in range(p_limit):
  ##    try:
  ##      Js[0,l,p] = J1(p,l)
  ##    except:
  ##      Js[0,l,p] = 0.12345678
  ##    try:
  ##      Js[1,l,p] = J2(p,l)
  ##    except:
  ##      Js[1,l,p] = 0.12345678
  ##    try:
  ##      Js[2,l,p] = J0_hat(p,l)
  ##    except:
  ##      Js[2,l,p] = 0.12345678
  ##    try:
  ##      Js[3,l,p] = J0_bar(p,l)
  ##    except:
  ##      Js[3,l,p] = 0.12345678
  ##string = "   p= "
  ##for p in range(p_limit):
  ##  string += "{:8d}".format(p)
  ##string += "\n"
  ##for l in range(l_limit):
  ##  string += "l = {:4d}: ".format(l)
  ##  for type in range(4):
  ##    for p in range(p_limit):
  ##      if Js[type,l,p] != 0.12345678:
  ##        string += "{:+7.3f} ".format(Js[type,l,p])
  ##      else:
  ##        string += " ------ "
  ##    if type == 0:
  ##      string += " J1"
  ##    elif type == 1:
  ##      string += " J2"
  ##    elif type == 2:
  ##      string += " J0_hat"
  ##    elif type == 3:
  ##      string += " J0_bar"
  ##    string += "\n"
  ##    if type != 3:
  ##      string += "          "
  ##print(string)
  ##print('***'*30)
  #HERE WE USE THE GENERAL EQUATION
  Js = np.zeros((4,l_limit,p_limit))
  M = 2+2+20
  Matrix00 = make_matrix_W(0,0, M)
  Matrix11 = make_matrix_W(1,1, M)
  Matrix22 = make_matrix_W(2,2, M)
  Matrix20 = make_matrix_W(2,0, M)
  for l in range(l_limit):
    for p in range(p_limit):
      Js[0,l,p] = J(1,p,1,l, Matrix11)
      Js[1,l,p] = J(2,p,2,l, Matrix22)
      Js[2,l,p] = 0.25 * J(2,p,0,l, Matrix20)
      Js[3,l,p] = J(0,p,0,l, Matrix00)
  string = "   p= "
  for p in range(p_limit):
    string += "{:8d}".format(p)
  string += "\n"
  for l in range(l_limit):
    string += "l = {:4d}: ".format(l)
    for type in range(4):
      for p in range(p_limit):
        string += "{:+7.3f} ".format(Js[type,l,p])
      if type == 0:
        string += " J(1,p,1,l) = J1"
      elif type == 1:
        string += " J(2,p,2,l) = J2"
      elif type == 2:
        string += " J(2,p,0,l) *0.25 = J0_hat"
      elif type == 3:
        string += " J(0,p,0,l) = J0_bar"
      string += "\n"
      if type != 3:
        string += "          "
    string += "\n"
  print(string)
      
def compare_Ks():
  p = 3
  l = 2
  b_ = b(1, 2, 42)
  #b_ = 20000
  k_ = k(200000000000000*Ryd_ener/64)
  #k_ = 4*b_
  np4 = b_/(k_*2)
  k1 = K(p, l, k_, b_)
  k2 = K_recursive(p, l, k_, b_)
  print("k1: {:.8e} k2: {:.8e}".format(k1,k2))

def test():      
  Z = 42
  n = 10
  #l = 2
  
  #E = 200000*Ryd_ener/n/n
  Ek = get_ionization_energy_1s(Z) #K-Absorption edge of Mo
  Ed = get_ionization_energy_2s(Z) #L3-Absorption edge of Mo
  Ec = get_ionization_energy_2p1_2(Z) #L2-Absorption edge of Mo
  Eb = get_ionization_energy_2p3_2(Z) #L1-Absorption edge of Mo
  p_limit = 5
  minimal = 0.99*Eb/h
  maximal = 1.2*Ed/h
  stepsize = 1.0E15
  l_start = 1
  l_end = 8
  import matplotlib.pyplot as plt
  fig, axs = plt.subplots(l_end-l_start,2,sharex=True)
  E=20 #this is E_gamma' in EM
  theta0 = math.pi/3
  alp = math.pi/1.823
  factor = h/(4*math.pi*math.pi*el_mass)
  for l in range(l_start,l_end):
    x = []
    ya = []
    yb = []
    yc = []
    yd = []
    yai = []
    ybi = []
    yci = []
    ydi = []
    for step in range(int(minimal), int(maximal), int(stepsize)):
      x.append(step)
      Aa  = A_a_1_product(Z, l, E, step, p_limit, theta0, alp)
      Ab  = A_b_1_product(Z, l, E, step, p_limit, theta0, alp)
      Ac  = A_c_0_product(Z, l, E, step, p_limit, theta0, alp)
      Ac2 = A_c_2_product(Z, l, E, step, p_limit, theta0, alp)
      Ad  = A_d_2_product(Z, l, E, step, p_limit, theta0, alp)
      ya.append( (factor/(((Eb+E)/h)-step))*(Aa[0].real+Aa[1].real)/2)
      yb.append( (factor/(((Ec+E)/h)-step))*(Ab[0].real+Ab[1].real)/2)
      yc.append( (factor/(((Ed+E)/h)-step))*(Ac[0].real+Ac[1].real+Ac2[0].real+Ac2[1].real)/2)
      yd.append( (factor/(((Ed+E)/h)-step))*(Ad[0].real+Ad[1].real)/2)
  
      yai.append( (factor/(((Ed+E)/h)-step))*(Aa[0].imag+Aa[1].imag)/2)
      ybi.append( (factor/(((Ec+E)/h)-step))*(Ab[0].imag+Ab[1].imag)/2)
      yci.append( (factor/(((Ed+E)/h)-step))*(Ac[0].imag+Ac[1].imag+Ac2[0].imag+Ac2[1].imag)/2)
      ydi.append( (factor/(((Ed+E)/h)-step))*(Ad[0].imag+Ad[1].imag)/2)
    
    axs[l-l_start,0].plot(x,ya,'+:')
    axs[l-l_start,0].plot(x,yb,'+:')
    axs[l-l_start,0].plot(x,yc,'+:')
    axs[l-l_start,0].plot(x,yd,'+:')
  
    axs[l-l_start,1].plot(x,yai,'+:',label="a")
    axs[l-l_start,1].plot(x,ybi,'+:',label="b")
    axs[l-l_start,1].plot(x,yci,'+:',label="c")
    axs[l-l_start,1].plot(x,ydi,'+:',label="d")
    #axs[l-4,0].axvline(x=nuk, color="gray", linestyle=":")
    #axs[l-4,1].axvline(x=nuk, color="gray", linestyle=":")
    axs[l-l_start,0].set(ylabel="l=%d"%l)
    axs[l-l_start,0].axvline(x=Eb/h, color="gray", linestyle=":")
    axs[l-l_start,1].axvline(x=Eb/h, color="gray", linestyle=":")
    axs[l-l_start,0].axvline(x=Ec/h, color="gray", linestyle=":")
    axs[l-l_start,1].axvline(x=Ec/h, color="gray", linestyle=":")
    axs[l-l_start,0].axvline(x=Ed/h, color="gray", linestyle=":")
    axs[l-l_start,1].axvline(x=Ed/h, color="gray", linestyle=":")

  axs[0, 1].set_title('Imag')
  axs[0, 0].set_title('Real')
  axs[0, 1].legend(loc='upper right')
  axs[l_end-2, 0].set(xlabel='Frequency [1/s]')
  axs[l_end-2, 1].set(xlabel='Frequency [1/s]')
  plt.show()

if __name__ == "__main__":
  #print_Js()
  test()