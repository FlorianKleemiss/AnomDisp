import math
import numpy as np
from constants_and_atomic_properties import *

def z_EE(E,E2):
  return (E+abs(E2))/abs(E2)

def z_kb(k,b):
  return pow(k,2)/pow(b,2) + 1

def z_nprime(n_prime, n_0):
  return n_0*n_0/pow(n_prime,2) + 1

def n_prime_from_z(z,n_0):
  return n_0/math.sqrt(z-1)

def z_nunu(nu_j, nu_2):
  return nu_j/nu_2

def sugiura_exps(z,n_0):
  return np.exp(-4*n_0/np.sqrt(z-1)*np.arctan(np.sqrt(z-1)))/(1-np.exp(-2*n_0*math.pi/np.sqrt(z-1)))

def b(n_0, l_0, Z):
  Z_eff = None
  if n_0 == 1:
    Z_eff = get_Zeff_1s(Z)
  elif n_0 == 2:
    if l_0 == 0:
      Z_eff = get_Zeff_2s(Z)
    elif l_0 == 1:
      Z_eff = get_Zeff_2p_1_2(Z)
    elif l_0 == 2:
      Z_eff = get_Zeff_2p_3_2(Z)
  elif n_0 == 3:
    if l_0 ==0:
      Z_eff = get_Zeff_3s(Z)
    elif l_0 == 1:
      Z_eff = get_Zeff_3p_1_2(Z)
    elif l_0 == 2:
      Z_eff = get_Zeff_3d_3_2(Z)
  return Z_eff/(n_0*a0)

def k(E):
  return 2*math.pi / h * math.sqrt(2*el_mass*E)

def n_prime(E, Z, n,l):
  return 2*b(n,l,Z)/k(E)

def N0_square(b_):
  return pow(b_,3)/math.pi

def N0(b_):
  return math.sqrt(N0_square(b_))

def product_n_prime_from_z(n_0,z,l):
  n_p = n_0/math.sqrt(z-1)
  fact = 1.0
  for nu in range(1,l+1):
    #print(nu)
    fact *= n_p * n_p + nu*nu
  denom = 1-math.exp(-2*math.pi*n_p)
  return fact/denom
  
def N_square_from_z(l, m, b_, n_0, z):
  if (m > l):
    return 0
  result = 2*(2*l+1)*math.factorial(l-m)/math.factorial(l+m) * 2 * n_0 * math.pi * el_mass/h/h * b_ * product_n_prime_from_z(n_0,z,l)
  if m >= 1:
    result *= 2
  return result

def N(l, m, b_, n_0, z):
  return math.sqrt(N_square_from_z(l,m,b_,n_0,z))

def N_square(l, m, b_, n_0, z):
  return N_square_from_z(l,m,b_,n_0,z)

def N_lm_from_z(l,m,z,b_,n_0):
  return N(l,m,b_,n_0, z)

def N_lm_square_from_z(l,m,z,b_,n_0):
  return N_square(l,m,b_,n_0,z)
  
  
def q(nu):
  return 2*math.pi*nu/speed_of_light

def delta(a,b):
  if a == b:
      return 1
  else:
      return 0

######################### BEGIN OF MAKING K_p,l #######################################################

def prepare_M(p,l,z,n_0):
  nprime = n_0/math.sqrt(z-1)
  M = np.zeros((p+2,p+2),dtype=complex)
  M.fill(complex(0.0,0.0))
  M[0,0] = complex(1.0,0)
  M[1,0] = complex(0.0,nprime)
  M[1,1] = -2*(l+1)
  for j in range(1,p+1):
    M[j+1,0] = -0.25*M[j,1] + complex(0,nprime)*M[j,0]
    for s in range(1,j):
      M[j+1,s] = -0.25*(s+1)*M[j,s+1] + complex(0,nprime)*M[j,s] + (s-1-2*(l+j+1))*M[j,s-1]
    M[j+1,j] = complex(0,nprime) * M[j,j] + (j-1-2*(l+j+1))*M[j,j-1]
    M[j+1,j+1] = (j-2*(l+j+1))*M[j,j]
  return M

def g_from_M(M, xhi, j):
  sum = M[j,j]
  for i in range(j-1,-1,-1):
    sum *= xhi
    sum += M[j,i]
  return sum
  
def K_recursive_from_z(p, l, b_, z, n_0):
  if l > p+1:
    return 0
  xhi = complex(0,1/(2*math.sqrt(z-1)))
  M = prepare_M(p,l,z, n_0)
  banane = g_from_M(M, xhi, p+1-l)
  zm1 = z-1
  part1 = -pow(2,p+3+l) \
    * complex(0,math.pi) \
    * pow(math.sqrt(zm1),p+2+l) \
    / pow(-z,p+2) \
    / pow(complex(0,b_),p+2-l)
  ex = exp_from_z(z,n_0)
  return part1 * banane * ex

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
  elif c == 3:
    Q[c,0] = 15
    Q[c,2] = -15
  elif c == 4:
    Q[c,0] = 105
    Q[c,2] = -210
    Q[c,4] = 105
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
  M = a+c+l
  if Matrix is None:
    Matrix = make_matrix_W(a,c, M)
  result = value_from_W(Matrix,b,l,M)
  return result


W20 = make_matrix_W(2, 0, 20)
def J0_hat(p,l):
  if p == 0:
    return 1/3 * delta(l,0) - 1/15 * delta(l,2)
  #elif p == 1:
  #  return 1/15 * delta(l,1) - 1/35 * delta(l,3)
  #elif p == 2:
  #  return 1/15 * delta(l,0) + 1/105 * delta(l,2) - 4/315 * delta(l,4)
  #elif p == 3:
  #  return 1/35 * delta(l,1) - 1/315 * delta(l,3) - 4/693 * delta(l,5)
  else:
    return 0.25 * J(2,p,0,l, W20)
 
W00 = make_matrix_W(0, 0, 20)    
def J0_bar(p,l):
  if p == 0:
    return 2 * delta(l,0)
  #elif p == 1:
  #  return 2/3 * delta(l,1)
  #elif p == 2:
  #  return 2/3 * delta(l,0) + 4/15 * delta(l,2)
  #elif p == 3:
  #  return 2/5 * delta(l,1) + 4/35 * delta(l,3)
  #elif p == 4:
  #  return 2/5 * delta(l,0) + 8/35 * delta(l,2) + 16/315 * delta(l,4)
  #elif p == 5:
  #  return 2/7 * delta(l,1) + 8/63 * delta(l,3) + 16/693 * delta(l,5)
  else:
    return J(0,p,0,l,W00)

W11 = make_matrix_W(1, 1, 20)
def J1(p,l):
  if p == 0:
    return 4/3 * delta(l,1)
  #elif p == 1:
  #  return 4/5 * delta(l,2)
  #elif p == 2:
  #  return 4/15 * delta(l,1) + 16/35 * delta(l,3)
  #elif p == 3:
  #  return 12/35 * delta(l,2) + 16/63 * delta(l,4)
  #elif p == 4:
  #  return 4/35 * delta(l,1) + 32/105 * delta(l,3) + 32/231 * delta(l,5)
  else:
    return J(1,p,1,l,W11)

W22 = make_matrix_W(2, 2, 20)  
def J2(p,l):
  if p == 0:
    return 16/5 * delta(l,2)
  #elif p == 1:
  #  return 16/7 * delta(l,3)
  #elif p == 2:
  #  return 16/35 * delta(l,2) + 32/21 * delta(l,4)
  #elif p == 3:
  #  return 16/21 * delta(l,3) + 32/33 * delta(l,5)
  else:
    return J(2,p,2,l,W22)
W31 = make_matrix_W(3, 1, 20)
W33 = make_matrix_W(3, 3, 20)
################################# END OF CALC Js ###################################    
if True:
  def C_l_from_z(b_, z, n_0, l, nu, p_limit):
    k_ = b_*math.sqrt(z-1)
    part1 = b_/pow(-2*k_,l+1)
    sum = 0
    for p in range(p_limit):
      n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
      J_ = J(1,p,1,l,W11)
      if (J_ == 0):
        continue
      K1 = K_recursive_from_z(p,l,b_,z, n_0)
      K2 = K_recursive_from_z(p+1,l,b_,z, n_0)
      K2_mod = b_/2*K2
      sum += n1 * J_ * (K1-K2_mod)
    return part1 * sum
  
  def A_l_from_z(b_, z, n_0, l, nu, p_limit):
    k_ = b_*math.sqrt(z-1)
    part1 = b_/2/pow(-2*k_,l+1)
    sum = 0
    for p in range(p_limit):
      n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
      J = J(1,p,1,l,W11)
      if (J == 0):
        continue
      K1 = K_recursive_from_z(p,l,b_,z, n_0)
      sum += n1 * J * K1
    return part1 * sum
    
  def B2_from_z(b_, z, n_0, l, nu, p_limit):
    k_ = b_*math.sqrt(z-1)
    part1 = b_*b_/(4*pow(-2*k_,l+1))
    sum = 0
    for p in range(p_limit):
      n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
      J = J2(p,l)
      if (J == 0):
        continue
      K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
      sum += n1 * J * K1
    return part1 * sum
  
  def B1_from_z(b_, z, n_0, l, nu, p_limit):
    k_ = b_*math.sqrt(z-1)
    part1 = b_*b_/(2*pow(-2*k_,l+1))
    sum = 0
    for p in range(p_limit):
      n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
      J = J1(p+1,l)
      if J == 0:
        continue
      K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
      sum += n1 * J * K1
    return part1 * sum
    
  def B0_from_z(b_, z, n_0, l, nu, p_limit):
    k_ = b_*math.sqrt(z-1)
    part1 = b_/pow(-2*k_,l+1)
    sum = 0
    for p in range(p_limit):
      n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
      J1 = J0_hat(p,l)
      J2 = J0_bar(p,l)
      if J1 == 0 and J2 == 0:
        continue
      K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
      K2 = K_recursive_from_z(p,l,b_,z,n_0)
      sum += n1 * (2*b_*J1 * K1 - J2 * K2)
    return part1 * sum

def A_l_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -b_/2/pow(-2*k_,l+1)
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J_ = J(1,p,1,l,W11)
  if (J_ == 0):
    return 0.0
  K1 = K_recursive_from_z(p,l,b_,z, n_0)
  return part1 * n1 * J_ * K1

def C_l_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -b_/pow(-2*k_,l+1)
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J_ = J(1,p,1,l,W11)
  if (J_ == 0):
    return 0.0
  K1 = K_recursive_from_z(p, l, b_, z, n_0)
  K2 = K_recursive_from_z(p+1, l, b_, z, n_0)
  K2_mod = b_/2*K2
  return part1 * n1 * J_ * (K1-K2_mod)

def E_l_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -b_/2/pow(-2*k_,l+1)
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J_ = J(1,p,1,l,W11)
  if (J_ == 0):
    return 0.0
  K1 = K_recursive_from_z(p, l, b_, z, n_0)
  K2 = K_recursive_from_z(p+1, l, b_, z, n_0)
  K3 = K_recursive_from_z(p+2, l, b_, z, n_0)
  return part1 * n1 * J_ * (3*K1-10*b_/3*K2+2*b_*b_/3*K3)

def B2_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -b_*b_/(4*pow(-2*k_,l+1))
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J_ = J(2,p,2,l,W22)
  if (J_ == 0):
    return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  return part1 * n1 * J_ * K1

def B1_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -b_*b_/(2*pow(-2*k_,l+1))
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J_ = J(1,p+1,1,l,W11)
  if J_ == 0:
    return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  return part1 * n1 * J_ * K1
  
def B0_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -b_/(2*pow(-2*k_,l+1))
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J1 = J(2,p,0,l,W20)
  J2 = J(0,p,0,l,W00)
  if J1 == 0 and J2 == 0:
    return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p,l,b_,z,n_0)
  return part1 * n1 * (-2 * J2 * K2 + b_* J1 * K1)

def D2_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -math.sqrt(2./3.)*b_*b_/(4*pow(-2*k_,l+1))
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J_ = J(2,p,2,l,W22)
  if (J_ == 0):
    return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+2,l,b_,z,n_0)
  return part1 * n1 * J_ * (3*K1 - b_ * K2)

def D1_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -math.sqrt(2./3.)*b_*b_/(2*pow(-2*k_,l+1))
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J_ = J(1,p+1,1,l,W11)
  if J_ == 0:
    return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+2,l,b_,z,n_0)
  return part1 * n1 * J_ * (3 * K1 - b_ * K2)
  
def D0_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -math.sqrt(2./3.)*b_/(2*pow(-2*k_,l+1))
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J1 = J(2,p,0,l,W20) 
  J2 = J(0,p,0,l,W00) 
  if J1 == 0 and J2 == 0:
    return 0.0
  K1 = K_recursive_from_z(p,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K3 = K_recursive_from_z(p+2,l,b_,z,n_0)
  return part1 * n1 * (J2*(-4*K1 + 2*b_*K2) + b_ * J1 * (3*K2-b_*K3))

def G0_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -math.sqrt(2./3.)*b_*b_/2/pow(-2*k_,l+1)
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J1 = J(2,p+1,0,l,W20) 
  J2 = J(0,p+1,0,l,W00)
  if J1 == 0 and J2 == 0:
    return 0.0
  if J1 > 0:
    K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  else:
    K1 = 0
  if J2 > 0:
    K2 = K_recursive_from_z(p+2,l,b_,z,n_0)
  else:
    K2 = 0
  return part1 * n1 * (-2*J2*K1 + b_*J1*K2)

def G1_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -math.sqrt(2./3.)*b_*b_/2/pow(-2*k_,l+1)
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J1 = J(1,p,1,l,W11)
  J2 = J(3,p,1,l,W31)
  if J1 == 0 and J2 == 0:
    return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+2,l,b_,z,n_0)
  return part1 * n1 * (-J1*K1 + b_/4.*J2*K2)

def G2_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -math.sqrt(2./3.)*b_*b_*b_/4/pow(-2*k_,l+1)
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J1 = J(2,p+1,2,l,W22)
  if J1 == 0:
    return 0.0
  K1 = K_recursive_from_z(p+2,l,b_,z,n_0)
  return part1 * n1 * (J1*K1)

def G3_from_z_for_p(b_, z, n_0, l, nu, p):
  k_ = b_*math.sqrt(z-1)
  part1 = -math.sqrt(2./3.)*b_*b_*b_/8/pow(-2*k_,l+1)
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J1 = J(3,p,3,l,W33)
  if J1 == 0:
    return 0.0
  K1 = K_recursive_from_z(p+2,l,b_,z,n_0)
  return part1 * n1 * (J1*K1)

def G4_from_z_for_p(b_, z, n_0, l, nu, p): #'This is G_tilde'
  k_ = b_*math.sqrt(z-1)
  part1 = -math.sqrt(1./2.)*b_*b_/2/pow(-2*k_,l+1)
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  J1 = J(1,p,1,l,W11)
  J2 = J(1,p+2,1,l,W11)
  if J1 == 0 and J2 == 0.0:
    return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+2,l,b_,z,n_0)
  return part1 * n1 * (1./3.*J1*(2*K1-b_*K2) + b_*J2*K2)

################## END of matrix element calculation

## Start of f functions for angle independant part of matrix products:

def f_a(Z,l,k,z,nu_in,n_0,p_limit):
  if z <= 1: return 0
  b_ = b(n_0,0,Z)
  prefactor = pow(N0(b_),2) * N_lm_square_from_z(l,1,z,b_,n_0)
  matrix_value = 0
  if n_0 == 1:
    matrix_value = A_l_from_z(b_,z,n_0,l,nu_in,p_limit)
  elif n_0 == 2:
    matrix_value = C_l_from_z(b_,z,n_0,l,nu_in,p_limit)
  postfactor = matrix_value * matrix_value.conjugate()
  return prefactor*postfactor

def f_b(Z,l,g_k,z,nu_in,n_0,p_limit):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = pow(N0(b_),2) * N_lm_square_from_z(l,1,z,b_,n_0)
  matrix_value1 = B1_from_z(b_, z, n_0, l, nu_in, p_limit)
  if g_k == 0: 
    conjugate_function = B0_from_z
  elif g_k == 1:
    conjugate_function = B1_from_z
  elif g_k == 2:
    conjugate_function = B2_from_z 
  matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p_limit).conjugate()
  return prefactor * matrix_value1 * matrix_value2

def f_c_0(Z,l,g_k,z,nu_in,n_0,p_limit):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  N0sq = pow(N0(b_),2)
  prefactor = N0sq * N_square(l,0,b_,n_0,z)
  matrix_value1 = B0_from_z(b_, z, n_0, l, nu_in, p_limit)
  if g_k == 0: 
    conjugate_function = B0_from_z
  elif g_k == 1:
    conjugate_function = B1_from_z
  elif g_k == 2:
    conjugate_function = B2_from_z 
  matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p_limit).conjugate()
  return prefactor * matrix_value1 * matrix_value2

def f_c_2(Z,l,g_k,z,nu_in,n_0,p_limit):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  N0sq = pow(N0(b_),2)
  prefactor = N0sq * N_square(l,2,b_,n_0,z)
  matrix_value1 = B2_from_z(b_, z, n_0, l, nu_in, p_limit)
  if g_k == 0: 
    conjugate_function = B0_from_z
  elif g_k == 1:
    conjugate_function = B1_from_z
  elif g_k == 2:
    conjugate_function = B2_from_z 
  matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p_limit).conjugate()
  return prefactor * matrix_value1 * matrix_value2

def f_d(Z,l,g_k,z,nu_in,n_0,p_limit):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = pow(N0(b_),2) * N_square(l,2,b_,n_0,z)
  matrix_value1 = B2_from_z(b_, z, n_0, l, nu_in, p_limit)
  if g_k == 0: 
    conjugate_function = B0_from_z
  elif g_k == 1:
    conjugate_function = B1_from_z
  elif g_k == 2:
    conjugate_function = B2_from_z 
  matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p_limit).conjugate()
  return prefactor * matrix_value1 * matrix_value2

def integrand_matrix(z,f_function,z0, Z, l, k, nu_in, n_0, p_limit):
  if z<1: return 0
  return 2*z/(z*z-z0*z0)*f_function(Z,l,k,z,nu_in,n_0,p_limit).real

def integrand_matrix_s(z,z0, Z, l,k, nu_in, n_0, p_limit):
  if z<1: return 0
  return 2*z/(z*z-z0*z0) * \
      f_a(Z,l,k,z,nu_in,n_0,p_limit).real 

def integrand_matrix_p(z,z0, Z, nu_in, n_0, p_limit):
  if z<1: return 0
  return 2*z/(z*z-z0*z0)\
    *1/3*(\
      f_c_0(Z,0,0,z,nu_in,n_0,p_limit).real - \
      f_c_2(Z,2,0,z,nu_in,n_0,p_limit).real * 20\
        )

def f_a_for_p(Z,l,k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0,0,Z)
  prefactor = N0_square(b_) * N_square(l,1,b_,n_0,z)
  result = []
  if n_0 == 1:
    func = A_l_from_z_for_p
  elif n_0 == 2:
    func = C_l_from_z_for_p
  elif n_0 == 3:
    func = E_l_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = func(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_b_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,1,b_,n_0,z)
  result = []
  if n_0 == 2:
    func = B1_from_z_for_p
    if g_k == 0: 
      conjugate_function = B0_from_z_for_p
    elif g_k == 1:
      conjugate_function = B1_from_z_for_p
    elif g_k == 2:
      conjugate_function = B2_from_z_for_p
  elif n_0 == 3:
    func = D1_from_z_for_p
    if g_k == 0: 
      conjugate_function = D0_from_z_for_p
    elif g_k == 1:
      conjugate_function = D1_from_z_for_p
    elif g_k == 2:
      conjugate_function = D2_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_c_0_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,0,b_,n_0,z)
  result = []
  if n_0 == 2:
    func = B0_from_z_for_p
    if g_k == 0: 
      conjugate_function = B0_from_z_for_p
    elif g_k == 1:
      conjugate_function = B1_from_z_for_p
    elif g_k == 2:
      conjugate_function = B2_from_z_for_p
  elif n_0 == 3:
    func = D0_from_z_for_p
    if g_k == 0: 
      conjugate_function = D0_from_z_for_p
    elif g_k == 1:
      conjugate_function = D1_from_z_for_p
    elif g_k == 2:
      conjugate_function = D2_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_c_2_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,2,b_,n_0,z)
  result = []
  if n_0 == 2:
    func = B2_from_z_for_p
    if g_k == 0: 
      conjugate_function = B0_from_z_for_p
    elif g_k == 1:
      conjugate_function = B1_from_z_for_p
    elif g_k == 2:
      conjugate_function = B2_from_z_for_p
  elif n_0 == 3:
    func = D2_from_z_for_p
    if g_k == 0: 
      conjugate_function = D0_from_z_for_p
    elif g_k == 1:
      conjugate_function = D1_from_z_for_p
    elif g_k == 2:
      conjugate_function = D2_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_d_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,2,b_,n_0,z)
  result = []
  if n_0 == 2:
    func = B2_from_z_for_p
    if g_k == 0: 
      conjugate_function = B0_from_z_for_p
    elif g_k == 1:
      conjugate_function = B1_from_z_for_p
    elif g_k == 2:
      conjugate_function = B2_from_z_for_p
  elif n_0 == 3:
    func = D2_from_z_for_p
    if g_k == 0: 
      conjugate_function = D0_from_z_for_p
    elif g_k == 1:
      conjugate_function = D1_from_z_for_p
    elif g_k == 2:
      conjugate_function = D2_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_e_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,2,b_,n_0,z)
  result = []
  if n_0 == 3:
    func = G2_from_z_for_p
    if g_k == 0: 
      conjugate_function = G0_from_z_for_p
    elif g_k == 1:
      conjugate_function = G1_from_z_for_p
    elif g_k == 2:
      conjugate_function = G2_from_z_for_p
    elif g_k == 3:
      conjugate_function = G3_from_z_for_p
    elif g_k == 4:
      conjugate_function = G4_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_f_0_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,0,b_,n_0,z)
  result = []
  if n_0 == 3:
    func = G0_from_z_for_p
    if g_k == 0: 
      conjugate_function = G0_from_z_for_p
    elif g_k == 1:
      conjugate_function = G1_from_z_for_p
    elif g_k == 2:
      conjugate_function = G2_from_z_for_p
    elif g_k == 3:
      conjugate_function = G3_from_z_for_p
    elif g_k == 4:
      conjugate_function = G4_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_f_2_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,2,b_,n_0,z)
  result = []
  if n_0 == 3:
    func = G2_from_z_for_p
    if g_k == 0: 
      conjugate_function = G0_from_z_for_p
    elif g_k == 1:
      conjugate_function = G1_from_z_for_p
    elif g_k == 2:
      conjugate_function = G2_from_z_for_p
    elif g_k == 3:
      conjugate_function = G3_from_z_for_p
    elif g_k == 4:
      conjugate_function = G4_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_g_1_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,1,b_,n_0,z)
  result = []
  if n_0 == 3:
    func = G1_from_z_for_p
    if g_k == 0: 
      conjugate_function = G0_from_z_for_p
    elif g_k == 1:
      conjugate_function = G1_from_z_for_p
    elif g_k == 2:
      conjugate_function = G2_from_z_for_p
    elif g_k == 3:
      conjugate_function = G3_from_z_for_p
    elif g_k == 4:
      conjugate_function = G4_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_g_3_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,3,b_,n_0,z)
  result = []
  if n_0 == 3:
    func = G3_from_z_for_p
    if g_k == 0: 
      conjugate_function = G0_from_z_for_p
    elif g_k == 1:
      conjugate_function = G1_from_z_for_p
    elif g_k == 2:
      conjugate_function = G2_from_z_for_p
    elif g_k == 3:
      conjugate_function = G3_from_z_for_p
    elif g_k == 4:
      conjugate_function = G4_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_h_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,1,b_,n_0,z)
  result = []
  if n_0 == 3:
    func = G4_from_z_for_p
    if g_k == 0: 
      conjugate_function = G0_from_z_for_p
    elif g_k == 1:
      conjugate_function = G1_from_z_for_p
    elif g_k == 2:
      conjugate_function = G2_from_z_for_p
    elif g_k == 3:
      conjugate_function = G3_from_z_for_p
    elif g_k == 4:
      conjugate_function = G4_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_i_1_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,1,b_,n_0,z)
  result = []
  if n_0 == 3:
    func = G1_from_z_for_p
    if g_k == 0: 
      conjugate_function = G0_from_z_for_p
    elif g_k == 1:
      conjugate_function = G1_from_z_for_p
    elif g_k == 2:
      conjugate_function = G2_from_z_for_p
    elif g_k == 3:
      conjugate_function = G3_from_z_for_p
    elif g_k == 4:
      conjugate_function = G4_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_i_3_for_p(Z,l,g_k,z,nu_in,n_0,p):
  if z <= 1: return 0
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,3,b_,n_0,z)
  result = []
  if n_0 == 3:
    func = G3_from_z_for_p
    if g_k == 0: 
      conjugate_function = G0_from_z_for_p
    elif g_k == 1:
      conjugate_function = G1_from_z_for_p
    elif g_k == 2:
      conjugate_function = G2_from_z_for_p
    elif g_k == 3:
      conjugate_function = G3_from_z_for_p
    elif g_k == 4:
      conjugate_function = G4_from_z_for_p
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

## end of angle independant part

def print_Js():
  p_limit = 6
  l_limit = 6
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