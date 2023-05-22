import math
import numpy as np
from constants_and_atomic_properties import *

def z_EE(E: float,E2: float) -> float:
  return (E+abs(E2))/abs(E2)

def z_kb(k: float,b: float) -> float:
  return pow(k,2)/pow(b,2) + 1

def z_nprime(n_prime:float, n_0:int) -> float:
  return n_0*n_0/pow(n_prime,2) + 1

def n_prime_from_z(z:float,n_0:int) -> float:
  if z>=1:
    return n_0/math.sqrt(z-1)
  else:
    return n_0/math.sqrt((1-z))

def z_nunu(nu_j: float, nu_2: float) -> float:
  return nu_j/nu_2

def b(n_0: int, l_0:int , Z: int) -> float:
  Z_eff = -20
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
  if Z_eff == -20:
    raise ValueError("Errr in calculation of b, unkown case")
  return Z_eff/(n_0*a0)

def k(E: float) -> float:
  return 2*math.pi / h * math.sqrt(2*el_mass*E)

def n_prime(E: float, Z: int, n: int,l: int) -> float:
  return 2*b(n,l,Z)/k(E)

#introduces m^3 / pi
def N0_square(b_: float) -> float:
  return pow(b_,3)/math.pi

def N0(b_: float) -> float:
  return math.sqrt(N0_square(b_))

def product_n_prime_from_z(n_0: int,z: float,l: int) -> float:
  n_p = n_prime_from_z(z,n_0)
  fact = 1.0
  for nu in range(1,l+1):
    fact *= n_p * n_p + nu*nu
  denom = 1-math.exp(-2*math.pi*n_p)
  return fact/denom

#Introduces factors 2pi m_e /h^2 and m
def N_square_from_z(l: int, m: int, b_: float, n_0: int, z: float) -> float:
  if (m > l):
    return 0
  result = (2*l+1)*math.factorial(l-m)/math.factorial(l+m) * constant_factor / math.pi * n_0 * b_ * product_n_prime_from_z(n_0,z,l)
  if m >= 1:
    result *= 2
  return result

def N(l: int, m: int, b_: float, n_0: int, z: float) -> float:
  return math.sqrt(N_square_from_z(l,m,b_,n_0,z))

def N_square(l: int, m: int, b_: float, n_0: int, z: float) -> float:
  return N_square_from_z(l,m,b_,n_0,z)

def N_lm_from_z(l: int,m: int,z: float,b_: float,n_0: int) -> float:
  return N(l,m,b_,n_0, z)

def N_lm_square_from_z(l: int,m: int,z: float,b_: float,n_0: int) -> float:
  return N_square(l,m,b_,n_0,z)
  
#2pi m/s^2 
def q(nu: float) -> float:
  return 2*math.pi*nu/speed_of_light

def delta(a: Union[int,float],b: Union[int,float]) -> int:
  if a == b:
      return 1
  else:
      return 0

######################### BEGIN OF MAKING K_p,l #######################################################
matrix_type_M = type(np.zeros((0,0),dtype=complex))
def prepare_M(p:int,l:int,z:float,n_0:int) -> matrix_type_M:
  nprime = n_prime_from_z(z,n_0)
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

def g_from_M(M: matrix_type_M, xhi:complex, j:int) -> complex:
  sum = M[j,j]
  for i in range(j-1,-1,-1):
    sum *= xhi
    sum += M[j,i]
  return sum
  
def K_recursive_from_z(p: int, l: int, b_: float, z: float, n_0: int) -> complex:
  if l > p+1:
    return 0
  zm1 = z-1
  if z<= 1:
    zm1 = abs(zm1)
  xhi = complex(0,1/(2*math.sqrt(zm1)))
  banane = g_from_M(prepare_M(p,l,z, n_0), xhi, p+1-l)
  part1 = -pow(2,p+3+l) \
    * complex(0,math.pi) \
    * pow(math.sqrt(zm1),p+2+l) \
    / pow(-z,p+2) \
    / pow(complex(0,b_),p+2-l)
  ex = exp_from_z(z,n_0)
  if z<= 1:
    return -part1 * banane * ex
  else:
    return part1 * banane * ex
################### END OF CALC K ############################################
################### BEGIN CALC J  ############################################
matrix_type_W = type(np.zeros((0,0)))
def make_matrix_W(a: int,c: int, M: int) -> matrix_type_W:
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

def value_from_W(W: matrix_type_W,b: int,l: int,M: int) -> float:
  sum = 0
  for s in range(0,M+1):
    x = b+s+1
    if x%2 != 0:
      sum += W[l,s] * 2/x
  return sum

def J(a: int,b: int,c: int,l: int, Matrix: matrix_type_W = np.zeros((0,0))) -> float:
  M = a+c+l
  if Matrix is np.zeros((0,0)):
    Matrix = make_matrix_W(a,c, M)
  result = value_from_W(Matrix,b,l,M)
  return result


W20 = make_matrix_W(2, 0, 20)
#def J0_hat(p,l):
#  if p == 0:
#    return 1/3 * delta(l,0) - 1/15 * delta(l,2)
#  #elif p == 1:
#  #  return 1/15 * delta(l,1) - 1/35 * delta(l,3)
#  #elif p == 2:
#  #  return 1/15 * delta(l,0) + 1/105 * delta(l,2) - 4/315 * delta(l,4)
#  #elif p == 3:
#  #  return 1/35 * delta(l,1) - 1/315 * delta(l,3) - 4/693 * delta(l,5)
#  else:
#    return 0.25 * J(2,p,0,l, W20)
 
W00 = make_matrix_W(0, 0, 20)    
#def J0_bar(p,l):
#  if p == 0:
#    return 2 * delta(l,0)
#  #elif p == 1:
#  #  return 2/3 * delta(l,1)
#  #elif p == 2:
#  #  return 2/3 * delta(l,0) + 4/15 * delta(l,2)
#  #elif p == 3:
#  #  return 2/5 * delta(l,1) + 4/35 * delta(l,3)
#  #elif p == 4:
#  #  return 2/5 * delta(l,0) + 8/35 * delta(l,2) + 16/315 * delta(l,4)
#  #elif p == 5:
#  #  return 2/7 * delta(l,1) + 8/63 * delta(l,3) + 16/693 * delta(l,5)
#  else:
#    return J(0,p,0,l,W00)

W11 = make_matrix_W(1, 1, 20)
#def J1(p,l):
#  if p == 0:
#    return 4/3 * delta(l,1)
#  #elif p == 1:
#  #  return 4/5 * delta(l,2)
#  #elif p == 2:
#  #  return 4/15 * delta(l,1) + 16/35 * delta(l,3)
#  #elif p == 3:
#  #  return 12/35 * delta(l,2) + 16/63 * delta(l,4)
#  #elif p == 4:
#  #  return 4/35 * delta(l,1) + 32/105 * delta(l,3) + 32/231 * delta(l,5)
#  else:
#    return J(1,p,1,l,W11)

W22 = make_matrix_W(2, 2, 20)  
#def J2(p,l):
#  if p == 0:
#    return 16/5 * delta(l,2)
#  #elif p == 1:
#  #  return 16/7 * delta(l,3)
#  #elif p == 2:
#  #  return 16/35 * delta(l,2) + 32/21 * delta(l,4)
#  #elif p == 3:
#  #  return 16/21 * delta(l,3) + 32/33 * delta(l,5)
#  else:
#    return J(2,p,2,l,W22)
W31 = make_matrix_W(3, 1, 20)
W33 = make_matrix_W(3, 3, 20)
################################# END OF CALC Js ###################################    
def A_l_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J_ = J(1,p,1,l,W11)
  if (J_ == 0): return 0.0
  K1 = K_recursive_from_z(p,l,b_,z, n_0)
  if K1 == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  if (l+1)%2 == 0:
    part1 = -b_/2/pow(4*b_**2*(z-1),(l+1)/2)
  else:
    part1 = -b_/2/pow(-2*b_*math.sqrt(z-1),l+1)
  return part1 * n1 * J_ * K1

def C_l_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J_ = J(1,p,1,l,W11)
  if (J_ == 0): return 0.0
  K1 = K_recursive_from_z(p, l, b_, z, n_0)
  K2 = K_recursive_from_z(p+1, l, b_, z, n_0)
  K2_mod = b_/2*K2
  Ktot = (K1-K2_mod)
  if Ktot == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  if (l+1)%2 == 0:
    part1 = -b_/pow(4*b_**2*(z-1),(l+1)/2)
  else:
    part1 = -b_/pow(-2*b_*math.sqrt(z-1),l+1)
  return part1 * n1 * J_ * Ktot

def E_l_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J_ = J(1,p,1,l,W11)
  if (J_ == 0): return 0.0
  K1 = K_recursive_from_z(p, l, b_, z, n_0)
  K2 = K_recursive_from_z(p+1, l, b_, z, n_0)
  K3 = K_recursive_from_z(p+2, l, b_, z, n_0)
  Ktot = (3*K1-10*b_/3*K2+2*b_*b_/3*K3)
  if Ktot == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  if (l+1)%2 == 0:
    part1 = -b_/2/pow(4*b_**2*(z-1),(l+1)/2)
  else:
    part1 = -b_/2/pow(-2*b_*math.sqrt(z-1),l+1)
  return part1 * n1 * J_ * Ktot

def B2_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J_ = J(2,p,2,l,W22)
  if (J_ == 0): return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  if K1 == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  if (l+1)%2 == 0:
    part1 = -b_*b_/(4*pow(4*b_**2*(z-1),(l+1)/2))
  else:
    part1 = -b_*b_/(4*pow(-2*b_*math.sqrt(z-1),l+1))
  return part1 * n1 * J_ * K1

def B1_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J_ = J(1,p+1,1,l,W11)
  if J_ == 0: return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  if K1 == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  if (l+1)%2 == 0:
    part1 = -b_*b_/(2*pow(4*b_**2*(z-1),(l+1)/2))
  else:
    part1 = -b_*b_/(2*pow(-2*b_*math.sqrt(z-1),l+1))
  return part1 * n1 * J_ * K1
  
def B0_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J1 = J(2,p,0,l,W20)
  J2 = J(0,p,0,l,W00)
  if J1 == 0 and J2 == 0: return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p,l,b_,z,n_0)
  tot_term = (-2 * J2 * K2 + b_* J1 * K1)
  if tot_term == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  if (l+1)%2 == 0:
    part1 = -b_/(2*pow(4*b_**2*(z-1),(l+1)/2))
  else:
    part1 = -b_/(2*pow(-2*b_*math.sqrt(z-1),l+1))
  return part1 * n1 * tot_term

def D2_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J_ = J(2,p,2,l,W22)
  if (J_ == 0): return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+2,l,b_,z,n_0)
  Ktot = (3*K1 - b_ * K2)
  if Ktot == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  part1 = -math.sqrt(2./3.)*b_*b_/(4*pow(-2*b_*math.sqrt(z-1),l+1))
  return part1 * n1 * J_ * Ktot

def D1_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J_ = J(1,p+1,1,l,W11)
  if J_ == 0: return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+2,l,b_,z,n_0)
  Ktot = (3 * K1 - b_ * K2)
  if Ktot == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  part1 = -math.sqrt(2./3.)*b_*b_/(2*pow(-2*b_*math.sqrt(z-1),l+1))
  return part1 * n1 * J_ * Ktot
  
def D0_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J1 = J(2,p,0,l,W20) 
  J2 = J(0,p,0,l,W00) 
  if J1 == 0 and J2 == 0: return 0.0
  K1 = K_recursive_from_z(p,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K3 = K_recursive_from_z(p+2,l,b_,z,n_0)
  tot_term = (J2*(-4*K1 + 2*b_*K2) + b_ * J1 * (3*K2-b_*K3))
  if tot_term == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  part1 = -math.sqrt(2./3.)*b_/(2*pow(-2*b_*math.sqrt(z-1),l+1))
  return part1 * n1 * tot_term

def G0_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J1 = J(2,p+1,0,l,W20) 
  J2 = J(0,p+1,0,l,W00)
  if J1 == 0 and J2 == 0: return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+2,l,b_,z,n_0)
  tot_term = (-2*J2*K1 + b_*J1*K2)
  if tot_term == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  part1 = -math.sqrt(2./3.)*b_*b_/2./pow(-2*b_*math.sqrt(z-1),l+1)
  return part1 * n1 * tot_term

def G1_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J1 = J(1,p,1,l,W11)
  J2 = J(3,p,1,l,W31)
  if J1 == 0 and J2 == 0: return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+2,l,b_,z,n_0)
  tot_term = (-J1*K1 + b_/4.*J2*K2)
  if tot_term == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  part1 = -math.sqrt(2./3.)*b_*b_/2./pow(-2*b_*math.sqrt(z-1),l+1)
  return part1 * n1 * tot_term

def G2_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J1 = J(2,p+1,2,l,W22)
  if J1 == 0: return 0.0
  K1 = K_recursive_from_z(p+2,l,b_,z,n_0)
  if K1 == 0: return 0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  part1 = -math.sqrt(2./3.)*b_*b_*b_/4./pow(-2*b_*math.sqrt(z-1),l+1)
  return part1 * n1 * (J1*K1)

def G3_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex:
  J1 = J(3,p,3,l,W33)
  if J1 == 0: return 0.0
  K1 = K_recursive_from_z(p+2,l,b_,z,n_0)
  if K1 == 0: return 0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  part1 = -math.sqrt(2./3.)*b_*b_*b_/8./pow(-2*b_*math.sqrt(z-1),l+1)
  return part1 * n1 * (J1*K1)

def G4_from_z_for_p(b_: float, z: float, n_0: int, l: int, nu: float, p: int) -> complex: #'This is G_tilde'
  J1 = J(1,p,1,l,W11)
  J2 = J(1,p+2,1,l,W11)
  if J1 == 0 and J2 == 0.0: return 0.0
  K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
  K2 = K_recursive_from_z(p+2,l,b_,z,n_0)
  tot_term = (1./3.*J1*(2*K1-b_*K2) + b_*J2*K2)
  if tot_term == 0: return 0.0
  n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
  part1 = -math.sqrt(1./2.)*b_*b_/2./pow(-2*b_*math.sqrt(z-1),l+1)
  return part1 * n1 * tot_term

################## END of matrix element calculation

## Start of f functions for angle independant part of matrix products:
def f_a_for_p(Z: int,l: int,k: int,z: float,nu_in:float,n_0:int,p: int) -> "list[complex]":
  if z <= 0: return [complex(0,0)] * (2*p+1)
  b_ = b(n_0,0,Z)
  prefactor = N0_square(b_) * N_square(l,1,b_,n_0,z)
  result = []
  func = dummy_func
  if n_0 == 1:
    func = A_l_from_z_for_p
  elif n_0 == 2:
    func = C_l_from_z_for_p
  elif n_0 == 3:
    func = E_l_from_z_for_p
  else:
    raise ValueError("undefined shell")
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = func(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_p_el_for_p(Z: int,l: int,g_m: int,g_k: int,z: float,nu_in: float,n_0: int,p: int) -> "list[complex]":
  if z <= 0: return [complex(0,0)] * (p+1)
  b_ = b(n_0, 1, Z)
  prefactor = N0_square(b_) * N_square(l,g_m,b_,n_0,z)
  result = []
  func = dummy_func
  conjugate_function = dummy_func
  if n_0 == 2:
    if g_m == 0: 
      func = B0_from_z_for_p
    elif g_m == 1:
      func = B1_from_z_for_p
    elif g_m == 2:
      func = B2_from_z_for_p
    if g_k == 0: 
      conjugate_function = B0_from_z_for_p
    elif g_k == 1:
      conjugate_function = B1_from_z_for_p
    elif g_k == 2:
      conjugate_function = B2_from_z_for_p
  elif n_0 == 3:
    if g_m == 0: 
      func = D0_from_z_for_p
    elif g_m == 1:
      func = D1_from_z_for_p
    elif g_m == 2:
      func = D2_from_z_for_p
    if g_k == 0: 
      conjugate_function = D0_from_z_for_p
    elif g_k == 1:
      conjugate_function = D1_from_z_for_p
    elif g_k == 2:
      conjugate_function = D2_from_z_for_p
  else:
    raise ValueError("undefined shell")
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

def f_d_el_for_p(Z: int,l: int,g_m: int,g_k: int,z: float,nu_in: float,n_0: int,p: int) -> "list[complex]":
  if z <= 0: return [complex(0,0)] * (p+1)
  b_ = b(n_0, 2, Z)
  if g_m <= 3: prefactor = N0_square(b_) * N_square(l,g_m,b_,n_0,z)
  else: prefactor = N0_square(b_) * N_square(l,1,b_,n_0,z)
  result = []
  func = dummy_func
  conjugate_function = dummy_func
  if n_0 == 3:
    if g_m == 0: 
      func = G0_from_z_for_p
    elif g_m == 1:
      func = G1_from_z_for_p
    elif g_m == 2:
      func = G2_from_z_for_p
    elif g_m == 3:
      func = G3_from_z_for_p
    elif g_m == 4:
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
  else:
    raise ValueError("undefined shell")
  for j in range(p+1):
    matrix_value1 = func(b_,z,n_0,l,nu_in,j)
    matrix_value2 = conjugate_function(b_,z,n_0,l,nu_in,p-j)
    postfactor = matrix_value1 * matrix_value2.conjugate()
    result.append(prefactor * postfactor)
  return result

## end of angle independant part

def print_Js() -> None:
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