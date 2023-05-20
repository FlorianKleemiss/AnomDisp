from legendre_plynomials import alpha_bar_coef,alpha_coef,beta_bar_coef,beta_coef
from matrix_coefficients_v2 import f_p_el_for_p,f_a_for_p,f_d_el_for_p,sugiura_exps,q,b
#from matrix_coefficients_v2 import f_a,f_c_0, f_c_2, integrand_matrix_p, integrand_matrix_s
from constants_and_atomic_properties import *
import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np
import multiprocessing
from itertools import repeat
from hoenl_like import calc_hoenllike, sugiura_k_purely_imag
#from mpl_toolkits import mplot3d

al00  = []
al11  = []
al22  = []
al33  = []
bbl22 = []
bbl11 = []
bbl33 = []
for l in range(10):
  al00.append(alpha_coef(l,0,0,0,0))
  al11.append(alpha_coef(l,1,1,0,0))
  al22.append(alpha_coef(l,2,2,0,0))
  al33.append(alpha_coef(l,3,3,0,0))
  bbl22.append(beta_bar_coef(l,2,2,0,0))
  bbl11.append(beta_bar_coef(l,1,1,0,0))
  bbl33.append(beta_bar_coef(l,3,3,0,0))

def f_s_1_hoenl(z):
  part1 = 64/(3*pow(z,3))
  part2 = sugiura_exps(z,1)
  return part1*part2

def f_s_2_1_hoenl(z):
  part1 = 256*(z-2)/(15*pow(z,5))
  part2 = sugiura_exps(z,1)
  return part1*part2

def f_s_2_2_hoenl(z):
  part1 = 256*(4*z-3)/(15*pow(z,5))
  part2 = sugiura_exps(z,1)
  return part1*part2

def f_s_1_EM(z):
  part1 = 512*(z+3)/(3*pow(z,4))
  part2 = sugiura_exps(z,2)
  return part1*part2

def f_s_2_1_EM(z):
  part1 = 512*(z-2)*(z+3)/(15*pow(z,6))
  part2 = sugiura_exps(z,2)
  return part1*part2

def f_s_2_2_EM(z):
  part1 = 2048*pow(z-1,2)*(z+3)/(15*pow(z,7))
  part2 = sugiura_exps(z,2)
  return part1*part2

def f_p_1_EM(z):
  part1 = 512*(z+8./3.)/(3*pow(z,5))
  part2 = sugiura_exps(z,2)
  return part1*part2

def f_p_2_1_EM(z):
  return (z-2)/z/z/5 * f_p_1_EM(z)
  part1 = 512*(z-2)*(z+8./3.)/(15*pow(z,7))
  part2 = sugiura_exps(z,2)
  return part1*part2

def f_p_2_2_EM(z):
  return 2*(11*z-6)*(z+3)/(z+8./3.)/z/z/15 * f_p_1_EM(z)
  part1 = 1024*(11*z-6)*(z+3)/(45*pow(z,7))
  part2 = sugiura_exps(z,2)
  return part1*part2

def f_s_1_WA(z):
  part1 = 64*pow(3*z+4,2)*(z+8)/(pow(z,6))
  part2 = sugiura_exps(z,3)
  return part1*part2

def f_s_2_1_WA(z):
  part1 = 256*pow(3*z+4,2)*(z+8)*(z-2)/(45*pow(z,8))
  part2 = sugiura_exps(z,3)
  return part1*part2

def f_s_2_2_WA(z):
  part1 = 256*(3*z-4)*(z+8)*(4*z+5)*(3*z-4)/(45*pow(z,8))
  part2 = sugiura_exps(z,3)
  return part1*part2

def f_p_1_WA(z):
  part1 = 512*(3*z*z+26*z+28)/(pow(z,6))
  part2 = sugiura_exps(z,3)
  return part1*part2

def f_d_1_WA(z):
  part1 = 1024/5.0*(5*z*z+46*z+48)/(pow(z,7))
  part2 = sugiura_exps(z,3)
  return part1*part2

def xn(nu_in,n_0, el_nr, l_0, n):
  return pow(n_0 * q(nu_in)/b(n_0, l_0, el_nr), n)
def x2(nu_in,n_0, el_nr, l_0):
  return xn(nu_in,n_0, el_nr, l_0, 2)

cf2 = 4*math.pi**3*el_mass/el_charge**2

def apply_angle_part_s_parallel(a_l,theta0,alpha, l):
  a_l *= alpha_coef(l,1,1,theta0,alpha)
  return a_l

def apply_angle_part_s_orthogonal(a_l,theta0,alpha, l):
  a_l *= beta_coef(l,1,1,theta0,alpha)
  return a_l

def apply_angle_part_p_parallel(b_l,c_0_l,c_2_l,d_l,_k, theta0, alpha, l):
  ct0 = np.cos(theta0)
  st0 = np.sin(theta0)
  ca = np.cos(alpha)
  sa = np.sin(alpha)
  if _k == 0:
    b_l   *= -st0 * alpha_coef(l,1,0,theta0,alpha)
    c_0_l *=  ct0 * alpha_coef(l,0,0,theta0,alpha) * ca
    c_2_l *=  ct0 * alpha_coef(l,2,0,theta0,alpha) * ca
    d_l   *=  ct0 * alpha_bar_coef(l,2,0,theta0,alpha) * sa
  elif _k == 1:
    b_l   *= ct0 * alpha_coef(l,1,1,theta0,alpha)
    c_0_l *= st0 * alpha_coef(l,0,1,theta0,alpha) * ca
    c_2_l *= st0 * alpha_coef(l,2,1,theta0,alpha) * ca
    d_l   *= st0 * alpha_bar_coef(l,2,1,theta0,alpha) * sa
  elif _k == 2:
    b_l   *= -st0 * alpha_coef(l,1,2,theta0,alpha)
    c_0_l *=  ct0 * alpha_coef(l,0,2,theta0,alpha) * ca - sa * beta_coef(l,0,2,theta0,alpha)
    c_2_l *=  ct0 * alpha_coef(l,2,2,theta0,alpha) * ca - sa * beta_coef(l,2,2,theta0,alpha)
    d_l   *=  ct0 * alpha_bar_coef(l,2,2,theta0,alpha) * sa + ca * beta_bar_coef(l,2,2,theta0,alpha)
  return b_l,c_0_l,c_2_l,d_l

def apply_angle_part_p_orthogonal(b_l,c_0_l,c_2_l,d_l,_k, theta0, alpha, l):
  ct0 = np.cos(theta0)
  st0 = np.sin(theta0)
  ca = np.cos(alpha)
  sa = np.sin(alpha)
  if _k == 0:
    b_l   *= 0
    c_0_l *= -alpha_coef(l,0,0,theta0,alpha) * sa
    c_2_l *= -alpha_coef(l,2,0,theta0,alpha) * sa
    d_l   *= alpha_bar_coef(l,2,0,theta0,alpha) * ca
  elif _k == 1:
    b_l   *= ct0 * beta_coef(l,1,1,theta0,alpha)
    c_0_l *= st0 * beta_coef(l,0,1,theta0,alpha) * ca
    c_2_l *= st0 * beta_coef(l,2,1,theta0,alpha) * ca
    d_l   *= st0 * beta_bar_coef(l,2,1,theta0,alpha) * sa
  elif _k == 2:
    b_l   *= -st0 * beta_coef(l,1,2,theta0,alpha)
    c_0_l *=  ct0 * beta_coef(l,0,2,theta0,alpha) * ca + sa * alpha_coef(l,0,2,theta0,alpha)
    c_2_l *=  ct0 * beta_coef(l,2,2,theta0,alpha) * ca + sa * alpha_coef(l,2,2,theta0,alpha)
    d_l   *=  ct0 * beta_bar_coef(l,2,2,theta0,alpha) * sa - ca * alpha_bar_coef(l,2,2,theta0,alpha)
  return b_l,c_0_l,c_2_l,d_l

def apply_angle_part_d_parallel(vals,_k, theta0, alpha, l):
  ct0 = np.cos(theta0)
  c2t0 = np.cos(2*theta0)
  st0 = np.sin(theta0)
  s2t0 = np.sin(2*theta0)
  ca = np.cos(alpha)
  c2a = np.cos(2*alpha)
  sa = np.sin(alpha)
  s2a = np.sin(2*alpha)
  if _k == 0:
    vals[0] *=  c2t0 * alpha_coef(l,0,0,theta0,alpha) * ca
    vals[1] *=  c2t0 * alpha_coef(l,2,0,theta0,alpha) * ca
    vals[2] *=  c2t0 * sa * alpha_bar_coef(l,2,0,theta0,alpha)
    vals[3] *=  0.5*s2t0 * alpha_bar_coef(l,1,0,theta0,alpha) * s2a
    vals[4] *=  0.5*s2t0 * alpha_bar_coef(l,3,0,theta0,alpha) * s2a
    vals[5] *=  -np.sqrt(3.0)/2.0 * s2t0 * alpha_coef(l,1,0,theta0,alpha)
    vals[6] *=  0.5*s2t0*c2a *alpha_coef(l,1,0,theta0,alpha)
    vals[7] *=  0.5*s2t0*c2a *alpha_coef(l,3,0,theta0,alpha)
  elif _k == 1:
    vals[0] *=  st0*sa * beta_coef(l,0,1,theta0,alpha)      - 0.5*s2t0*ca*alpha_coef(l,0,1,theta0,alpha)
    vals[1] *=  st0*sa * beta_coef(l,2,1,theta0,alpha)      - 0.5*s2t0*ca*alpha_coef(l,2,1,theta0,alpha)
    vals[2] *= -st0*ca * beta_bar_coef(l,2,1,theta0,alpha)  - 0.5*s2t0*sa*alpha_bar_coef(l,2,1,theta0,alpha)
    vals[3] *=  ct0*c2a * beta_bar_coef(l,1,1,theta0,alpha) + (1-0.5*st0*st0)*s2a*alpha_bar_coef(l,1,1,theta0,alpha)
    vals[4] *=  ct0*c2a * beta_bar_coef(l,3,1,theta0,alpha) + (1-0.5*st0*st0)*s2a*alpha_bar_coef(l,3,1,theta0,alpha)
    vals[5] *=                                                np.sqrt(3.0)*0.5*st0*st0 * alpha_coef(l,1,1,theta0,alpha)
    vals[6] *= -ct0*s2a * beta_coef(l,1,1,theta0,alpha)     + (1-0.5*st0*st0)*c2a*alpha_coef(l,1,1,theta0,alpha)
    vals[7] *= -ct0*s2a * beta_coef(l,3,1,theta0,alpha)     + (1-0.5*st0*st0)*c2a*alpha_coef(l,3,1,theta0,alpha)
  elif _k == 2:
    vals[0] *=  c2t0*ca * alpha_coef(l,0,2,theta0,alpha)          - ct0*sa*beta_coef(l,0,2,theta0,alpha)
    vals[1] *=  c2t0*ca * alpha_coef(l,2,2,theta0,alpha)          - ct0*sa*beta_coef(l,2,2,theta0,alpha)
    vals[2] *=  c2t0*sa * alpha_bar_coef(l,2,2,theta0,alpha)      + ct0*ca*beta_bar_coef(l,2,2,theta0,alpha)
    vals[3] *=  0.5*s2t0*s2a * alpha_bar_coef(l,1,2,theta0,alpha) + st0*c2a*beta_bar_coef(l,1,2,theta0,alpha)
    vals[4] *=  0.5*s2t0*s2a * alpha_bar_coef(l,3,2,theta0,alpha) + st0*c2a*beta_bar_coef(l,3,2,theta0,alpha)
    vals[5] *=                                                    - np.sqrt(3.0)*0.5*s2t0 * alpha_coef(l,1,2,theta0,alpha)
    vals[6] *=  0.5*s2t0*c2a * alpha_coef(l,1,2,theta0,alpha)     - st0*s2a*beta_coef(l,1,2,theta0,alpha)
    vals[7] *=  0.5*s2t0*c2a * alpha_coef(l,3,2,theta0,alpha)     - st0*s2a*beta_coef(l,3,2,theta0,alpha)
  elif _k == 3:
    vals[0] *=  st0*sa * beta_coef(l,0,3,theta0,alpha)      - 0.5*s2t0*ca*alpha_coef(l,0,3,theta0,alpha)
    vals[1] *=  st0*sa * beta_coef(l,2,3,theta0,alpha)      - 0.5*s2t0*ca*alpha_coef(l,2,3,theta0,alpha)
    vals[2] *= -st0*ca * beta_bar_coef(l,2,3,theta0,alpha)  - 0.5*s2t0*sa*alpha_bar_coef(l,2,3,theta0,alpha)
    vals[3] *=  ct0*c2a * beta_bar_coef(l,1,3,theta0,alpha) + (1-0.5*st0*st0)*s2a*alpha_bar_coef(l,1,3,theta0,alpha)
    vals[4] *=  ct0*c2a * beta_bar_coef(l,3,3,theta0,alpha) + (1-0.5*st0*st0)*s2a*alpha_bar_coef(l,3,3,theta0,alpha)
    vals[5] *=                                                np.sqrt(3.0)*0.5*st0*st0 * alpha_coef(l,1,3,theta0,alpha)
    vals[6] *= -ct0*s2a * beta_coef(l,1,3,theta0,alpha)     + (1-0.5*st0*st0)*c2a*alpha_coef(l,1,3,theta0,alpha)
    vals[7] *= -ct0*s2a * beta_coef(l,3,3,theta0,alpha)     + (1-0.5*st0*st0)*c2a*alpha_coef(l,3,3,theta0,alpha)
  elif _k == 4:
    vals[0] *=  0.5*np.sqrt(3.0) *s2t0*ca * alpha_coef(l,0,1,theta0,alpha) 
    vals[1] *=  0.5*np.sqrt(3.0) *s2t0*ca * alpha_coef(l,2,1,theta0,alpha) 
    vals[2] *=  np.sqrt(3.0)*0.5*s2t0*sa * alpha_bar_coef(l,2,1,theta0,alpha)
    vals[3] *=  np.sqrt(3.0)*0.5*st0*st0*s2a * alpha_bar_coef(l,1,1,theta0,alpha)
    vals[4] *=  np.sqrt(3.0)*0.5*st0*st0*s2a * alpha_bar_coef(l,3,1,theta0,alpha)
    vals[5] *=  (1.5*ct0*ct0-0.5) * alpha_coef(l,1,1,theta0,alpha)
    vals[6] *=  np.sqrt(3.0)*0.5*st0*st0*c2a *alpha_coef(l,1,1,theta0,alpha)
    vals[7] *=  np.sqrt(3.0)*0.5*st0*st0*c2a *alpha_coef(l,3,1,theta0,alpha)
  return vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7]

def apply_angle_part_d_orthogonal(vals,_k, theta0, alpha, l):
  ct0 = np.cos(theta0)
  c2t0 = np.cos(2*theta0)
  st0 = np.sin(theta0)
  s2t0 = np.sin(2*theta0)
  ca = np.cos(alpha)
  c2a = np.cos(2*alpha)
  sa = np.sin(alpha)
  s2a = np.sin(2*alpha)
  if _k == 0:
    vals[0] *=  -ct0*sa * alpha_coef(l,0,0,theta0,alpha) 
    vals[1] *=  -ct0*sa * alpha_coef(l,2,0,theta0,alpha)
    vals[2] *=  ct0 * ca * alpha_bar_coef(l,2,0,theta0,alpha)  #
    vals[3] *=  st0 * c2a * alpha_bar_coef(l,1,0,theta0,alpha) #
    vals[4] *=  st0 * c2a * alpha_bar_coef(l,3,0,theta0,alpha) 
    vals[5] *=  0
    vals[6] *=  -st0*s2a * alpha_coef(l,1,0,theta0,alpha)
    vals[7] *=  -st0*s2a * alpha_coef(l,3,0,theta0,alpha)
  elif _k == 1:
    vals[0] *= 0.5*s2t0*ca * beta_coef(l,0,1,theta0,alpha)              + st0*sa*alpha_coef(l,0,1,theta0,alpha)      #1
    vals[1] *= 0.5*s2t0*ca * beta_coef(l,2,1,theta0,alpha)              + st0*sa*alpha_coef(l,2,1,theta0,alpha)      #1
    vals[2] *= 0.5*s2t0*sa * beta_bar_coef(l,2,1,theta0,alpha)          - st0*ca*alpha_bar_coef(l,2,1,theta0,alpha)  #2
    vals[3] *= -(1-0.5*st0*st0)*s2a * beta_bar_coef(l,1,1,theta0,alpha) + ct0*c2a*alpha_bar_coef(l,1,1,theta0,alpha) #2
    vals[4] *= -(1-0.5*st0*st0)*s2a * beta_bar_coef(l,3,1,theta0,alpha) + ct0*c2a*alpha_bar_coef(l,3,1,theta0,alpha) #2
    vals[5] *= -np.sqrt(3.0)*0.5*st0*st0 * beta_coef(l,1,1,theta0,alpha)                                             #
    vals[6] *= -(1-0.5*st0*st0)*c2a * beta_coef(l,1,1,theta0,alpha)     - ct0*s2a*alpha_coef(l,1,1,theta0,alpha)     #1
    vals[7] *= -(1-0.5*st0*st0)*c2a * beta_coef(l,3,1,theta0,alpha)     - ct0*s2a*alpha_coef(l,3,1,theta0,alpha)     #1
  elif _k == 2:
    vals[0] *= c2t0*ca * beta_coef(l,0,2,theta0,alpha)           + ct0*sa*alpha_coef(l,0,2,theta0,alpha)      #1
    vals[1] *= c2t0*ca * beta_coef(l,2,2,theta0,alpha)           + ct0*sa*alpha_coef(l,2,2,theta0,alpha)      #1
    vals[2] *= -ct0*ca * alpha_bar_coef(l,2,2,theta0,alpha)      + c2t0*sa*beta_bar_coef(l,2,2,theta0,alpha)  #1
    vals[3] *= -st0*c2a * alpha_bar_coef(l,1,2,theta0,alpha) + 0.5*s2t0*s2a*beta_bar_coef(l,1,2,theta0,alpha) #1
    vals[4] *= -st0*c2a * alpha_bar_coef(l,3,2,theta0,alpha) + 0.5*s2t0*s2a*beta_bar_coef(l,3,2,theta0,alpha) #1
    vals[5] *= -np.sqrt(3.0)*0.5*s2t0 * beta_coef(l,1,2,theta0,alpha)                                         #
    vals[6] *= 0.5*s2t0*c2a * beta_coef(l,1,2,theta0,alpha)     + st0*s2a*alpha_coef(l,1,2,theta0,alpha)      #1
    vals[7] *= 0.5*s2t0*c2a * beta_coef(l,3,2,theta0,alpha)     + st0*s2a*alpha_coef(l,3,2,theta0,alpha)      #1
  elif _k == 3:
    vals[0] *= -0.5*s2t0*ca * beta_coef(l,0,3,theta0,alpha)      - st0*sa*alpha_coef(l,0,3,theta0,alpha)
    vals[1] *= -0.5*s2t0*ca * beta_coef(l,2,3,theta0,alpha)      - st0*sa*alpha_coef(l,2,3,theta0,alpha)
    vals[2] *= -0.5*s2t0*sa * beta_bar_coef(l,2,3,theta0,alpha)  + st0*ca*alpha_bar_coef(l,2,3,theta0,alpha)
    vals[3] *= (1-0.5*st0*st0)*s2a * beta_bar_coef(l,1,3,theta0,alpha) - ct0*c2a*alpha_bar_coef(l,1,3,theta0,alpha)
    vals[4] *= (1-0.5*st0*st0)*s2a * beta_bar_coef(l,3,3,theta0,alpha) - ct0*c2a*alpha_bar_coef(l,3,3,theta0,alpha)
    vals[5] *= np.sqrt(3.0)*0.5*st0*st0 * beta_coef(l,1,3,theta0,alpha)
    vals[6] *= (1-0.5*st0*st0)*c2a * beta_coef(l,1,3,theta0,alpha)     + ct0*s2a*alpha_coef(l,1,3,theta0,alpha)
    vals[7] *= (1-0.5*st0*st0)*c2a * beta_coef(l,3,3,theta0,alpha)     + ct0*s2a*alpha_coef(l,3,3,theta0,alpha)
  elif _k == 4:
    vals[0] *= 0.5*np.sqrt(3.0)*s2t0*ca * beta_coef(l,0,1,theta0,alpha) 
    vals[1] *= 0.5*np.sqrt(3.0)*s2t0*ca * beta_coef(l,2,1,theta0,alpha) 
    vals[2] *= 0.5*np.sqrt(3.0)*s2t0*sa * beta_bar_coef(l,2,1,theta0,alpha)
    vals[3] *= 0.5*np.sqrt(3.0)*st0*st0*s2a * beta_bar_coef(l,1,1,theta0,alpha)
    vals[4] *= 0.5*np.sqrt(3.0)*st0*st0*s2a * beta_bar_coef(l,3,1,theta0,alpha)
    vals[5] *= (1.5*ct0*ct0-0.5) * beta_coef(l,1,1,theta0,alpha)
    vals[6] *= np.sqrt(3.0)*0.5*st0*st0*c2a * beta_coef(l,1,1,theta0,alpha)
    vals[7] *= np.sqrt(3.0)*0.5*st0*st0*c2a * beta_coef(l,3,1,theta0,alpha)
  return vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7]

def calc_Intensity_s_orbital(alpha_loc, nu_in, t0, l_max, p_max, n0, el_nr):
  if n0 == 1:
    z_temp = nu_in / (get_ionization_energy_1s(el_nr) / h)
    de = delta(el_nr,1,0,0)
  elif n0 == 2:
    z_temp = nu_in / (get_ionization_energy_2s(el_nr) / h)
    de = delta(el_nr,2,0,0)
  elif n0 == 3:
    z_temp = nu_in / (get_ionization_energy_3s(el_nr) / h)
    de = delta(el_nr,3,0,0)
  if z_temp < 1: return 0
  z_temp *= (1-de)

  par = 0
  orth = 0
  fac = []
  for p in range(int(p_max/2+1)):
    fac.append(p+1)
  for p in reversed(range(int(p_max/2+1))):
    fac.append(p+1)

  for l in range(l_max+1):
    temp = 0
    for p in range(0,p_max+1,2):
      for r in f_a_for_p(el_nr, l, 0, z_temp, nu_in, n0, p):
        temp += r.real
    par += alpha_coef(l,1,1,t0,alpha_loc) * temp
    orth += beta_coef(l,1,1,t0,alpha_loc) * temp
  
  return (par**2 + orth**2)

def calc_Intensity_p_orbital(alpha_loc, nu_in, t0, l_max, p_max, n0, subshell, el_nr):
  if n0 == 1:
    print("This shell doesn't have p orbitals")
    exit(-1)
  elif n0 == 2:
    if subshell == 1:
      z_temp = nu_in / (get_ionization_energy_2p1_2(el_nr) / h)
      de = delta(el_nr,2,1,0.5)
    else:
      z_temp = nu_in / (get_ionization_energy_2p3_2(el_nr) / h)
      de = delta(el_nr,2,1,1.5)
  elif n0 == 3:
    if subshell == 1:
      z_temp = nu_in / (get_ionization_energy_3p_1_2(el_nr) / h)
      de = delta(el_nr,3,1,0.5)
    else:
      z_temp = nu_in / (get_ionization_energy_3p_3_2(el_nr) / h)
      de = delta(el_nr,3,1,1.5)
  if z_temp < 1: return 0
  z_temp *= (1-de)

  par = 0
  orth = 0
  fac = []
  for p in range(int(p_max/2+1)):
    fac.append(p+1)
  for p in reversed(range(int(p_max/2+1))):
    fac.append(p+1)

  ms = [0,1,2]
  for l in range(l_max+1):
    mults = [al00[l],al11[l],al22[l]+bbl22[l]]
    for p in range(0,p_max+1,2):
      for nummy in range(len(ms)):
        for num,r in enumerate(f_p_el_for_p(el_nr, l, ms[nummy], ms[nummy], z_temp, nu_in, n0, p)):
          
          par += apply_angle_part_s_parallel(r.real * mults[nummy], t0, alpha_loc, fac[num])
          orth += apply_angle_part_s_orthogonal(r.real * mults[nummy], t0, alpha_loc, fac[num])
  
  #Has division by 3**2 to avergae over all 3 p orbitals a.k.a. the p-wonder
  return (par**2 + orth**2) / 9

def calc_Intensity_d_orbital(alpha_loc, nu_in, t0, l_max, p_max, n0, subshell, el_nr):
  if n0 == 1:
    print("This shell doesn't have d orbitals")
    exit(-1)
  elif n0 == 2:
    print("This shell doesn't have d orbitals")
    exit(-1)
  elif n0 == 3:
    if subshell == 1:
      z_temp = nu_in / (get_ionization_energy_3d_3_2(el_nr) / h)
      de = delta(el_nr,3,2,1.5)
    else:
      z_temp = nu_in / (get_ionization_energy_3d_5_2(el_nr) / h)
      de = delta(el_nr,3,2,2.5)
  if z_temp < 1: return 0
  z_temp *= (1-de)

  par = 0
  orth = 0
  fac = []
  for p in range(int(p_max/2+1)):
    fac.append(p+1)
  for p in reversed(range(int(p_max/2+1))):
    fac.append(p+1)

  ms = [0,1,2,3,4]
  for l in range(l_max+1):
    mults = [al00[l],bbl11[l]+al11[l],al22[l]+bbl22[l],bbl33[l]+al33[l],al11[l]]
    for p in range(0,p_max+1,2):
      for nummy in range(len(ms)):
        for num,r in enumerate(f_d_el_for_p(el_nr, l, ms[nummy], ms[nummy], z_temp, nu_in, n0, p)):
          
          par += apply_angle_part_s_parallel(r.real * mults[nummy], t0, alpha_loc, fac[num])
          orth += apply_angle_part_s_orthogonal(r.real * mults[nummy], t0, alpha_loc, fac[num])
  #Has division by 5**2 to accomodate the averaging over 5 d orbtials  a.k.a. the d-wonder
  return (par**2 + orth**2) / 25.0

def calc_stuff(nu_in, t0, l_max, p_max, el_nr):
  #Integrate Intensity over all poalrization angles alpha and then divide by 2pi and constant_coef, to get rid of physical prefactor
  k_s = integrate.quad(calc_Intensity_s_orbital, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 1, el_nr))[0] / (2*math.pi*constant_factor**2)
  l_s = integrate.quad(calc_Intensity_s_orbital, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 2, el_nr))[0] / (2*math.pi*constant_factor**2)
  m_s = integrate.quad(calc_Intensity_s_orbital, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 3, el_nr))[0] / (2*math.pi*constant_factor**2)

  l_p1 = integrate.quad(calc_Intensity_p_orbital, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 2, 1, el_nr))[0] / (2*math.pi*constant_factor**2)
  l_p3 = integrate.quad(calc_Intensity_p_orbital, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 2, 2, el_nr))[0] / (2*math.pi*constant_factor**2)
  m_p1 = integrate.quad(calc_Intensity_p_orbital, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 3, 1, el_nr))[0] / (2*math.pi*constant_factor**2)
  m_p3 = integrate.quad(calc_Intensity_p_orbital, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 3, 2, el_nr))[0] / (2*math.pi*constant_factor**2)

  m_d3 = integrate.quad(calc_Intensity_d_orbital, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 3, 1, el_nr))[0] / (2*math.pi*constant_factor**2)
  m_d5 = integrate.quad(calc_Intensity_d_orbital, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 3, 2, el_nr))[0] / (2*math.pi*constant_factor**2)

  #returns the evaluation of f(z*(1-delta)) averaged over incident and outgoing polarization angle
  return math.sqrt(k_s), math.sqrt(l_s), math.sqrt(m_s), math.sqrt(l_p1), math.sqrt(l_p3), \
         math.sqrt(m_p1), math.sqrt(m_p3), math.sqrt(m_d3), math.sqrt(m_d5)

def calc_stuff_with_brennan(nu_in, t0, l_max, p_max, el_nr, br, e):
  wavelength = speed_of_light / nu_in * 1E10

  brennan_fdp = br.at_angstrom(wavelength, e)[1]
  hönl_K,hönl_L,hönl_M,hönl = calc_hoenllike(h*nu_in, el_nr)
  t1,t2,t3,t4,t5,t6,t7,t8,t9 = calc_stuff(nu_in, t0, l_max, p_max, el_nr)
  
  return t1,t2,t3,t4,t5,t6,t7,t8,t9, brennan_fdp, hönl_K, hönl_L, hönl_M, hönl

def calc_stuff_only_sums(nu_in, t0, l_max, p_max, el_nr):
  t1,t2,t3,t4,t5,t6,t7,t8,t9 = calc_stuff(nu_in, t0, l_max, p_max, el_nr)

  return t1,t2,t3,t4+2*t5,t6+2*t7,3*t8+2*t9

def fdp_from_z(z, delta, osz_func):
  if z<1:  return 0
  return osz_func((1 - delta) * z)*2*math.pi
  #return 2**7/3.0 *sugiura_exps((1 - delta) * z,n0)/((1 - delta) * z)**3*math.pi

def eta_K(z,n0,nu_0,delta_K):
  if z<1: return 0
  emd = 1 - delta_K
  z_eff = emd * z
  nu = z * nu_0
  #return 2**7/3/z_eff**4*sugiura_exps(z_eff,n0)/nu_0/nu*math.pi/2*emd
  return f_s_1_hoenl(z_eff)/nu**2*math.pi

def eta_K_unsers(z, nu_0, el_nr, delta_K):
  #emd = 1 - delta_K
  nu_in = z * nu_0
  #z_eff = emd * z
  f = calc_stuff_only_sums(nu_in,0,7,1,el_nr)[0]
  #return 2**7/3/z_eff**4*sugiura_exps(z_eff,n0)/nu_0/nu*math.pi/2*emd
  return f/nu_in**2*math.pi

if __name__ == "__main__":
  test_1s = False
  higher_p_test = False
  M_shell_test = False
  M_shell_test_orthogonal = False
  form_factor_test = False
  numplot=False
  cmaps = ['blue','orange','green','red','black','blue','orange','green','red','black']
  t0 = 0
  alp = 0
  l_max = 7
  p_max = 1
  el_nr = 35
  lam_min = 3.00
  lam_max = 0.15
  steps = 300
  x = None

  import tkinter as tk
  from tkinter import ttk
  root = tk.Tk()
  root.title("Calculation Setup Selection")
  root.resizable(True,False)
  v = tk.IntVar()
  v.set(0)
  v2 = tk.StringVar()
  v2.set(str(el_nr))
  v3 = tk.StringVar()
  v3.set(str(l_max))
  v4 = tk.StringVar()
  v4.set(str(p_max))
  v5 = tk.StringVar()
  v5.set(str(lam_max))
  v6 = tk.StringVar()
  v6.set(str(lam_min))
  v7 = tk.StringVar()
  v7.set(str(steps))
  options = [
        ('1s-test', 0),
        ('p test', 2),
        ('M miracle', 3),
        ('M miracle orthogonal', 5),
        ('Dispersion Coefficient Plot', 4),
        ('Numerical_plot', 6),
        ]
  global result
  result = 0
  def set_result():
    global result
    result = v.get()
    root.destroy()
  tk.Label(root,text="Choose Z of element",state="active",justify=tk.CENTER,padx=20).pack()
  tk.Entry(root, textvariable=v2).pack()
  tk.Label(root,text="Choose l_max",state="active",justify=tk.CENTER,padx=20).pack()
  tk.Entry(root, textvariable=v3).pack()
  tk.Label(root,text="Choose nr of steps",state="active",justify=tk.CENTER,padx=20).pack()
  tk.Entry(root, textvariable=v7).pack()
  tk.Label(root,text="Below only for disp calc",state="active",justify=tk.CENTER,padx=20).pack()
  ttk.Separator(root,orient="horizontal").pack(fill='x')
  tk.Label(root,text="Choose p_max",state="active",justify=tk.CENTER,padx=20).pack()
  tk.Entry(root, textvariable=v4).pack()
  tk.Label(root,text="Choose lambda max",state="active",justify=tk.CENTER,padx=20).pack()
  tk.Entry(root, textvariable=v5).pack()
  tk.Label(root,text="Choose lambda min",state="active",justify=tk.CENTER,padx=20).pack()
  tk.Entry(root, textvariable=v6).pack()
  tk.Label(root,text="Choose what to run",state="active",justify=tk.CENTER,padx=20).pack()
  for l,val in options:
    tk.Radiobutton(root,text=l,padx=20,indicatoron=0, width=30, height=2, variable=v,command=set_result,value=val).pack()
  root.mainloop()
  el_nr = int(v2.get())
  l_max = int(v3.get())
  p_max = int(v4.get())
  lam_max = float(v5.get())
  lam_min = float(v6.get())
  steps = int(v7.get())
  if result == 0:
    test_1s = True
  elif result == 2:
    higher_p_test = True
  elif result == 3:
    M_shell_test = True
  elif result == 4:
    form_factor_test = True
  elif result == 5:
    M_shell_test_orthogonal = True
  elif result == 6:
    numplot = True

  if test_1s == True:
    b_ = b(1,0,1)
    lambd = 0.7E-10
    blambd = b_ * lambd
    pi2_blambd = 2 *math.pi / blambd
    pi2_b = 2 *math.pi / b_
    lam = lambd*1E10 #this is in angstrom
    lam_bohr = lambd/a0
    def final(p):
      return np.vectorize(lambda n: 
        np.sum(np.fromiter((pow(-1,x) * pow(pi2_b * n * 1E10, 2*x) * (x+1) for x in range(p,-1,-1)), dtype=np.float64))
      )
    def sin_t_l_function(t,l):
      return np.round(np.sin(t/2)/l,4)
    def stl_to_t0(stl,l):
      return np.arcsin(stl*l)*2

    stl = np.array([0.01,0.02,0.03,0.04,0.05, 
           0.06,0.07,0.08,0.09,0.10, 
           0.11,0.12,0.13,0.14,0.15, 
           0.16,0.17,0.18,0.19,0.20, 
           0.22,0.24,0.25,0.26,0.28,0.30, 
           0.32,0.34,0.35,0.36,0.38,0.40, 
           0.42,0.44,0.45,0.46,0.48,0.50,
           0.55,0.60,0.65,0.70,0.80,0.90,1.00])

    scat = np.array([0.998,0.991,0.980,0.966,0.947, 
           0.925,0.900,0.872,0.842,0.811, 
           0.778,0.744,0.710,0.676,0.641, 
           0.608,0.574,0.542,0.511,0.481, 
           0.424,0.373,0.350,0.328,0.287,0.251, 
           0.220,0.193,0.180,0.169,0.148,0.130,
           0.115,0.101,0.095,0.090,0.079,0.071,
           0.053,0.040,0.031,0.024,0.015,0.010,0.007])

    stl2 = np.array([0.0   ,0.0215,0.0429,0.0644,0.0859,
                     0.1073,0.1288,0.1503,0.1718,0.1932,
                     0.2147,0.2576,0.3006,0.3435,0.3864,
                     0.4294,0.4723,0.5153,0.5582,0.6011,
                     0.6441,0.6870,0.7300,0.7729,0.8158,
                     0.8588,0.9017,0.9447,0.9876,1.0305])

    scat2 = np.array([1.0  ,0.9924,0.9704,0.9352,0.8892,
                     0.8350,0.7752,0.7125,0.6492,0.5871,
                     0.5277,0.4201,0.3301,0.2573,0.1998,
                     0.1552,0.1208,0.0945,0.0744,0.0592,
                     0.0474,0.0383,0.0311,0.0254,0.0208,
                     0.0171,0.0140,0.0116,0.0096,0.0080])

    x = np.linspace(0,1.5,steps)
    fig = plt.figure()
    axes = fig.add_subplot(111)
    axes.plot(stl,scat,"--",label="ITC")
    axes.plot(stl2,scat2,"--",label="ITC_bound")
    Z2 = ( 0.493002*np.exp(-10.5109 * pow(x, 2))
          +0.322912*np.exp(-26.1257 * pow(x, 2))
          +0.140191*np.exp(-3.14236 * pow(x, 2))
          +0.040810*np.exp(-57.7997 * pow(x, 2))
          +0.003038)
    axes.plot(x,Z2,".",label="SFAC")
    for p in range(50,301,150):
      print(p)
      Z = final(p)(x)
      axes.plot(x,Z, label="t = %d"%p)
    x_in_meter = x*1E10
    temp_k = 4*math.pi*a0*x_in_meter
    Z3 = 16/pow(4+pow(temp_k,2),2)
    axes.plot(x,Z3,"-.",label="Thakkar Model")
    temp_z = 1.15
    Z4 = 16*pow(temp_z,4)/pow(pow(2*temp_z,2)+pow(temp_k,2),2)
    axes.scatter(x,Z4,label="Thakkar Model Extended")
    axes.legend()
    axes.set_xlabel('sin (theta) / lambda [Ang^-1]')
    axes.set_ylabel('Scattering power')
    axes.set_ylim(0,1.1)
    axes.set_xlim(0,1.5)
    plt.show()
    exit()
  
  if higher_p_test == True:
    print("Performing Test for higher values of P")
    nu_in = 0.8*get_ionization_energy_1s(el_nr)/h
    x = np.linspace(1.0001,3.0001,steps)
    k_ = 0
    x1_2 = xn(nu_in, 1, el_nr, 0, 2)
    x20_2 = xn(nu_in, 2, el_nr, 0, 2)
    x21_2 = xn(nu_in, 2, el_nr, 1, 2)
    x30_2 = xn(nu_in, 3, el_nr, 0, 2)
    x31_2 = xn(nu_in, 3, el_nr, 1, 2)
    x32_2 = xn(nu_in, 3, el_nr, 2, 2)
    from time import time
    numpy.random.seed(int(np.ceil(time())))
    t0 = numpy.random.random(1)[0] * math.pi
    #alp = numpy.random.random(1)[0] * math.pi * 2
    third = 1./3.
    res_2 = [0,0,0]
    res_4 = [0,0,0,0,0]
    res_6 = [0,0,0,0,0,0,0]
  
    if True:
      a_result  = np.zeros_like(x)
      b_result  = np.zeros_like(x)
      c_result  = np.zeros_like(x)
      b0_result = np.zeros_like(x)
      b1_result = np.zeros_like(x)
      b2_result = np.zeros_like(x)
      k_s       = np.zeros_like(x)
      d_result  = np.zeros_like(x)
      e_result  = np.zeros_like(x)
      e0_result = np.zeros_like(x)
      e1_result = np.zeros_like(x)
      e2_result = np.zeros_like(x)
      l_s       = np.zeros_like(x)
      f_result  = np.zeros_like(x)
    
      p0_result   = np.zeros_like(x)
      p2_0_result = np.zeros_like(x)
      p2_1_result = np.zeros_like(x)
      p4_0_result = np.zeros_like(x)
      p4_1_result = np.zeros_like(x)
      p4_2_result = np.zeros_like(x)
      p6_0_result = np.zeros_like(x)
      p6_1_result = np.zeros_like(x)
      p6_2_result = np.zeros_like(x)
      p6_3_result = np.zeros_like(x)
      l_p         = np.zeros_like(x)
      
      p0_result_M   = np.zeros_like(x)
      p2_0_result_M = np.zeros_like(x)
      p2_1_result_M = np.zeros_like(x)
      p4_0_result_M = np.zeros_like(x)
      p4_1_result_M = np.zeros_like(x)
      p4_2_result_M = np.zeros_like(x)
      p6_0_result_M = np.zeros_like(x)
      p6_1_result_M = np.zeros_like(x)
      p6_2_result_M = np.zeros_like(x)
      p6_3_result_M = np.zeros_like(x)
      l_p_M         = np.zeros_like(x)
    
      g_result = np.zeros_like(x)
      h_result = np.zeros_like(x)
      i_result = np.zeros_like(x)
      j_result = np.zeros_like(x)
      k_result = np.zeros_like(x)
      l_result = np.zeros_like(x)
    
      m_result = np.zeros_like(x)
      n_result = np.zeros_like(x)
      o_result = np.zeros_like(x)
    for i,z in enumerate(x):
      for l in range(l_max+1):
        res_0 = f_a_for_p(el_nr, l, k_, z, nu_in, 1, 0)[0]
        res_0 = apply_angle_part_s_parallel(res_0.real, t0, alp, l)
        res_2 = f_a_for_p(el_nr, l, k_, z, nu_in, 1, 2)
        for runny in range(len(res_2)):
          res_2[runny] = apply_angle_part_s_parallel(res_2[runny].real,t0,alp, l)
        res_4 = f_a_for_p(el_nr,l,k_,z,nu_in,1,4)
        for runny in range(len(res_4)):
          res_4[runny] = apply_angle_part_s_parallel(res_4[runny].real,t0,alp, l)
        a_result[i] += (res_0)
        b_result[i] += ((res_2[0]+res_2[2]))
        c_result[i] += (res_2[1])
        b0_result[i] += ((res_4[0]+res_4[-1]))
        b1_result[i] += ((res_4[1]+res_4[-2]))
        b2_result[i] += ((res_4[2]))
        res_0 = f_a_for_p(el_nr, l, k_, z, nu_in, 2, 0)[0]
        res_0 = apply_angle_part_s_parallel(res_0.real, t0, alp, l)
        res_2 = f_a_for_p(el_nr, l, k_, z, nu_in, 2, 2)
        for runny in range(len(res_2)):
          res_2[runny] = apply_angle_part_s_parallel(res_2[runny].real, t0, alp, l)
        res_4 = f_a_for_p(el_nr,l,k_,z,nu_in,2,4)
        for runny in range(len(res_4)):
          res_4[runny] = apply_angle_part_s_parallel(res_4[runny].real,t0,alp, l)
        d_result[i] += (res_0)
        e_result[i] += ((res_2[0]+res_2[2]))
        f_result[i] += (res_2[1])
        e0_result[i] += ((res_4[0]+res_4[-1]))
        e1_result[i] += ((res_4[1]+res_4[-2]))
        e2_result[i] += ((res_4[2]))
        
        for _k in range(3):
          b_l   = f_p_el_for_p(el_nr, l,1, _k, z, nu_in, 2, 0)[0]
          c_0_l = f_p_el_for_p(el_nr, l,0, _k, z, nu_in, 2, 0)[0] 
          c_2_l = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 0)[0]
          d_l   = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 0)[0]
          b_l,c_0_l,c_2_l,d_l = apply_angle_part_p_parallel(b_l.real, c_0_l.real, c_2_l.real, d_l.real, _k, t0, alp, l)
          p0_result[i] += third * (b_l + c_0_l + c_2_l + d_l)
          
          b_l   = f_p_el_for_p(el_nr, l,1, _k, z, nu_in, 2, 2)
          c_0_l = f_p_el_for_p(el_nr, l,0, _k, z, nu_in, 2, 2) 
          c_2_l = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 2)
          d_l   = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 2)
          for runny in range(len(b_l)):
            b_l[runny],c_0_l[runny],c_2_l[runny],d_l[runny] = apply_angle_part_p_parallel(b_l[runny].real, 
                                                                                 c_0_l[runny].real,
                                                                                 c_2_l[runny].real,
                                                                                 d_l[runny].real,
                                                                                 _k, t0, alp, l)
          p2_0_result[i] += third * (b_l[0] + c_0_l[0] + c_2_l[0] + d_l[0] + b_l[2] + c_0_l[2] + c_2_l[2] + d_l[2])
          p2_1_result[i] += third * (b_l[1] + c_0_l[1] + c_2_l[1] + d_l[1])
          
          b_l   = f_p_el_for_p(el_nr,l,1,_k,z,nu_in,2,4)
          c_0_l = f_p_el_for_p(el_nr,l,0,_k,z,nu_in,2,4) 
          c_2_l = f_p_el_for_p(el_nr,l,2,_k,z,nu_in,2,4)
          d_l   = f_p_el_for_p(el_nr,l,2,_k,z,nu_in,2,4)
          for runny in range(len(b_l)):
            b_l[runny],c_0_l[runny],c_2_l[runny],d_l[runny] = apply_angle_part_p_parallel(b_l[runny].real,
                                                                                 c_0_l[runny].real,
                                                                                 c_2_l[runny].real,
                                                                                 d_l[runny].real,
                                                                                 _k, t0, alp, l)
          p4_0_result[i] += third * (b_l[0] + c_0_l[0] + c_2_l[0] + d_l[0] + b_l[-1] + c_0_l[-1] + c_2_l[-1] + d_l[-1])
          p4_1_result[i] += third * (b_l[1] + c_0_l[1] + c_2_l[1] + d_l[1] + b_l[-2] + c_0_l[-2] + c_2_l[-2] + d_l[-2])
          p4_2_result[i] += third * (b_l[2] + c_0_l[2] + c_2_l[2] + d_l[2])

          #b_l   =   f_b_for_p(el_nr,l,_k,z,nu_in,2,6)
          #c_0_l = f_c_0_for_p(el_nr,l,_k,z,nu_in,2,6) 
          #c_2_l = f_c_2_for_p(el_nr,l,_k,z,nu_in,2,6)
          #d_l   =   f_d_for_p(el_nr,l,_k,z,nu_in,2,6)
          #for runny in range(len(b_l)):
          #  b_l[runny],c_0_l[runny],c_2_l[runny],d_l[runny] = apply_angle_part_p_parallel(b_l[runny].real,
          #                                                                       c_0_l[runny].real,
          #                                                                       c_2_l[runny].real,
          #                                                                       d_l[runny].real,
          #                                                                       _k, t0, alp, l)
          #p6_0_result[i] += third * (b_l[0] + c_0_l[0] + c_2_l[0] + d_l[0] + b_l[-1] + c_0_l[-1] + c_2_l[-1] + d_l[-1])
          #p6_1_result[i] += third * (b_l[1] + c_0_l[1] + c_2_l[1] + d_l[1] + b_l[-2] + c_0_l[-2] + c_2_l[-2] + d_l[-2])
          #p6_1_result[i] += third * (b_l[2] + c_0_l[2] + c_2_l[2] + d_l[2] + b_l[-3] + c_0_l[-3] + c_2_l[-3] + d_l[-3])
          #p6_3_result[i] += third * (b_l[3] + c_0_l[3] + c_2_l[3] + d_l[3])
          
        #for p = 0 using the 0 angle applying a common angle function later
        b_l   = f_p_el_for_p(el_nr,l,1,1,z,nu_in,2,0)[0].real
        c_0_l = f_p_el_for_p(el_nr,l,0,0,z,nu_in,2,0)[0].real
        c_2_l = f_p_el_for_p(el_nr,l,2,2,z,nu_in,2,0)[0].real
        p0_result_M[i] += apply_angle_part_s_parallel(third * (  b_l   * alpha_coef(l, 1, 1, 0, 0)
                                                      + c_0_l * alpha_coef(l, 0, 0, 0, 0)
                                                      + c_2_l * (alpha_coef(l, 2, 2, 0, 0) + beta_bar_coef(l, 2, 2, 0, 0)))
                                             , t0, alp, 1)
        
        #for p = 2 using the 0 angle applying a common angle function later
        b_l   = f_p_el_for_p(el_nr,l,1,1,z,nu_in,2,2) 
        c_0_l = f_p_el_for_p(el_nr,l,0,0,z,nu_in,2,2) 
        c_2_l = f_p_el_for_p(el_nr,l,2,2,z,nu_in,2,2) 
        fac = [1,2,1]
        for runny in range(len(b_l)):
          res_2[runny] = apply_angle_part_s_parallel(third * (  b_l[runny].real * alpha_coef(l, 1, 1, 0, 0) 
                                               + c_0_l[runny].real * alpha_coef(l, 0, 0, 0, 0)
                                               + c_2_l[runny].real * (alpha_coef(l,2,2,0,0) + beta_bar_coef(l,2,2,0,0)))
                                            , t0, alp, fac[runny])
        p2_0_result_M[i] += ((res_2[0]+res_2[-1]))
        p2_1_result_M[i] += ((res_2[1]))
        
        #for p = 4 using the 0 angle applying a common angle function later
        b_l   = f_p_el_for_p(el_nr,l,1,1,z,nu_in,2,4) 
        c_0_l = f_p_el_for_p(el_nr,l,0,0,z,nu_in,2,4) 
        c_2_l = f_p_el_for_p(el_nr,l,2,2,z,nu_in,2,4) 
        fac = [1,2,3,2,1]
        for runny in range(len(b_l)):
          res_4[runny] = apply_angle_part_s_parallel(third * (  b_l[runny].real* alpha_coef(l,1,1,0,0) 
                                                     + c_0_l[runny].real * alpha_coef(l,0,0,0,0)
                                                     + c_2_l[runny].real * (alpha_coef(l,2,2,0,0) + beta_bar_coef(l,2,2,0,0)))
                                            , t0, alp, fac[runny])
        p4_0_result_M[i] += ((res_4[0]+res_4[-1]))
        p4_1_result_M[i] += ((res_4[1]+res_4[-2]))
        p4_2_result_M[i] += ((res_4[2]))
      
        ##for p = 6 using the 0 angle applying a common angle function later
        #b_l   = f_p_el_for_p(el_nr,l,1,1,z,nu_in,2,6) 
        #c_0_l = f_p_el_for_p(el_nr,l,0,0,z,nu_in,2,6) 
        #c_2_l = f_p_el_for_p(el_nr,l,2,2,z,nu_in,2,6) 
        #fac = [0,1,2,3,2,1,0]
        #for runny in range(len(b_l)):
        #  res_6[runny] = apply_angle_part_s_parallel(third * (  b_l[runny].real * alpha_coef(l,1,1,0,0) 
        #                                             + c_0_l[runny].real * alpha_coef(l,0,0,0,0)
        #                                             + c_2_l[runny].real * (alpha_coef(l,2,2,0,0) + beta_bar_coef(l,2,2,0,0)))
        #                                    , t0, alp, fac[runny])
        #p6_0_result_M[i] += ((res_6[0]+res_6[-1]))
        #p6_1_result_M[i] += ((res_6[1]+res_6[-2]))
        #p6_2_result_M[i] += ((res_6[2]+res_6[-3]))
        #p6_3_result_M[i] += ((res_6[3]))
        
  
      k_s[i] = a_result[i] + b_result[i] + c_result[i] + b0_result[i] + b1_result[i] + b2_result[i]
      l_s[i] = d_result[i] + e_result[i] + f_result[i] + e0_result[i] + e1_result[i] + e2_result[i]
      l_p[i] = p0_result[i] + p2_0_result[i] + p2_1_result[i] + p4_0_result[i] + p4_1_result[i] + p4_2_result[i]
      l_p_M[i] = p0_result_M[i] + p2_0_result_M[i] + p2_1_result_M[i] + p4_0_result_M[i] + p4_1_result_M[i] + p4_2_result_M[i]
      
      g_result[i] += apply_angle_part_s_parallel(f_s_1_hoenl(z)  *constant_factor     , t0, alp, 1)
      h_result[i] += apply_angle_part_s_parallel(f_s_2_1_hoenl(z)*constant_factor*x1_2, t0, alp, 1)
      i_result[i] += apply_angle_part_s_parallel(f_s_2_2_hoenl(z)*constant_factor*x1_2, t0, alp, 2)
  
      j_result[i] += apply_angle_part_s_parallel(f_s_1_EM(z)  *constant_factor      , t0, alp, 1)
      k_result[i] += apply_angle_part_s_parallel(f_s_2_1_EM(z)*constant_factor*x20_2, t0, alp, 1)
      l_result[i] += apply_angle_part_s_parallel(f_s_2_2_EM(z)*constant_factor*x20_2, t0, alp, 2)
  
      m_result[i] += apply_angle_part_s_parallel(f_p_1_EM(z)  *constant_factor      , t0, alp, 1)
      n_result[i] += apply_angle_part_s_parallel(f_p_2_1_EM(z)*constant_factor*x21_2, t0, alp, 1)
      o_result[i] += apply_angle_part_s_parallel(f_p_2_2_EM(z)*constant_factor*x21_2, t0, alp, 2)
    
    fig, axes = plt.subplots(2,2)
    axes[0,0].scatter(x,g_result,s=20,facecolors='none',edgecolors='b',marker='^',label="Hönl fs_1^(0)")
    axes[0,0].scatter(x,h_result,s=20,facecolors='none',edgecolors='g',marker='^',label="Hönl fs_1^(2)")
    axes[0,0].scatter(x,i_result,s=20,facecolors='none',edgecolors='r',marker='^',label="Hönl fs_2^(2)")
    axes[0,0].plot(x,a_result,color='b')
    axes[0,0].plot(x,b_result,color='g')
    axes[0,0].plot(x,c_result,color='r')
    axes[0,0].plot(x,b0_result,linestyle='dotted')
    axes[0,0].plot(x,b1_result,linestyle='dotted')
    axes[0,0].plot(x,b2_result,linestyle='dotted')
    axes[0,0].plot(x,k_s,color='black')
    axes[0,0].legend()
    axes[0,0].axhline(y=0,linestyle='dashed',color='gray')
  
    axes[1,0].scatter(x,j_result,s=20,facecolors='none',edgecolors='b',marker='^',label="EM fs_1^(0)")
    axes[1,0].scatter(x,k_result,s=20,facecolors='none',edgecolors='g',marker='^',label="EM fs_1^(2)")
    axes[1,0].scatter(x,l_result,s=20,facecolors='none',edgecolors='r',marker='^',label="EM fs_2^(2)")
    axes[1,0].plot(x,d_result,color='b')
    axes[1,0].plot(x,e_result,color='g')
    axes[1,0].plot(x,f_result,color='r')
    axes[1,0].plot(x,e0_result,linestyle='dotted')
    axes[1,0].plot(x,e1_result,linestyle='dotted')
    axes[1,0].plot(x,e2_result,linestyle='dotted')
    axes[1,0].plot(x,l_s,color='black')
    axes[1,0].legend()
    axes[1,0].axhline(y=0,linestyle='dashed',color='gray')
    
    axes[0,1].scatter(x,m_result,s=20,facecolors='none',edgecolors='b',marker='^',label="EM fp_1^(0)").legend_elements()
    axes[0,1].scatter(x,n_result,s=20,facecolors='none',edgecolors='g',marker='^',label="EM fp_1^(2)").legend_elements()
    axes[0,1].scatter(x,o_result,s=20,facecolors='none',edgecolors='r',marker='^',label="EM fp_2^(2)").legend_elements()
    ours1, = axes[0,1].plot(x,p0_result,color='b')
    ours2, = axes[0,1].plot(x,p2_0_result,color='g')
    ours3, = axes[0,1].plot(x,p2_1_result,color='r')
    ours4, = axes[0,1].plot(x,p4_0_result,linestyle='dotted')
    ours5, = axes[0,1].plot(x,p4_1_result,linestyle='dotted')
    ours6, = axes[0,1].plot(x,p4_2_result,linestyle='dotted')
    ours7, = axes[0,1].plot(x,p6_0_result,linestyle='dashed')
    ours8, = axes[0,1].plot(x,p6_1_result,linestyle='dashed')
    ours9, = axes[0,1].plot(x,p6_2_result,linestyle='dashed')
    ours10, = axes[0,1].plot(x,p6_3_result,linestyle='dashed')
    ours11, = axes[0,1].plot(x,l_p,color='black')
    axes[0,1].legend()
    axes[0,1].axhline(y=0,linestyle='dashed',color='gray')
    legend1 = axes[1,1].legend(
                    handles=[ours1,
                             ours2,ours3,
                             ours4,ours5,ours6,
                             ours7,ours8,ours9,ours10,
                             ours11],
                    labels=["0/0",
                            "0/2 + 2/0","1/1",
                            "0/4 + 4/0","1/3 + 3/1","2/2",
                            "0/6 + 6/0","1/5 + 5/1","2/4 + 4/2","3/3",
                            "sum"],
                    loc = 'upper left',
                    bbox_to_anchor=(-0.15,1.0))
    
    axes[1,1].plot(x,p0_result,color='b')
    axes[1,1].plot(x,p2_0_result,color='g')
    axes[1,1].plot(x,p2_1_result,color='r')
    axes[1,1].plot(x,p4_0_result,linestyle='dotted')
    axes[1,1].plot(x,p4_1_result,linestyle='dotted')
    axes[1,1].plot(x,p4_2_result,linestyle='dotted')
    axes[1,1].plot(x,l_p,color='k')
    axes[1,1].scatter(x,p0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 with M")
    axes[1,1].scatter(x,p2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")
    axes[1,1].scatter(x,p2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")
    axes[1,1].scatter(x,p4_0_result_M,s=20,facecolors='none',edgecolor='b',marker='o',label="p=4 0/4 with M")
    axes[1,1].scatter(x,p4_1_result_M,s=20,facecolors='none',edgecolor='orange',marker='o',label="p=4 1/3 with M")
    axes[1,1].scatter(x,p4_2_result_M,s=20,facecolors='none',edgecolor='g',marker='o',label="p=4 2/2 with M")
    axes[1,1].scatter(x,l_p_M,s=20,facecolors='none',edgecolor='k',marker='*',label="sum")
    axes[1,1].legend(loc='upper right')
    axes[1,1].add_artist(legend1)
  
    plt.subplots_adjust(left=0.025, bottom=0.04, right=1.0, top=1.0, wspace=0.15, hspace=0.05)
    fig.suptitle("alpha = {:4.2f}, theta_0 = {:4.2f}".format(alp,t0))

    print(p4_2_result/p4_2_result_M)
    print(p4_1_result/p4_1_result_M)
    plt.show()
    exit()
  
  if M_shell_test == True:
    print("Performing M parallel test")
    nu_in = 1.2 * get_ionization_energy_1s(el_nr) / h
    x = np.linspace(1.0001, 3.0001, steps)
    k_ = 0
    x1_2 = xn(nu_in, 1, el_nr, 0, 2)
    x20_2 = xn(nu_in, 2, el_nr, 0, 2)
    x21_2 = xn(nu_in, 2, el_nr, 1, 2)
    x30_2 = xn(nu_in, 3, el_nr, 0, 2)
    x31_2 = xn(nu_in, 3, el_nr, 1, 2)
    x32_2 = xn(nu_in, 3, el_nr, 2, 2)
    from time import time
    numpy.random.seed(int(np.ceil(time())))
    t0 = numpy.random.random(1)[0] * math.pi
    alp = numpy.random.random(1)[0] * math.pi * 2
    third = 1./3.
    res_2 = [0, 0, 0]
  
    if True:
      #K-shell
      a_result = np.zeros_like(x)
      b_result = np.zeros_like(x)
      c_result = np.zeros_like(x)
      k_s      = np.zeros_like(x)
      
      #L-shell
      l_s      = np.zeros_like(x)
      d_result = np.zeros_like(x)
      e_result = np.zeros_like(x)
      f_result = np.zeros_like(x)
    
      p0_result   = np.zeros_like(x)
      p2_0_result = np.zeros_like(x)
      p2_1_result = np.zeros_like(x)
      l_p         = np.zeros_like(x)
      
      p0_result_M   = np.zeros_like(x)
      p2_0_result_M = np.zeros_like(x)
      p2_1_result_M = np.zeros_like(x)
      l_p_M         = np.zeros_like(x)
      
      #M-shell
      Ms0_result = np.zeros_like(x)
      Ms1_result = np.zeros_like(x)
      Ms2_result = np.zeros_like(x)
      m_s        = np.zeros_like(x)
      
      M_p0_result   = np.zeros_like(x)
      M_p2_0_result = np.zeros_like(x)
      M_p2_1_result = np.zeros_like(x)
      m_p           = np.zeros_like(x)
      
      M_p0_result_M   = np.zeros_like(x)
      M_p2_0_result_M = np.zeros_like(x)
      M_p2_1_result_M = np.zeros_like(x)
      m_p_M           = np.zeros_like(x)

      M_d0_result = np.zeros_like(x)
      M_d2_0_result = np.zeros_like(x)
      M_d2_1_result = np.zeros_like(x)
      m_d = np.zeros_like(x)

      M_d0_result_M = np.zeros_like(x)
      M_d2_0_result_M = np.zeros_like(x)
      M_d2_1_result_M = np.zeros_like(x)
      m_d_M = np.zeros_like(x)
    
      #For comparison to Hönl etc
      g_result = np.zeros_like(x)
      h_result = np.zeros_like(x)
      i_result = np.zeros_like(x)
      j_result = np.zeros_like(x)
      k_result = np.zeros_like(x)
      l_result = np.zeros_like(x)
    
      m_result = np.zeros_like(x)
      n_result = np.zeros_like(x)
      o_result = np.zeros_like(x)
  
      x_result = np.zeros_like(x)
      y_result = np.zeros_like(x)
      z_result = np.zeros_like(x)
      WA_p0_result = np.zeros_like(x)
      WA_d0_result = np.zeros_like(x)
    for i,z in enumerate(x):
      for l in range(l_max+1):
        al00 = alpha_coef(l,0,0,0,0)
        al11 = alpha_coef(l,1,1,0,0)
        al22 = alpha_coef(l,2,2,0,0)
        al33 = alpha_coef(l,3,3,0,0)
        bbl22 = beta_bar_coef(l,2,2,0,0)
        bbl11 = beta_bar_coef(l,1,1,0,0)
        bbl33 = beta_bar_coef(l,3,3,0,0)

        #K-Shell
        res_0 = f_a_for_p(el_nr, l, k_, z, nu_in, 1, 0)[0]
        res_0 = apply_angle_part_s_parallel(res_0.real, t0, alp, l)
        res_2 = f_a_for_p(el_nr, l, k_, z, nu_in, 1, 2)
        for runny in range(len(res_2)):
          res_2[runny] = apply_angle_part_s_parallel(res_2[runny].real, t0, alp, l)
  
        a_result[i] += (res_0) /constant_factor
        b_result[i] += ((res_2[0]+res_2[2])) /constant_factor
        c_result[i] += (res_2[1]) /constant_factor
  
        #L-Shell
        #S-orbital
        res_0 = f_a_for_p(el_nr, l, k_, z, nu_in, 2, 0)[0]
        res_0 = apply_angle_part_s_parallel(res_0.real, t0, alp, l)
        res_2 = f_a_for_p(el_nr, l, k_, z, nu_in, 2, 2)
        for runny in range(len(res_2)):
          res_2[runny] = apply_angle_part_s_parallel(res_2[runny].real, t0, alp, l)
  
        d_result[i] += (res_0) /constant_factor
        e_result[i] += ((res_2[0]+res_2[2])) /constant_factor
        f_result[i] += (res_2[1]) /constant_factor
        
        #p-orbital
        for _k in range(3):
          b_l   = f_p_el_for_p(el_nr, l,1, _k, z, nu_in, 2, 0)[0]
          c_0_l = f_p_el_for_p(el_nr, l,0, _k, z, nu_in, 2, 0)[0]
          c_2_l = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 0)[0]
          d_l   = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 0)[0]
          b_l,c_0_l,c_2_l,d_l = apply_angle_part_p_parallel(b_l.real, 
                                                   c_0_l.real, 
                                                   c_2_l.real, 
                                                   d_l.real, 
                                                   _k, t0, alp, l)
          p0_result[i] += third * (b_l + c_0_l + c_2_l + d_l) /constant_factor
          
          b_l   = f_p_el_for_p(el_nr, l,1, _k, z, nu_in, 2, 2)
          c_0_l = f_p_el_for_p(el_nr, l,0, _k, z, nu_in, 2, 2)
          c_2_l = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 2)
          d_l   = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 2)
          for runny in range(len(b_l)):
            b_l[runny],c_0_l[runny],c_2_l[runny],d_l[runny] = apply_angle_part_p_parallel(b_l[runny].real, 
                                                                                 c_0_l[runny].real,
                                                                                 c_2_l[runny].real,
                                                                                 d_l[runny].real,
                                                                                 _k, t0, alp, l)
          p2_0_result[i] += third * (b_l[0] + c_0_l[0] + c_2_l[0] + d_l[0] + b_l[2] + c_0_l[2] + c_2_l[2] + d_l[2]) /constant_factor
          p2_1_result[i] += third * (b_l[1] + c_0_l[1] + c_2_l[1] + d_l[1]) /constant_factor
          
        #for p = 0 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l, 1,1, z, nu_in, 2, 0)[0].real
        c_0_l   = f_p_el_for_p(el_nr, l, 0,0, z, nu_in, 2, 0)[0].real
        c_2_l   = f_p_el_for_p(el_nr, l, 2,2, z, nu_in, 2, 0)[0].real
        p0_result_M[i] += apply_angle_part_s_parallel(third * (  b_l   * al11
                                                      + c_0_l * al00
                                                      + c_2_l * (al22 + bbl22)
                                                     )
                                             , t0, alp, 1) /constant_factor
        
        #for p = 2 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l, 1,1, z, nu_in, 2, 2)
        c_0_l   = f_p_el_for_p(el_nr, l, 0,0, z, nu_in, 2, 2)
        c_2_l   = f_p_el_for_p(el_nr, l, 2,2, z, nu_in, 2, 2)
        fac = [1,2,1]
        for num in range(3):
          res_2[num] = apply_angle_part_s_parallel(third * (  b_l[num].real * al11
                                                 + c_0_l[num].real * al00
                                                 + c_2_l[num].real * (al22+bbl22)
                                                  )
                                                 , t0, alp, fac[num])
        p2_0_result_M[i] += ((res_2[0]+res_2[-1])) /constant_factor
        p2_1_result_M[i] += ((res_2[1])) /constant_factor
        
        #M-Shell
        #S-orbital
        res_0 = f_a_for_p(el_nr, l, k_, z, nu_in, 3, 0)[0]
        res_0 = apply_angle_part_s_parallel(res_0.real, t0, alp, l)
        res_2 = f_a_for_p(el_nr, l, k_, z, nu_in, 3, 2)
        for runny in range(len(res_2)):
          res_2[runny] = apply_angle_part_s_parallel(res_2[runny].real, t0, alp, l)
  
        Ms0_result[i] += (res_0) /constant_factor
        Ms1_result[i] += ((res_2[0]+res_2[2])) /constant_factor
        Ms2_result[i] += (res_2[1]) /constant_factor
        
        #p-orbital
        #Calculate all angle dependant parts
        for _k in range(3):
          b_l   = f_p_el_for_p(el_nr, l,1, _k, z, nu_in, 3, 0)[0]
          c_0_l = f_p_el_for_p(el_nr, l,0, _k, z, nu_in, 3, 0)[0]
          c_2_l = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 3, 0)[0]
          d_l   = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 3, 0)[0]
          b_l,c_0_l,c_2_l,d_l = apply_angle_part_p_parallel(b_l.real, 
                                                   c_0_l.real, 
                                                   c_2_l.real, 
                                                   d_l.real, 
                                                   _k, t0, alp, l)
          M_p0_result[i] += third * (b_l + c_0_l + c_2_l + d_l) /constant_factor
          
          b_l   = f_p_el_for_p(el_nr, l,1, _k, z, nu_in, 3, 2)
          c_0_l = f_p_el_for_p(el_nr, l,0, _k, z, nu_in, 3, 2)
          c_2_l = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 3, 2)
          d_l   = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 3, 2)
          for runny in range(len(b_l)):
            b_l[runny],c_0_l[runny],c_2_l[runny],d_l[runny] = apply_angle_part_p_parallel(b_l[runny].real, 
                                                                                 c_0_l[runny].real,
                                                                                 c_2_l[runny].real,
                                                                                 d_l[runny].real,
                                                                                 _k, t0, alp, l)
          M_p2_0_result[i] += third * (b_l[0] + c_0_l[0] + c_2_l[0] + d_l[0] + b_l[2] + c_0_l[2] + c_2_l[2] + d_l[2]) /constant_factor
          M_p2_1_result[i] += third * (b_l[1] + c_0_l[1] + c_2_l[1] + d_l[1]) /constant_factor
          
        #for p = 0 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l, 1,1, z, nu_in, 3, 0)[0].real
        c_0_l   = f_p_el_for_p(el_nr, l, 0,0, z, nu_in, 3, 0)[0].real
        c_2_l   = f_p_el_for_p(el_nr, l, 2,2, z, nu_in, 3, 0)[0].real
        M_p0_result_M[i] += apply_angle_part_s_parallel(third * (  b_l   * al11
                                                        + c_0_l * al00
                                                        + c_2_l * (al22 + bbl22)
                                                        )
                                             , t0, alp, 1) /constant_factor
        
        #for p = 2 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l, 1,1, z, nu_in, 3, 2)
        c_0_l   = f_p_el_for_p(el_nr, l, 0,0, z, nu_in, 3, 2)
        c_2_l   = f_p_el_for_p(el_nr, l, 2,2, z, nu_in, 3, 2)
        fac = [1,2,1]
        for num in range(3):
          res_2[num] = apply_angle_part_s_parallel(third * (  b_l[num].real * al11
                                                 + c_0_l[num].real * al00
                                                 + c_2_l[num].real * (al22+bbl22)
                                                  )
                                                 , t0, alp, fac[num])
        M_p2_0_result_M[i] += ((res_2[0]+res_2[-1])) /constant_factor
        M_p2_1_result_M[i] += ((res_2[1])) /constant_factor

        #d-orbital
        #Calculate all angle dependant parts
        ms = [0,2,2,1,3,4,1,3]
        for _k in range(5):
          res = [0,0,0,0,0,0,0,0]
          for index in range(len(ms)):
            res[index] = f_d_el_for_p(el_nr, l, ms[index], _k, z, nu_in, 3, 0)[0].real
          res = apply_angle_part_d_parallel(res, _k, t0, alp, l)
          M_d0_result[i] += 0.2 * sum(res) /constant_factor

          res = [[0,0,0,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0]]
          for index in range(len(ms)):
            temp = f_d_el_for_p(el_nr, l, ms[index], _k, z, nu_in, 3, 2)
            num = 0
            for v in temp:
              res[num][index] = v.real
              num +=1

          for runny in range(len(res)):
            res[runny] = apply_angle_part_d_parallel(res[runny], _k, t0, alp, l)
          M_d2_0_result[i] += 0.2 * (sum(res[0]) + sum(res[2])) /constant_factor
          M_d2_1_result[i] += 0.2 * (sum(res[1])) /constant_factor

        #for p = 0 using the 0 angle applying a common angle function later
        ms = [0,1,2,3,4]
        mults = [al00,bbl11+al11,al22+bbl22,bbl33+al33,al11]
        res = 0
        for nummy in range(len(ms)):
          r = f_d_el_for_p(el_nr, l, ms[nummy], ms[nummy], z, nu_in, 3, 0)
          res += r[0].real * mults[nummy]
        M_d0_result_M[i] += apply_angle_part_s_parallel(0.2 * res, t0, alp, 1)  /constant_factor
        
        #for p = 2 using the 0 angle applying a common angle function later
        res_2 = [0,0,0]
        fac = [1,2,1]
        for nummy in range(len(ms)):
          r = f_d_el_for_p(el_nr, l, ms[nummy], ms[nummy], z, nu_in, 3, 2)
          for index in range(3):
            res_2[index] += r[index].real * mults[nummy]
        for num in range(3):
          res_2[num] = apply_angle_part_s_parallel(0.2 * res_2[num], t0, alp, fac[num])
        M_d2_0_result_M[i] += ((res_2[0]+res_2[-1])) /constant_factor
        M_d2_1_result_M[i] += ((res_2[1])) /constant_factor

      k_s[i]   = (a_result[i] + b_result[i] + c_result[i]             )
      l_s[i]   = (d_result[i] + e_result[i] + f_result[i]             )
      l_p[i]   = (p0_result[i] + p2_0_result[i] + p2_1_result[i]      )
      l_p_M[i] = (p0_result_M[i] + p2_0_result_M[i] + p2_1_result_M[i])
      
      m_s[i]   = (Ms0_result[i] + Ms1_result[i] + Ms2_result[i]             ) 
      m_p[i]   = (M_p0_result[i] + M_p2_0_result[i] + M_p2_1_result[i]      )
      m_p_M[i] = (M_p0_result_M[i] + M_p2_0_result_M[i] + M_p2_1_result_M[i])

      m_d[i]   = (M_d0_result[i] + M_d2_0_result[i] + M_d2_1_result[i]      )
      m_d_M[i] = (M_d0_result_M[i] + M_d2_0_result_M[i] + M_d2_1_result_M[i])
      
      g_result[i] += apply_angle_part_s_parallel(f_s_1_hoenl(z)       ,t0, alp, 1)
      h_result[i] += apply_angle_part_s_parallel(f_s_2_1_hoenl(z)*x1_2,t0, alp, 1)
      i_result[i] += apply_angle_part_s_parallel(f_s_2_2_hoenl(z)*x1_2,t0, alp, 2)
  
      j_result[i] += apply_angle_part_s_parallel(f_s_1_EM(z)        , t0, alp, 1)
      k_result[i] += apply_angle_part_s_parallel(f_s_2_1_EM(z)*x20_2, t0, alp, 1)
      l_result[i] += apply_angle_part_s_parallel(f_s_2_2_EM(z)*x20_2, t0, alp, 2)
  
      m_result[i] += apply_angle_part_s_parallel(f_p_1_EM(z)        , t0, alp, 1)
      n_result[i] += apply_angle_part_s_parallel(f_p_2_1_EM(z)*x21_2, t0, alp, 1)
      o_result[i] += apply_angle_part_s_parallel(f_p_2_2_EM(z)*x21_2, t0, alp, 2)
  
      x_result[i] += apply_angle_part_s_parallel(f_s_1_WA(z)        ,t0, alp, 1)
      y_result[i] += apply_angle_part_s_parallel(f_s_2_1_WA(z)*x30_2,t0, alp, 1)
      z_result[i] += apply_angle_part_s_parallel(f_s_2_2_WA(z)*x30_2,t0, alp, 2)

      WA_p0_result[i] += apply_angle_part_s_parallel(f_p_1_WA(z) ,t0, alp, 1)
      WA_d0_result[i] += apply_angle_part_s_parallel(f_d_1_WA(z) ,t0, alp, 1)
    
    fig, axes = plt.subplots(3,3)
    #K-shell s-orbitals
    axes[0,0].scatter(x,g_result,s=20,facecolors='none',edgecolors='b',marker='^',label="Hönl fs_1^(0)")
    axes[0,0].scatter(x,h_result,s=20,facecolors='none',edgecolors='g',marker='^',label="Hönl fs_1^(2)")
    axes[0,0].scatter(x,i_result,s=20,facecolors='none',edgecolors='r',marker='^',label="Hönl fs_2^(2)")
    axes[0,0].plot(x,a_result,color='b')
    axes[0,0].plot(x,b_result,color='g')
    axes[0,0].plot(x,c_result,color='r')
    axes[0,0].plot(x,k_s,color='black')
    axes[0,0].legend()
    axes[0,0].axhline(y=0,linestyle='dashed',color='gray')
    axes[0,0].set_title("K-shell s-electrons", y=1.0, pad=-14)
    #L-shell s-orbital
    axes[1,0].scatter(x,j_result,s=20,facecolors='none',edgecolors='b',marker='^',label="EM fs_1^(0)")
    axes[1,0].scatter(x,k_result,s=20,facecolors='none',edgecolors='g',marker='^',label="EM fs_1^(2)")
    axes[1,0].scatter(x,l_result,s=20,facecolors='none',edgecolors='r',marker='^',label="EM fs_2^(2)")
    axes[1,0].plot(x,d_result,color='b')
    axes[1,0].plot(x,e_result,color='g')
    axes[1,0].plot(x,f_result,color='r')
    axes[1,0].plot(x,l_s,color='black')
    axes[1,0].legend()
    axes[1,0].axhline(y=0,linestyle='dashed',color='gray')
    axes[1,0].set_title("L-shell s-electrons", y=1.0, pad=-14)
    #M-shell s-orbital
    axes[2,0].scatter(x,x_result,s=20,facecolors='none',edgecolors='b',marker='^',label="WA fs_1^(0)")
    axes[2,0].scatter(x,y_result,s=20,facecolors='none',edgecolors='g',marker='^',label="WA fs_1^(2)")
    axes[2,0].scatter(x,z_result,s=20,facecolors='none',edgecolors='r',marker='^',label="WA fs_2^(2)")
    axes[2,0].plot(x,Ms0_result,color='b')
    axes[2,0].plot(x,Ms1_result,color='g')
    axes[2,0].plot(x,Ms2_result,color='r')
    axes[2,0].plot(x,m_s,color='black')
    axes[2,0].legend()
    axes[2,0].axhline(y=0,linestyle='dashed',color='gray')
    axes[2,0].set_title("M-shell s-electrons", y=1.0, pad=-14)
    
    ours1, = axes[1,1].plot(x,p0_result,color='b')
    ours2, = axes[1,1].plot(x,p2_0_result,color='g')
    ours3, = axes[1,1].plot(x,p2_1_result,color='r')
    ours11, = axes[1,1].plot(x,l_p,color='black')
    legend1 = axes[1,1].legend(
                    handles=[ours1,
                             ours2,ours3,
                             ours11],
                    labels=["0/0",
                            "0/2 + 2/0","1/1",
                            "sum"],
                    loc = 'upper left',
                    bbox_to_anchor=(0.5,1.5))
    
    axes[1,1].plot(x,p0_result,  color='b')
    axes[1,1].plot(x,p2_0_result,color='g')
    axes[1,1].plot(x,p2_1_result,color='r')
    axes[1,1].scatter(x,m_result,s=40,facecolors='none',edgecolors='b',marker='o',label="EM fp_1^(0)")
    axes[1,1].scatter(x,n_result,s=40,facecolors='none',edgecolors='g',marker='o',label="EM fp_1^(2)")
    axes[1,1].scatter(x,o_result,s=40,facecolors='none',edgecolors='r',marker='o',label="EM fp_2^(2)")
    axes[1,1].plot(x,l_p,color='black')
    axes[1,1].scatter(x,p0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")
    axes[1,1].scatter(x,p2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")
    axes[1,1].scatter(x,p2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")
    axes[1,1].scatter(x,l_p_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")
    axes[1,1].legend()
    axes[1,1].add_artist(legend1)
    axes[1,1].axhline(y=0,linestyle='dashed',color='gray')
    axes[1,1].set_title("L-shell p-electrons", y=1.0, pad=-14)
    
    axes[2,1].plot(x,M_p0_result,  color='b')
    axes[2,1].plot(x,M_p2_0_result,color='g')
    axes[2,1].plot(x,M_p2_1_result,color='r')
    axes[2,1].plot(x,m_p,color='black')
    axes[2,1].scatter(x,M_p0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")
    axes[2,1].scatter(x,M_p2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")
    axes[2,1].scatter(x,M_p2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")
    axes[2,1].scatter(x,WA_p0_result,s=40,facecolors='none',edgecolor='b',marker='o',label="Wagenfeld p-els D")
    axes[2,1].scatter(x,m_p_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")
    axes[2,1].legend()
    axes[2,1].axhline(y=0,linestyle='dashed',color='gray')
    axes[2,1].set_title("M-shell p-electrons", y=1.0, pad=-14)

    axes[0,1].set_axis_off()
    axes[0,2].set_axis_off()
    axes[1,2].set_axis_off()

    axes[2,2].plot(x,M_d0_result,  color='b')
    axes[2,2].plot(x,M_d2_0_result,color='g')
    axes[2,2].plot(x,M_d2_1_result,color='r')
    axes[2,2].plot(x,m_d,color='black')
    axes[2,2].scatter(x,M_d0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")
    axes[2,2].scatter(x,M_d2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")
    axes[2,2].scatter(x,M_d2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")
    axes[2,2].scatter(x,WA_d0_result,s=40,facecolors='none',edgecolor='b',marker='o',label="Wagenfeld d-els D")
    axes[2,2].scatter(x,m_d_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")
    axes[2,2].legend()
    axes[2,2].axhline(y=0,linestyle='dashed',color='gray')
    axes[2,2].set_title("M-shell d-electrons", y=1.0, pad=-14)
  
    plt.subplots_adjust(left=0.04, bottom=0.04, right=1.0, top=0.95, wspace=0.15, hspace=0.1)
    fig.suptitle("PARA alpha = {:4.2f}, theta_0 = {:4.2f}".format(alp,t0))
    plt.show()
    exit()

  if M_shell_test_orthogonal == True:
    print("Performing M Orthogonal test")
    nu_in = 1.2 * get_ionization_energy_1s(el_nr) / h
    x = np.linspace(1.0001, 3.0001, steps)
    k_ = 0
    x1_2 = xn(nu_in, 1, el_nr, 0, 2)
    x20_2 = xn(nu_in, 2, el_nr, 0, 2)
    x21_2 = xn(nu_in, 2, el_nr, 1, 2)
    x30_2 = xn(nu_in, 3, el_nr, 0, 2)
    x31_2 = xn(nu_in, 3, el_nr, 1, 2)
    x32_2 = xn(nu_in, 3, el_nr, 2, 2)
    from time import time
    numpy.random.seed(int(np.ceil(time())))
    t0 = numpy.random.random(1)[0] * math.pi
    #t0 = math.pi / 2
    alp = numpy.random.random(1)[0] * math.pi * 2
    #alp = math.pi
    third = 1./3.
    res_2 = [0, 0, 0]
  
    if True:
      #K-shell
      a_result = np.zeros_like(x)
      b_result = np.zeros_like(x)
      c_result = np.zeros_like(x)
      k_s      = np.zeros_like(x)
      
      #L-shell
      l_s      = np.zeros_like(x)
      d_result = np.zeros_like(x)
      e_result = np.zeros_like(x)
      f_result = np.zeros_like(x)
    
      p0_result   = np.zeros_like(x)
      p2_0_result = np.zeros_like(x)
      p2_1_result = np.zeros_like(x)
      l_p         = np.zeros_like(x)
      
      p0_result_M   = np.zeros_like(x)
      p2_0_result_M = np.zeros_like(x)
      p2_1_result_M = np.zeros_like(x)
      l_p_M         = np.zeros_like(x)
      
      #M-shell
      Ms0_result = np.zeros_like(x)
      Ms1_result = np.zeros_like(x)
      Ms2_result = np.zeros_like(x)
      m_s        = np.zeros_like(x)
      
      M_p0_result   = np.zeros_like(x)
      M_p2_0_result = np.zeros_like(x)
      M_p2_1_result = np.zeros_like(x)
      m_p           = np.zeros_like(x)
      
      M_p0_result_M   = np.zeros_like(x)
      M_p2_0_result_M = np.zeros_like(x)
      M_p2_1_result_M = np.zeros_like(x)
      m_p_M           = np.zeros_like(x)

      M_d0_result = np.zeros_like(x)
      M_d2_0_result = np.zeros_like(x)
      M_d2_1_result = np.zeros_like(x)
      m_d = np.zeros_like(x)

      M_d0_result_M = np.zeros_like(x)
      M_d2_0_result_M = np.zeros_like(x)
      M_d2_1_result_M = np.zeros_like(x)
      m_d_M = np.zeros_like(x)
    
      #For comparison to Hönl etc
      g_result = np.zeros_like(x)
      h_result = np.zeros_like(x)
      i_result = np.zeros_like(x)
      j_result = np.zeros_like(x)
      k_result = np.zeros_like(x)
      l_result = np.zeros_like(x)
    
      m_result = np.zeros_like(x)
      n_result = np.zeros_like(x)
      o_result = np.zeros_like(x)
  
      x_result = np.zeros_like(x)
      y_result = np.zeros_like(x)
      z_result = np.zeros_like(x)
      WA_p0_result = np.zeros_like(x)
      WA_d0_result = np.zeros_like(x)
    for i,z in enumerate(x):
      for l in range(l_max+1):
        al00 = alpha_coef(l,0,0,0,0)
        al11 = alpha_coef(l,1,1,0,0)
        al22 = alpha_coef(l,2,2,0,0)
        al33 = alpha_coef(l,3,3,0,0)
        bbl22 = beta_bar_coef(l,2,2,0,0)
        bbl11 = beta_bar_coef(l,1,1,0,0)
        bbl33 = beta_bar_coef(l,3,3,0,0)

        #K-Shell
        res_0 = f_a_for_p(el_nr, l, k_, z, nu_in, 1, 0)[0]
        res_0 = apply_angle_part_s_orthogonal(res_0.real, t0, alp, l)
        res_2 = f_a_for_p(el_nr, l, k_, z, nu_in, 1, 2)
        for runny in range(len(res_2)):
          res_2[runny] = apply_angle_part_s_orthogonal(res_2[runny].real, t0, alp, l)
  
        a_result[i] += (res_0)
        b_result[i] += ((res_2[0]+res_2[2]))
        c_result[i] += (res_2[1])
  
        #L-Shell
        #S-orbital
        res_0 = f_a_for_p(el_nr, l, k_, z, nu_in, 2, 0)[0]
        res_0 = apply_angle_part_s_orthogonal(res_0.real, t0, alp, l)
        res_2 = f_a_for_p(el_nr, l, k_, z, nu_in, 2, 2)
        for runny in range(len(res_2)):
          res_2[runny] = apply_angle_part_s_orthogonal(res_2[runny].real, t0, alp, l)
  
        d_result[i] += (res_0)
        e_result[i] += ((res_2[0]+res_2[2]))
        f_result[i] += (res_2[1])
        
        #p-orbital
        for _k in range(3):
          b_l   = f_p_el_for_p(el_nr, l,1, _k, z, nu_in, 2, 0)[0]
          c_0_l = f_p_el_for_p(el_nr, l,0, _k, z, nu_in, 2, 0)[0]
          c_2_l = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 0)[0]
          d_l   = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 0)[0]
          b_l,c_0_l,c_2_l,d_l = apply_angle_part_p_orthogonal(b_l.real, 
                                                   c_0_l.real, 
                                                   c_2_l.real, 
                                                   d_l.real, 
                                                   _k, t0, alp, l)
          p0_result[i] += third * (b_l + c_0_l + c_2_l + d_l)
          
          b_l   = f_p_el_for_p(el_nr, l,1, _k, z, nu_in, 2, 2)
          c_0_l = f_p_el_for_p(el_nr, l,0, _k, z, nu_in, 2, 2)
          c_2_l = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 2)
          d_l   = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 2, 2)
          for runny in range(len(b_l)):
            b_l[runny],c_0_l[runny],c_2_l[runny],d_l[runny] = apply_angle_part_p_orthogonal(b_l[runny].real, 
                                                                                 c_0_l[runny].real,
                                                                                 c_2_l[runny].real,
                                                                                 d_l[runny].real,
                                                                                 _k, t0, alp, l)
          p2_0_result[i] += third * (b_l[0] + c_0_l[0] + c_2_l[0] + d_l[0] + b_l[2] + c_0_l[2] + c_2_l[2] + d_l[2])
          p2_1_result[i] += third * (b_l[1] + c_0_l[1] + c_2_l[1] + d_l[1])
          
        #for p = 0 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l,1, 1, z, nu_in, 2, 0)[0].real
        c_0_l   = f_p_el_for_p(el_nr, l,0, 0, z, nu_in, 2, 0)[0].real
        c_2_l   = f_p_el_for_p(el_nr, l,2, 2, z, nu_in, 2, 0)[0].real
        #d_l     =   f_d_for_p(el_nr, l, 2, z, nu_in, 2, 0)[0].real
        p0_result_M[i] += apply_angle_part_s_orthogonal(third * (  b_l   * al11
                                                      + c_0_l * al00
                                                      + c_2_l * (al22 + bbl22)
                                                     )
                                                     # + d_l   * bbl22)
                                             , t0, alp, 1)
        
        #for p = 2 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l,1, 1, z, nu_in, 2, 2)
        c_0_l   = f_p_el_for_p(el_nr, l,0, 0, z, nu_in, 2, 2)
        c_2_l   = f_p_el_for_p(el_nr, l,2, 2, z, nu_in, 2, 2)
        #d_l     =   f_d_for_p(el_nr, l, 2, z, nu_in, 2, 2)
        fac = [1,2,1]
        for num in range(3):
          res_2[num] = apply_angle_part_s_orthogonal(third * (  b_l[num].real * al11
                                                 + c_0_l[num].real * al00
                                                 + c_2_l[num].real * (al22+bbl22)
                                                 #+ d_l[num].real   * bbl22)
                                                  )
                                                 , t0, alp, fac[num])
        p2_0_result_M[i] += ((res_2[0]+res_2[-1]))
        p2_1_result_M[i] += ((res_2[1]))
        
        #M-Shell
        #S-orbital
        res_0 = f_a_for_p(el_nr, l, k_, z, nu_in, 3, 0)[0]
        res_0 = apply_angle_part_s_orthogonal(res_0.real, t0, alp, l)
        res_2 = f_a_for_p(el_nr, l, k_, z, nu_in, 3, 2)
        for runny in range(len(res_2)):
          res_2[runny] = apply_angle_part_s_orthogonal(res_2[runny].real, t0, alp, l)
  
        Ms0_result[i] += (res_0)
        Ms1_result[i] += ((res_2[0]+res_2[2]))
        Ms2_result[i] += (res_2[1])
        
        #p-orbital
        #Calculate all angle dependant parts
        for _k in range(3):
          b_l   = f_p_el_for_p(el_nr, l,1, _k, z, nu_in, 3, 0)[0]
          c_0_l = f_p_el_for_p(el_nr, l,0, _k, z, nu_in, 3, 0)[0]
          c_2_l = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 3, 0)[0]
          d_l   = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 3, 0)[0]
          b_l,c_0_l,c_2_l,d_l = apply_angle_part_p_orthogonal(b_l.real, 
                                                   c_0_l.real, 
                                                   c_2_l.real, 
                                                   d_l.real, 
                                                   _k, t0, alp, l)
          M_p0_result[i] += third * (b_l + c_0_l + c_2_l + d_l)
          
          b_l   = f_p_el_for_p(el_nr, l,1, _k, z, nu_in, 3, 2)
          c_0_l = f_p_el_for_p(el_nr, l,0, _k, z, nu_in, 3, 2)
          c_2_l = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 3, 2)
          d_l   = f_p_el_for_p(el_nr, l,2, _k, z, nu_in, 3, 2)
          for runny in range(len(b_l)):
            b_l[runny],c_0_l[runny],c_2_l[runny],d_l[runny] = apply_angle_part_p_orthogonal(b_l[runny].real, 
                                                                                 c_0_l[runny].real,
                                                                                 c_2_l[runny].real,
                                                                                 d_l[runny].real,
                                                                                 _k, t0, alp, l)
          M_p2_0_result[i] += third * (b_l[0] + c_0_l[0] + c_2_l[0] + d_l[0] + b_l[2] + c_0_l[2] + c_2_l[2] + d_l[2])
          M_p2_1_result[i] += third * (b_l[1] + c_0_l[1] + c_2_l[1] + d_l[1])
          
        #for p = 0 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l,1, 1, z, nu_in, 3, 0)[0].real
        c_0_l   = f_p_el_for_p(el_nr, l,0, 0, z, nu_in, 3, 0)[0].real
        c_2_l   = f_p_el_for_p(el_nr, l,2, 2, z, nu_in, 3, 0)[0].real
        #d_l     =   f_d_for_p(el_nr, l, 2, z, nu_in, 3, 0)[0].real
        M_p0_result_M[i] += apply_angle_part_s_orthogonal(third * (  b_l   * al11
                                                        + c_0_l * al00
                                                        + c_2_l * (al22 + bbl22)
                                                        #+ d_l   * bbl22)
                                                        )
                                             , t0, alp, 1)
        
        #for p = 2 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l,1, 1, z, nu_in, 3, 2)
        c_0_l   = f_p_el_for_p(el_nr, l,0, 0, z, nu_in, 3, 2)
        c_2_l   = f_p_el_for_p(el_nr, l,2, 2, z, nu_in, 3, 2)
        #d_l     =   f_d_for_p(el_nr, l, 2, z, nu_in, 3, 2)
        fac = [1,2,1]
        for num in range(3):
          res_2[num] = apply_angle_part_s_orthogonal(third * (  b_l[num].real * al11
                                                 + c_0_l[num].real * al00
                                                 + c_2_l[num].real * (al22+bbl22)
                                                 #+ d_l[num].real   * bbl22)
                                                  )
                                                 , t0, alp, fac[num])
        M_p2_0_result_M[i] += ((res_2[0]+res_2[-1]))
        M_p2_1_result_M[i] += ((res_2[1]))

        #d-orbital
        #Calculate all angle dependant parts
        ms = [0,2,2,1,3,4,1,3]
        for _k in range(5):
          #f_0_l = f_d_el_for_p(el_nr, l, 0, _k, z, nu_in, 3, 2)
          #f_2_l = f_d_el_for_p(el_nr, l, 2, _k, z, nu_in, 3, 2)
          #e_l   = f_d_el_for_p(el_nr, l, 2, _k, z, nu_in, 3, 2)
          #g_1_l = f_d_el_for_p(el_nr, l, 1, _k, z, nu_in, 3, 2)
          #g_3_l = f_d_el_for_p(el_nr, l, 3, _k, z, nu_in, 3, 2)
          #h_l   = f_d_el_for_p(el_nr, l, 4, _k, z, nu_in, 3, 2)
          #i_1_l = f_d_el_for_p(el_nr, l, 1, _k, z, nu_in, 3, 2)
          #i_3_l = f_d_el_for_p(el_nr, l, 3, _k, z, nu_in, 3, 2)
          res = [0,0,0,0,0,0,0,0]
          for index in range(len(ms)):
            res[index] = f_d_el_for_p(el_nr, l, ms[index], _k, z, nu_in, 3, 0)[0].real
          res = apply_angle_part_d_orthogonal(res, _k, t0, alp, l)
          M_d0_result[i] += 0.2 * sum(res)

          res = [[0,0,0,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0],
                 [0,0,0,0,0,0,0,0]]
          for index in range(len(ms)):
            temp = f_d_el_for_p(el_nr, l, ms[index], _k, z, nu_in, 3, 2)
            num = 0
            for v in temp:
              res[num][index] = v.real
              num +=1

          for runny in range(len(res)):
            res[runny] = apply_angle_part_d_orthogonal(res[runny], _k, t0, alp, l)
          M_d2_0_result[i] += 0.2 * (sum(res[0]) + sum(res[2]))
          M_d2_1_result[i] += 0.2 * (sum(res[1]))

        #for p = 0 using the 0 angle applying a common angle function later
        ms = [0,1,2,3,4]
        mults = [al00,bbl11+al11,al22+bbl22,bbl33+al33,al11]
        res = 0
        for nummy in range(len(ms)):
          r = f_d_el_for_p(el_nr, l, ms[nummy], ms[nummy], z, nu_in, 3, 0)
          res += r[0].real * mults[nummy]
        M_d0_result_M[i] += apply_angle_part_s_orthogonal(0.2 * res, t0, alp, 1)
        
        #for p = 2 using the 0 angle applying a common angle function later
        res_2 = [0,0,0]
        fac = [1,2,1]
        for nummy in range(len(ms)):
          r = f_d_el_for_p(el_nr, l, ms[nummy], ms[nummy], z, nu_in, 3, 2)
          for index in range(3):
            res_2[index] += r[index].real * mults[nummy]
        for num in range(3):
          res_2[num] = apply_angle_part_s_orthogonal(0.2 * res_2[num], t0, alp, fac[num])
        M_d2_0_result_M[i] += ((res_2[0]+res_2[-1]))
        M_d2_1_result_M[i] += ((res_2[1]))

      k_s[i] = a_result[i] + b_result[i] + c_result[i]
      l_s[i] = d_result[i] + e_result[i] + f_result[i]
      l_p[i] = p0_result[i] + p2_0_result[i] + p2_1_result[i]
      l_p_M[i] = p0_result_M[i] + p2_0_result_M[i] + p2_1_result_M[i]
      
      m_s[i] = Ms0_result[i] + Ms1_result[i] + Ms2_result[i]
      m_p[i] = M_p0_result[i] + M_p2_0_result[i] + M_p2_1_result[i]
      m_p_M[i] = M_p0_result_M[i] + M_p2_0_result_M[i] + M_p2_1_result_M[i]

      m_d[i] = M_d0_result[i] + M_d2_0_result[i] + M_d2_1_result[i]
      m_d_M[i] = M_d0_result_M[i] + M_d2_0_result_M[i] + M_d2_1_result_M[i]
      
      g_result[i] += apply_angle_part_s_orthogonal(f_s_1_hoenl(z)  *constant_factor     ,t0, alp, 1)
      h_result[i] += apply_angle_part_s_orthogonal(f_s_2_1_hoenl(z)*constant_factor*x1_2,t0, alp, 1)
      i_result[i] += apply_angle_part_s_orthogonal(f_s_2_2_hoenl(z)*constant_factor*x1_2,t0, alp, 2)
  
      j_result[i] += apply_angle_part_s_orthogonal(f_s_1_EM(z)  * constant_factor     , t0, alp, 1)
      k_result[i] += apply_angle_part_s_orthogonal(f_s_2_1_EM(z)* constant_factor*x20_2, t0, alp, 1)
      l_result[i] += apply_angle_part_s_orthogonal(f_s_2_2_EM(z)* constant_factor*x20_2, t0, alp, 2)
  
      m_result[i] += apply_angle_part_s_orthogonal(f_p_1_EM(z)  * constant_factor     , t0, alp, 1)
      n_result[i] += apply_angle_part_s_orthogonal(f_p_2_1_EM(z)* constant_factor*x21_2, t0, alp, 1)
      o_result[i] += apply_angle_part_s_orthogonal(f_p_2_2_EM(z)* constant_factor*x21_2, t0, alp, 2)
  
      x_result[i] += apply_angle_part_s_orthogonal(f_s_1_WA(z)  * constant_factor     ,t0, alp, 1)
      y_result[i] += apply_angle_part_s_orthogonal(f_s_2_1_WA(z)* constant_factor*x30_2,t0, alp, 1)
      z_result[i] += apply_angle_part_s_orthogonal(f_s_2_2_WA(z)* constant_factor*x30_2,t0, alp, 2)

      WA_p0_result[i] += apply_angle_part_s_orthogonal(f_p_1_WA(z) * constant_factor     ,t0, alp, 1)
      WA_d0_result[i] += apply_angle_part_s_orthogonal(f_d_1_WA(z) * constant_factor     ,t0, alp, 1)
    
    fig, axes = plt.subplots(3,3)
    #K-shell s-orbitals
    axes[0,0].scatter(x,g_result,s=20,facecolors='none',edgecolors='b',marker='^',label="Hönl fs_1^(0)")
    axes[0,0].scatter(x,h_result,s=20,facecolors='none',edgecolors='g',marker='^',label="Hönl fs_1^(2)")
    axes[0,0].scatter(x,i_result,s=20,facecolors='none',edgecolors='r',marker='^',label="Hönl fs_2^(2)")
    axes[0,0].plot(x,a_result,color='b')
    axes[0,0].plot(x,b_result,color='g')
    axes[0,0].plot(x,c_result,color='r')
    axes[0,0].plot(x,k_s,color='black')
    axes[0,0].legend()
    axes[0,0].axhline(y=0,linestyle='dashed',color='gray')
    axes[0,0].set_title("K-shell s-electrons", y=1.0, pad=-14)
    #L-shell s-orbital
    axes[1,0].scatter(x,j_result,s=20,facecolors='none',edgecolors='b',marker='^',label="EM fs_1^(0)")
    axes[1,0].scatter(x,k_result,s=20,facecolors='none',edgecolors='g',marker='^',label="EM fs_1^(2)")
    axes[1,0].scatter(x,l_result,s=20,facecolors='none',edgecolors='r',marker='^',label="EM fs_2^(2)")
    axes[1,0].plot(x,d_result,color='b')
    axes[1,0].plot(x,e_result,color='g')
    axes[1,0].plot(x,f_result,color='r')
    axes[1,0].plot(x,l_s,color='black')
    axes[1,0].legend()
    axes[1,0].axhline(y=0,linestyle='dashed',color='gray')
    axes[1,0].set_title("L-shell s-electrons", y=1.0, pad=-14)
    #M-shell s-orbital
    axes[2,0].scatter(x,x_result,s=20,facecolors='none',edgecolors='b',marker='^',label="WA fs_1^(0)")
    axes[2,0].scatter(x,y_result,s=20,facecolors='none',edgecolors='g',marker='^',label="WA fs_1^(2)")
    axes[2,0].scatter(x,z_result,s=20,facecolors='none',edgecolors='r',marker='^',label="WA fs_2^(2)")
    axes[2,0].plot(x,Ms0_result,color='b')
    axes[2,0].plot(x,Ms1_result,color='g')
    axes[2,0].plot(x,Ms2_result,color='r')
    axes[2,0].plot(x,m_s,color='black')
    axes[2,0].legend()
    axes[2,0].axhline(y=0,linestyle='dashed',color='gray')
    axes[2,0].set_title("M-shell s-electrons", y=1.0, pad=-14)
    
    ours1, = axes[1,1].plot(x,p0_result,color='b')
    ours2, = axes[1,1].plot(x,p2_0_result,color='g')
    ours3, = axes[1,1].plot(x,p2_1_result,color='r')
    ours11, = axes[1,1].plot(x,l_p,color='black')
    legend1 = axes[1,1].legend(
                    handles=[ours1,
                             ours2,ours3,
                             ours11],
                    labels=["0/0",
                            "0/2 + 2/0","1/1",
                            "sum"],
                    loc = 'upper left',
                    bbox_to_anchor=(0.5,1.5))
    
    axes[1,1].plot(x,p0_result,  color='b')
    axes[1,1].plot(x,p2_0_result,color='g')
    axes[1,1].plot(x,p2_1_result,color='r')
    axes[1,1].scatter(x,m_result,s=40,facecolors='none',edgecolors='b',marker='o',label="EM fp_1^(0)")
    axes[1,1].scatter(x,n_result,s=40,facecolors='none',edgecolors='g',marker='o',label="EM fp_1^(2)")
    axes[1,1].scatter(x,o_result,s=40,facecolors='none',edgecolors='r',marker='o',label="EM fp_2^(2)")
    axes[1,1].plot(x,l_p,color='black')
    axes[1,1].scatter(x,p0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")
    axes[1,1].scatter(x,p2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")
    axes[1,1].scatter(x,p2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")
    axes[1,1].scatter(x,l_p_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")
    axes[1,1].legend()
    axes[1,1].add_artist(legend1)
    axes[1,1].axhline(y=0,linestyle='dashed',color='gray')
    axes[1,1].set_title("L-shell p-electrons", y=1.0, pad=-14)
    
    axes[2,1].plot(x,M_p0_result,  color='b')
    axes[2,1].plot(x,M_p2_0_result,color='g')
    axes[2,1].plot(x,M_p2_1_result,color='r')
    axes[2,1].plot(x,m_p,color='black')
    axes[2,1].scatter(x,M_p0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")
    axes[2,1].scatter(x,M_p2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")
    axes[2,1].scatter(x,M_p2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")
    axes[2,1].scatter(x,WA_p0_result,s=40,facecolors='none',edgecolor='b',marker='o',label="Wagenfeld p-els D")
    axes[2,1].scatter(x,m_p_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")
    axes[2,1].legend()
    axes[2,1].axhline(y=0,linestyle='dashed',color='gray')
    axes[2,1].set_title("M-shell p-electrons", y=1.0, pad=-14)

    axes[0,1].set_axis_off()
    axes[0,2].set_axis_off()
    axes[1,2].set_axis_off()

    axes[2,2].plot(x,M_d0_result,  color='b')
    axes[2,2].plot(x,M_d2_0_result,color='g')
    axes[2,2].plot(x,M_d2_1_result,color='r')
    axes[2,2].plot(x,m_d,color='black')
    axes[2,2].scatter(x,M_d0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")
    axes[2,2].scatter(x,M_d2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")
    axes[2,2].scatter(x,M_d2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")
    axes[2,2].scatter(x,WA_d0_result,s=40,facecolors='none',edgecolor='b',marker='o',label="Wagenfeld d-els D")
    axes[2,2].scatter(x,m_d_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")
    axes[2,2].legend()
    axes[2,2].axhline(y=0,linestyle='dashed',color='gray')
    axes[2,2].set_title("M-shell d-electrons", y=1.0, pad=-14)
  
    plt.subplots_adjust(left=0.04, bottom=0.04, right=1.0, top=0.95, wspace=0.15, hspace=0.1)
    fig.suptitle("ORTHO alpha = {:4.2f}, theta_0 = {:4.2f}".format(alp,t0))
    plt.show()
    exit()
    
  if form_factor_test:
    print("Performing Disp Correction Calculation")
    #axis x is in nu
    x = np.logspace(math.log10(speed_of_light*1E10/lam_min), 
                    math.log10(speed_of_light*1E10/lam_max), 
                    steps)

    from time import time
    numpy.random.seed(int(np.ceil(time())))
    #t0 = numpy.random.random(1)[0] * math.pi
    t0 = 0
    TS = (1+pow(np.cos(t0),2))/2
    
    factor_x = 2*math.pi

    from brennan import brennan
    br = brennan()

    #K-shell
    k_s    = np.zeros_like(x)
    
    #L-shell
    l_s    = np.zeros_like(x)
    l_p1   = np.zeros_like(x)
    l_p3   = np.zeros_like(x)
    
    #M-shell
    m_s    = np.zeros_like(x)
    m_p1   = np.zeros_like(x)
    m_p3   = np.zeros_like(x)
    m_d3   = np.zeros_like(x)
    m_d5   = np.zeros_like(x)

    zero_angle = np.zeros_like(x)
    hönl_K  = np.zeros_like(x)
    hönl_L  = np.zeros_like(x)
    hönl_M  = np.zeros_like(x)
    hönl    = np.zeros_like(x)
    brennan_fdp = np.zeros_like(x)
    e = elements[el_nr]
    
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
      res = pool.starmap(calc_stuff_with_brennan, zip(x, repeat(t0), repeat(l_max), repeat(p_max), repeat(el_nr), repeat(br), repeat(e)))
      for i,r in enumerate(res):
        k_s[i],\
        l_s[i],\
        m_s[i],\
        l_p1[i],l_p3[i],\
        m_p1[i],m_p3[i],\
        m_d3[i],m_d5[i], brennan_fdp[i],\
        hönl_K[i], hönl_L[i], hönl_M[i], hönl[i] = r
      
      #res = pool.starmap(calc_stuff_only_sums, zip(x, repeat(0), repeat(l_max), repeat(p_max), repeat(el_nr)))
      #for i,r in enumerate(res):
      #  zero_angle[i] = sum(np.sqrt(r))
    
    x2 = speed_of_light / x * 1E10
    achse_z   = x / (get_ionization_energy_1s(el_nr)     / h)
    achse_z1  = x / (get_ionization_energy_2s(el_nr)     / h)
    achse_z2  = x / (get_ionization_energy_2p1_2(el_nr)  / h)
    achse_z3  = x / (get_ionization_energy_2p3_2(el_nr)  / h)
    achse_zm1 = x / (get_ionization_energy_3s(el_nr)     / h)
    achse_zm2 = x / (get_ionization_energy_3p_1_2(el_nr) / h)
    achse_zm3 = x / (get_ionization_energy_3p_3_2(el_nr) / h)
    achse_zm4 = x / (get_ionization_energy_3d_3_2(el_nr) / h)
    achse_zm5 = x / (get_ionization_energy_3d_5_2(el_nr) / h)
    fdp_nach_hoenl = np.zeros_like(achse_z)
    fdp_nach_EM = np.zeros_like(achse_z)
    fdp_nach_Wa = np.zeros_like(achse_z)
    delta_K  = delta(el_nr,1,0,0)
    delta_l1 = delta(el_nr,2,0,0)
    delta_l2 = delta(el_nr,2,1,0.5)
    delta_l3 = delta(el_nr,2,1,1.5)
    delta_m1 = delta(el_nr,3,0,0)
    delta_m2 = delta(el_nr,3,1,0.5)
    delta_m3 = delta(el_nr,3,1,1.5)
    delta_m4 = delta(el_nr,3,2,1.5)
    delta_m5 = delta(el_nr,3,2,2.5)
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
      res = pool.starmap(fdp_from_z, zip(achse_z, repeat(delta_K), repeat(f_s_1_hoenl)))
      for i,r in enumerate(res):
        fdp_nach_hoenl[i] = r
      
      res = pool.starmap(fdp_from_z, zip(achse_z1, repeat(delta_l1), repeat(f_s_1_EM)))
      for i,r in enumerate(res):
        fdp_nach_EM[i] += r
      res = pool.starmap(fdp_from_z, zip(achse_z2, repeat(delta_l2), repeat(f_p_1_EM)))
      for i,r in enumerate(res):
        fdp_nach_EM[i] += r
      res = pool.starmap(fdp_from_z, zip(achse_z3, repeat(delta_l3), repeat(f_p_1_EM)))
      for i,r in enumerate(res):
        fdp_nach_EM[i] += 2*r
      
      res = pool.starmap(fdp_from_z, zip(achse_zm1, repeat(delta_m1), repeat(f_s_1_WA)))
      for i,r in enumerate(res):
        fdp_nach_Wa[i] += r
      res = pool.starmap(fdp_from_z, zip(achse_zm2, repeat(delta_m2), repeat(f_p_1_WA)))
      for i,r in enumerate(res):
        fdp_nach_Wa[i] += r
      res = pool.starmap(fdp_from_z, zip(achse_zm3, repeat(delta_m3), repeat(f_p_1_WA)))
      for i,r in enumerate(res):
        fdp_nach_Wa[i] += 2*r
      res = pool.starmap(fdp_from_z, zip(achse_zm4, repeat(delta_m4), repeat(f_d_1_WA)))
      for i,r in enumerate(res):
        fdp_nach_Wa[i] += 3*r
      res = pool.starmap(fdp_from_z, zip(achse_zm5, repeat(delta_m5), repeat(f_d_1_WA)))
      for i,r in enumerate(res):
        fdp_nach_Wa[i] += 2*r
      
    fig, axes = plt.subplots(3,3)

    def plot_stuff(ax,o,p,t,name):
      ax.scatter(x2,o,s=20,facecolors='none',edgecolors='b',marker='^',label="orthogonal")
      ax.scatter(x2,p,s=20,facecolors='none',edgecolors='g',marker='^',label="parallel")
      ax.plot(x2,t,color='r',label="total")
      ax.axhline(y=0,linestyle='dashed',color='gray')
      ax.set_title(name, y=1.0, pad=-14)

    def plot(ax,t,p1=None,p2=None,name=""):
      if p1 is not None:
        ax.scatter(x2,p1,s=20,facecolors='none',edgecolors='b',marker='^',label="subshell1")
      if p2 is not None:
        ax.scatter(x2,p2,s=20,facecolors='none',edgecolors='g',marker='^',label="subshell2")
      ax.plot(x2,t,color='r',label="total")
      ax.axhline(y=0,linestyle='dashed',color='gray')
      ax.set_title(name, y=1.0, pad=-14)
    
    #plot(axes[0,0],np.sqrt(k_s)*factor_x,name="K-shell s-electrons")
    #plot(axes[1,0],np.sqrt(l_s)*factor_x,name="L-shell s-electrons")
    #plot(axes[2,0],np.sqrt(m_s)*factor_x,name="M-shell s-electrons")

    #r1, r2 = np.sqrt(l_p1)*factor_x, np.sqrt(l_p3)*factor_x
    #plot(axes[1,1],r1+2*r2,r1,r2,"L-shell p-electrons")
    #r1, r2 = np.sqrt(m_p1)*factor_x, np.sqrt(m_p3)*factor_x
    #plot(axes[2,1],r1+2*r2,r1,r2,"M-shell p-electrons")

    #r1, r2 = np.sqrt(m_d3)*factor_x, np.sqrt(m_d5)*factor_x
    #plot(axes[2,2],3*r1+2*r2,r1,r2,"M-shell d-electrons")
    #del r1
    #del r2

    res_K = k_s*factor_x
    res_L = (l_s+l_p1+2*l_p3)*factor_x
    res_M = (m_s+m_p1+2*m_p3+3*m_d3+2*m_d5)*factor_x
    resy = res_K+res_L+res_M
    #resy = (np.sqrt(k_s))
    #resy_br = 2*resy/constant_factor
    #ausdruck = (resy)/(1-4.0*math.pi/3*resy)*constant_factor
    resy_br = resy

    axes[0,2].plot(x2,resy_br,color='k',label="all")
    axes[0,2].axhline(y=0,linestyle='dashed',color='gray')
    axes[0,2].plot(x2,brennan_fdp,color='r',label="brennan")
    axes[0,2].plot(x2,hönl,color='b',label="hönl")

    axes[0,1].scatter(x2,fdp_nach_hoenl,facecolors='none',edgecolors='g',marker='^',label="K purely paper")
    axes[0,1].scatter(x2,fdp_nach_EM,facecolors='none',edgecolors='b',marker='^',label="L purely paper")
    axes[0,1].scatter(x2,fdp_nach_Wa,facecolors='none',edgecolors='r',marker='^',label="M purely paper")
    
    axes[0,2].scatter(x2,fdp_nach_Wa+fdp_nach_EM+fdp_nach_hoenl,facecolors='none',edgecolors='k',marker='^',label="total purely paper")

    axes[1,2].plot(x2,resy_br/brennan_fdp,color='k',label="ratio all/brennan")
    axes[1,2].plot(x2,hönl/brennan_fdp,color='b',label="ratio hönl/brennan")
    axes[1,2].plot(x2,resy_br/hönl,color='g',label="ratio all/hönl")
    axes[1,2].plot(x2,(fdp_nach_Wa+fdp_nach_EM+fdp_nach_hoenl)/hönl,color='r',label="ratio hönl/hönl_integral")
    axes[1,2].axhline(y=1,linestyle='dashed',color='gray')

    

    axes[0,1].plot(x2,res_K,label="K-shell electrons")
    axes[0,1].plot(x2,res_L,label="L-shell electrons")
    axes[0,1].plot(x2,res_M,label="M-shell electrons")

    axes[0,1].plot(x2,hönl_K,linestyle='dashed',label="Hönl K-shell")
    axes[0,1].plot(x2,hönl_L,linestyle='dashed',label="Hönl L-shell")
    axes[0,1].plot(x2,hönl_M,linestyle='dashed',label="Hönl M-shell")

    #axes[0,1].axhline(y=1,linestyle='dashed',color='gray')
    #steps = 1
    #print("Starting angle calculations!")
    #for i in range(0,steps+1):
    #  theta = i*(1./steps)*math.pi
    #  temp = np.zeros_like(zero_angle)
    #  print("{}/{} steps".format(i,steps))
    #  with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
    #    res = pool.starmap(calc_stuff_only_sums, zip(x, repeat(theta), repeat(l_max), repeat(p_max), repeat(el_nr)))
    #    for j,r in enumerate(res):
    #      temp[j] = sum(np.sqrt(r))
    #    axes[0,1].plot(x2,zero_angle/temp*(1+pow(np.cos(theta),2))/2,label=r"$\theta_0 = ${:4.2f}$\pi$".format(theta/math.pi))

    #axes[0,1].plot(x2,zero_angle/resy * TS,color='k',label=r"$\theta_0 = ${:5.3f}$\pi$".format(t0/math.pi))
    #axes[0,1].set_title(r"Ratio between $\theta_0=0$ and $\theta_0$", y=1.0, pad=-14)

    ticks = np.array([0.2,0.3,0.5,0.7,1.0,1.5,2.0,3.0])

    for ax in fig.get_axes():
      ax.set_xscale('log')
      ax.set_xlim(lam_min,lam_max)
      ax.set_xticks(ticks)
      ax.get_xaxis().get_major_formatter().labelOnlyBase = False
      ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
      ax.legend()

    ##This adds the frequency in ExaHz on top of the 1,2 axes
    def tickfunc(x):
      N = speed_of_light / x * 1E-8
      return  ["%.3f" %z for z in N]
    ax2 = axes[1,2].twiny()
    ax2.set_xscale('log')
    ax2.set_xlim(lam_min,lam_max)
    ax2.set_xticks(ticks)
    ax2.set_xticklabels(tickfunc(ticks))

    
    fig.suptitle(r"$\theta_0$ = {:4.2f} $\pi$".format(t0/math.pi))
    plt.subplots_adjust(left=0.05, bottom=0.04, right=1.0, top=0.95, wspace=0.10, hspace=0.175)
    plt.show()
    exit()
  
  if numplot == True: 
    x = numpy.linspace(0.9,1.3,2*steps)
    eta_nach_hoenl = np.zeros_like(x)
    eta_nach_integral = np.zeros_like(x)
    eta_nach_unserem = np.zeros_like(x)
    nprime = np.zeros_like(x)
    K_recursive = np.zeros_like(x)
    SE = np.zeros_like(x)
    
    nu_0 = (get_ionization_energy_1s(el_nr) / h)
    Z_s_sq = pow(get_Zeff_1s(el_nr),2)
    e_ion = get_ionization_energy_1s(el_nr)
    nu_k = e_ion / h
    delta_K = 1 + alpha_sq * Z_s_sq / 4 - e_ion/(Ryd_ener * Z_s_sq) #21a in Hoenl

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
      res = pool.starmap(eta_K, zip(x, repeat(1), repeat(nu_0),repeat(delta_K)))
      for i,r in enumerate(res):
        eta_nach_hoenl[i] = r
      res = pool.starmap(sugiura_k_purely_imag, zip(repeat(el_nr), x*get_ionization_energy_1s(el_nr)))
      for i,r in enumerate(res):
        eta_nach_integral[i] = r
      res = pool.starmap(eta_K_unsers, zip(x,repeat(nu_0), repeat(el_nr), repeat(delta_K)))
      for i,r in enumerate(res):
        eta_nach_unserem[i] = r
        
      from matrix_coefficients_v2 import n_prime_from_z, K_recursive_from_z
      res = pool.starmap(n_prime_from_z, zip(x,repeat(1)))
      for i,r in enumerate(res):
        nprime[i] = r
        
      res = pool.starmap(K_recursive_from_z, zip(repeat(0),repeat(1), repeat(b(1,0,el_nr)), x,repeat(1)))
      for i,r in enumerate(res):
        K_recursive[i] = r
      res = pool.starmap(sugiura_exps, zip(x,repeat(1)))
      for i,r in enumerate(res):
        SE[i] = r
      

    fig, axes = plt.subplots(2,2)
    axes[0,0].plot(x,eta_nach_hoenl,label="without integral")
    axes[0,0].plot(x,eta_nach_unserem,label="unsers")
    axes[0,0].plot(x,eta_nach_integral,label="with integral")
    axes[0,0].legend()

    ratio1 = numpy.zeros_like(eta_nach_hoenl)
    ratio2 = numpy.zeros_like(eta_nach_hoenl)
    for i,e in enumerate(eta_nach_hoenl):
      ratio1[i] = eta_nach_hoenl[i]  /eta_nach_integral[i]
      ratio2[i] = eta_nach_unserem[i]/eta_nach_integral[i]
      if x[i] <= 1:
        ratio1[i] += 1
        ratio2[i] += 1

    axes[1,0].plot(x,ratio1,label="ratio")
    axes[1,0].plot(x,ratio2,label="ratio unsers")
    axes[1,0].legend()

    axes[0,0].axvline(x=1/(1-delta_K),linestyle='dashed',color='gray')
    axes[1,0].axvline(x=1/(1-delta_K),linestyle='dashed',color='gray')
    axes[1,1].axvline(x=1,linestyle='dashed',color='gray')
    axes[0,1].axvline(x=1,linestyle='dashed',color='gray')
    
    axes[0,1].plot(x,nprime,label="nprime(z)")
    #axes[1,1].plot(x,SE,label="sugiura(z)")
    axes[1,1].plot(x,K_recursive,label="K(z)")
    print(eta_nach_hoenl[-1]/eta_nach_integral[-1], eta_nach_unserem[-1]/eta_nach_integral[-1])
    
    plt.show()
    exit()
