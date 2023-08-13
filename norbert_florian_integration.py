from constants_and_atomic_properties import h,speed_of_light,prefactor
from constants_and_atomic_properties import one_minus_delta_edge, elements
from constants_and_atomic_properties import get_ionization_energy_1s,get_ionization_energy_2s,get_ionization_energy_3s,get_ionization_energy_2p1_2,get_ionization_energy_2p3_2
from constants_and_atomic_properties import get_ionization_energy_3p_1_2, get_ionization_energy_3p_3_2, get_ionization_energy_3d_5_2, get_ionization_energy_3d_3_2
from constants_and_atomic_properties import get_line_width_K, get_line_width_Ls, get_line_width_Lp_1_2, get_line_width_Lp_3_2
import scipy.integrate as integrate
import math
from hoenl_like import perform_real_integration, perform_imag_integration
from typing import overload, Union, Callable
import numpy as np
import numpy.typing as npt

from legendre_plynomials import alpha_coef, beta_bar_coef, beta_coef, alpha_bar_coef
from matrix_coefficients_v2 import f_p_el_for_p,f_a_for_p,f_d_el_for_p
from matrix_coefficients_v2 import f_p_el_for_p_gM,f_a_for_p_gM,f_d_el_for_p_gM

al00  = []
al11  = []
al22  = []
al33  = []
bbl22 = []
bbl11 = []
bbl33 = []
for l in range(11):
  al00.append(alpha_coef(l,0,0,0,0))
  al11.append(alpha_coef(l,1,1,0,0))
  al22.append(alpha_coef(l,2,2,0,0))
  al33.append(alpha_coef(l,3,3,0,0))
  bbl22.append(beta_bar_coef(l,2,2,0,0))
  bbl11.append(beta_bar_coef(l,1,1,0,0))
  bbl33.append(beta_bar_coef(l,3,3,0,0))

@overload
def apply_angle_part_s_parallel(a_l: float, theta0:float, alpha:float, l:int) -> float: ...

@overload
def apply_angle_part_s_parallel(a_l: npt.NDArray, theta0:float, alpha:float, l:int) -> npt.NDArray: ...

def apply_angle_part_s_parallel(a_l: Union[float,npt.NDArray], theta0:float, alpha:float, l:int) -> Union[float,npt.NDArray]:
  a_l *= alpha_coef(l,1,1,theta0,alpha)
  return a_l

@overload
def apply_angle_part_s_orthogonal(a_l: float, theta0:float, alpha:float, l:int) -> float: ...

@overload
def apply_angle_part_s_orthogonal(a_l: npt.NDArray, theta0:float, alpha:float, l:int) -> npt.NDArray: ...

def apply_angle_part_s_orthogonal(a_l: Union[float,npt.NDArray],theta0: float,alpha:float, l: int) -> Union[float,npt.NDArray]:
  a_l *= beta_coef(l,1,1,theta0,alpha)
  return a_l

def apply_angle_part_p_parallel(b_l: float,c_0_l: float,c_2_l: float,d_l: float,_k: int, theta0: float, alpha: float, l: int) -> "list[float]":
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
  return [b_l,c_0_l,c_2_l,d_l]

def apply_angle_part_p_orthogonal(b_l: float,c_0_l: float,c_2_l: float,d_l: float,_k: int, theta0: float, alpha: float, l: int) -> "list[float]":
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
  return [b_l,c_0_l,c_2_l,d_l]

def apply_angle_part_d_parallel(vals: "list[float]",_k: int, theta0: float, alpha:float, l: int) -> "list[float]":
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
  return vals

def apply_angle_part_d_orthogonal(vals: "list[float]",_k: int, theta0: float, alpha:float, l: int) -> "list[float]":
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
  return vals


def calc_Intensity_s_orbital(alpha_loc: float, nu_in: float, t0: float, l_max: int, p_max: int, n0: int, el_nr: int) -> float:
  z_temp = None
  if n0 == 1:
    z_temp = nu_in / (get_ionization_energy_1s(el_nr) / h)
    de = one_minus_delta_edge(el_nr,1,0,0)
  elif n0 == 2:
    z_temp = nu_in / (get_ionization_energy_2s(el_nr) / h)
    de = one_minus_delta_edge(el_nr,2,0,0)
  elif n0 == 3:
    z_temp = nu_in / (get_ionization_energy_3s(el_nr) / h)
    de = one_minus_delta_edge(el_nr,3,0,0)
  else:
    z_temp = 0
    de = 0
  if z_temp < 1: return 0
  z_temp *= de

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

def calc_Intensity_p_orbital(alpha_loc: float, nu_in: float, t0: float, l_max: int, p_max: int, n0: int, subshell: int, el_nr: int) -> float:
  if n0 == 1:
    print("This shell doesn't have p orbitals")
    exit(-1)
  elif n0 == 2:
    if subshell == 1:
      z_temp = nu_in / (get_ionization_energy_2p1_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,2,1,0.5)
    else:
      z_temp = nu_in / (get_ionization_energy_2p3_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,2,1,1.5)
  elif n0 == 3:
    if subshell == 1:
      z_temp = nu_in / (get_ionization_energy_3p_1_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,3,1,0.5)
    else:
      z_temp = nu_in / (get_ionization_energy_3p_3_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,3,1,1.5)
  else:
    z_temp = 0
    de = 0
  if z_temp < 1: return 0
  z_temp *= de

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

def calc_Intensity_d_orbital(alpha_loc: float, nu_in: float, t0: float, l_max: int, p_max: int, n0: int, subshell: int, el_nr: int) -> float:
  if n0 == 1:
    print("This shell doesn't have d orbitals")
    exit(-1)
  elif n0 == 2:
    print("This shell doesn't have d orbitals")
    exit(-1)
  elif n0 == 3:
    if subshell == 1:
      z_temp = nu_in / (get_ionization_energy_3d_3_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,3,2,1.5)
    else:
      z_temp = nu_in / (get_ionization_energy_3d_5_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,3,2,2.5)
  else:
    z_temp = 0
    de = 0
  if z_temp < 1: return 0
  z_temp *= de

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

def calc_Intensity_s_orbital_wM(alpha_loc: float, nu_in: float, t0: float, l_max: int, p_max: int, n0: int, el_nr: int, M_matrizes) -> float:
  l_M = []
  if n0 == 1:
    l_M = M_matrizes[1][1]
    z_temp = nu_in / (get_ionization_energy_1s(el_nr) / h)
    de = one_minus_delta_edge(el_nr,1,0,0)
  elif n0 == 2:
    l_M = M_matrizes[1][2]
    z_temp = nu_in / (get_ionization_energy_2s(el_nr) / h)
    de = one_minus_delta_edge(el_nr,2,0,0)
  elif n0 == 3:
    l_M = M_matrizes[1][3]
    z_temp = nu_in / (get_ionization_energy_3s(el_nr) / h)
    de = one_minus_delta_edge(el_nr,3,0,0)
  else:
    z_temp = 0
    de = 0
  if z_temp < 1: return 0
  z_temp *= de

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
      for r in f_a_for_p_gM(el_nr, l, 0, z_temp, nu_in, n0, p, l_M):
        temp += r.real
    par += alpha_coef(l,1,1,t0,alpha_loc) * temp
    orth += beta_coef(l,1,1,t0,alpha_loc) * temp
  
  return (par**2 + orth**2)

def calc_Intensity_p_orbital_wM(alpha_loc: float, nu_in: float, t0: float, l_max: int, p_max: int, n0: int, subshell: int, el_nr: int, M_matrizes) -> float:
  l_M = []
  if n0 == 1:
    print("This shell doesn't have p orbitals")
    exit(-1)
  elif n0 == 2:
    if subshell == 1:
      l_M = M_matrizes[2][2]
      z_temp = nu_in / (get_ionization_energy_2p1_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,2,1,0.5)
    else:
      l_M = M_matrizes[3][2]
      z_temp = nu_in / (get_ionization_energy_2p3_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,2,1,1.5)
  elif n0 == 3:
    if subshell == 1:
      l_M = M_matrizes[2][3]
      z_temp = nu_in / (get_ionization_energy_3p_1_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,3,1,0.5)
    else:
      l_M = M_matrizes[3][3]
      z_temp = nu_in / (get_ionization_energy_3p_3_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,3,1,1.5)
  else:
    z_temp = 0
    de = 0
  if z_temp < 1: return 0
  z_temp *= de

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
        for num,r in enumerate(f_p_el_for_p_gM(el_nr, l, ms[nummy], ms[nummy], z_temp, nu_in, n0, p, l_M)):
          par += apply_angle_part_s_parallel(r.real * mults[nummy], t0, alpha_loc, fac[num])
          orth += apply_angle_part_s_orthogonal(r.real * mults[nummy], t0, alpha_loc, fac[num])
  
  #Has division by 3**2 to avergae over all 3 p orbitals a.k.a. the p-wonder
  return (par**2 + orth**2) / 9

def calc_Intensity_d_orbital_wM(alpha_loc: float, nu_in: float, t0: float, l_max: int, p_max: int, n0: int, subshell: int, el_nr: int, M_matrizes) -> float:
  l_M = []
  if n0 == 1:
    print("This shell doesn't have d orbitals")
    exit(-1)
  elif n0 == 2:
    print("This shell doesn't have d orbitals")
    exit(-1)
  elif n0 == 3:
    if subshell == 1:
      l_M = M_matrizes[4][3]
      z_temp = nu_in / (get_ionization_energy_3d_3_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,3,2,1.5)
    else:
      l_M = M_matrizes[5][3]
      z_temp = nu_in / (get_ionization_energy_3d_5_2(el_nr) / h)
      de = one_minus_delta_edge(el_nr,3,2,2.5)
  else:
    z_temp = 0
    de = 0
  if z_temp < 1: return 0
  z_temp *= de

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
        for num,r in enumerate(f_d_el_for_p_gM(el_nr, l, ms[nummy], ms[nummy], z_temp, nu_in, n0, p, l_M)):
          
          par += apply_angle_part_s_parallel(r.real * mults[nummy], t0, alpha_loc, fac[num])
          orth += apply_angle_part_s_orthogonal(r.real * mults[nummy], t0, alpha_loc, fac[num])
  #Has division by 5**2 to accomodate the averaging over 5 d orbtials  a.k.a. the d-wonder
  return (par**2 + orth**2) / 25.0

from constants_and_atomic_properties import constant_factor
def k_shell_s(nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return math.sqrt(integrate.quad(calc_Intensity_s_orbital_wM, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 1,    el_nr, M_matrizes))[0] / (2*math.pi*constant_factor**2))
def l_shell_s(nu_in, t0, l_max, p_max, el_nr, M_matrizes):   
  return math.sqrt(integrate.quad(calc_Intensity_s_orbital_wM, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 2,    el_nr, M_matrizes))[0] / (2*math.pi*constant_factor**2))
def m_shell_s(nu_in, t0, l_max, p_max, el_nr, M_matrizes):   
  return math.sqrt(integrate.quad(calc_Intensity_s_orbital_wM, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 3,    el_nr, M_matrizes))[0] / (2*math.pi*constant_factor**2))
def l_shell_p1(nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return math.sqrt(integrate.quad(calc_Intensity_p_orbital_wM, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 2, 1, el_nr, M_matrizes))[0] / (2*math.pi*constant_factor**2))
def l_shell_p3(nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return math.sqrt(integrate.quad(calc_Intensity_p_orbital_wM, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 2, 2, el_nr, M_matrizes))[0] / (2*math.pi*constant_factor**2))
def m_shell_p1(nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return math.sqrt(integrate.quad(calc_Intensity_p_orbital_wM, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 3, 1, el_nr, M_matrizes))[0] / (2*math.pi*constant_factor**2))
def m_shell_p3(nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return math.sqrt(integrate.quad(calc_Intensity_p_orbital_wM, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 3, 2, el_nr, M_matrizes))[0] / (2*math.pi*constant_factor**2))
def m_shell_d3(nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return math.sqrt(integrate.quad(calc_Intensity_d_orbital_wM, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 3, 1, el_nr, M_matrizes))[0] / (2*math.pi*constant_factor**2))
def m_shell_d5(nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return math.sqrt(integrate.quad(calc_Intensity_d_orbital_wM, 0, 2*math.pi, args=(nu_in, t0, l_max, p_max, 3, 2, el_nr, M_matrizes))[0] / (2*math.pi*constant_factor**2))

def integrand_for_k_s(nu_j, nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return k_shell_s(nu_in, t0, l_max, p_max, el_nr, M_matrizes) * 2/nu_j
def integrand_for_l_s(nu_j, nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return l_shell_s(nu_in, t0, l_max, p_max, el_nr, M_matrizes) * 2/nu_j
def integrand_for_l_p1(nu_j, nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return l_shell_p1(nu_in, t0, l_max, p_max, el_nr, M_matrizes) * 2/nu_j
def integrand_for_l_p3(nu_j, nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return l_shell_p3(nu_in, t0, l_max, p_max, el_nr, M_matrizes) * 2/nu_j
def integrand_for_m_s(nu_j, nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return m_shell_s(nu_in, t0, l_max, p_max, el_nr, M_matrizes) * 2/nu_j
def integrand_for_m_p1(nu_j, nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return m_shell_p1(nu_in, t0, l_max, p_max, el_nr, M_matrizes) * 2/nu_j
def integrand_for_m_p3(nu_j, nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return m_shell_p3(nu_in, t0, l_max, p_max, el_nr, M_matrizes) * 2/nu_j
def integrand_for_m_d3(nu_j, nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return m_shell_d3(nu_in, t0, l_max, p_max, el_nr, M_matrizes) * 2/nu_j
def integrand_for_m_d5(nu_j, nu_in, t0, l_max, p_max, el_nr, M_matrizes):
  return m_shell_d5(nu_in, t0, l_max, p_max, el_nr, M_matrizes) * 2/nu_j

def real_NF_integration(n_j:float,chi_j:float,nu:float,t0:float,l_max:int,p_max:int,el_nr:int,oscillator_density_function:Callable, M_matrizes) -> float:
  res = (pow(n_j,2)-pow(nu,2))/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*oscillator_density_function(n_j, nu, t0, l_max, p_max, el_nr, M_matrizes)

def imag_NF_integration(n_j:float,chi_j:float,nu:float,t0:float,l_max:int,p_max:int,el_nr:int,oscillator_density_function:Callable, M_matrizes) -> float:
  res = -nu*chi_j/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*oscillator_density_function(n_j, nu, t0, l_max, p_max, el_nr, M_matrizes)

def perform_i_integration(nu: float, nu_edge: float, x_j:float, t0:float,l_max:int,p_max:int,el_nr:int, func: Callable, M_matrizes) -> float:
  epsilon2 = 0.1
  lower_limit = nu_edge
  if (1+epsilon2)*nu < nu_edge:
    upper_limit = 1000*nu_edge
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral = integrate.quad(imag_NF_integration,
                                 lower_limit,
                                 upper_limit,
                                 args=(x_j,nu,t0,l_max,p_max,el_nr,func, M_matrizes),
                                 points=nu,
                                 limit=200000,
                                 epsabs=1E-55,
                                 epsrel=1E-10)
  if imag_integral[1] > 0.5*abs(imag_integral[0]):
    print("!!! inaccurate IMAG integral at nu_edge/nu = " + "{:6.4f}".format(nu_edge/nu) + " ratio:" + str(imag_integral[1]/abs(imag_integral[0])) + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))
  return imag_integral[0]

def perform_r_integration(nu: float, nu_edge:float, x_j:float, t0:float,l_max:int,p_max:int,el_nr:int,func: Callable, M_matrizes) -> float:
  inte1 = integrate.quad(real_NF_integration,
                         nu_edge,
                         1000*nu,
                         args=(x_j,nu,t0,l_max,p_max,el_nr,func, M_matrizes),
                         points=nu,
                         limit=200000,
                         epsabs=1E-50,
                         epsrel=1E-10)
  if inte1[1] > 0.5*abs(inte1[0]):
    print("!!! inaccurate REAL1 integral at nu_edge/nu = " + "{:8.4f}".format(nu_edge/nu) + " " + str(inte1) + str(inte1[1]/abs(inte1[0])))
  return inte1[0]

def NF_K(Z: int = -1, nu:float =8047.8, t0:float = 0.0, l_max:int = 8, p_max:int = 4, disp_only: bool=False, M_matrizes=[]):
  nu_k = get_ionization_energy_1s(Z) / h
  #omdelta_K = one_minus_delta_edge(Z,1,0,0) #21a in Hoenl
  #n_0 = nu_k/omdelta_K #22 in Hoenl
  x_j = get_line_width_K(Z)/h
  
  if disp_only == True:
    n_disp_el = integrate.quad_vec(integrand_for_k_s,nu_k,20000*nu_k,args=(nu,t0,l_max,p_max,Z),limit=200000,epsabs=1E-20,workers=16)[0]

    return 2*n_disp_el

  #REAL PART
  integral = perform_r_integration(nu,nu_k, x_j, t0,l_max,p_max,Z, integrand_for_k_s, M_matrizes)
  #IMAG PART
  imag_integral = perform_i_integration(nu,nu_k, x_j, t0,l_max,p_max,Z, integrand_for_k_s, M_matrizes)

  alpha_K_sugiura_damp = complex(integral,imag_integral)
  ausdruck = alpha_K_sugiura_damp/(1-4.0*math.pi/3*alpha_K_sugiura_damp*prefactor)
  factor = pow(nu,2)
  return complex(-ausdruck.real*factor, -ausdruck.imag*factor)

def NF_L(Z: int = -1, nu:float =8047.8, t0:float = 0.0, l_max:int = 8, p_max:int = 4, disp_only: bool=False, M_matrizes=[]):
  nu_l1 = get_ionization_energy_2s(Z)   /h
  nu_l2 = get_ionization_energy_2p1_2(Z)/h
  nu_l3 = get_ionization_energy_2p3_2(Z)/h
  #omdelta_l1 = one_minus_delta_edge(Z,2,0,0)
  #omdelta_l2 = one_minus_delta_edge(Z,2,1,0.5)
  #omdelta_l3 = one_minus_delta_edge(Z,2,1,1.5)
  x_j1 = get_line_width_Ls(Z)/h
  x_j2 = get_line_width_Lp_1_2(Z)/h
  x_j3 = get_line_width_Lp_3_2(Z)/h
  #n_0_1 = nu_l1/omdelta_l1 #22 in Hoenl
  #n_0_2 = nu_l2/omdelta_l2 #22 in Hoenl
  #n_0_3 = nu_l3/omdelta_l3 #22 in Hoenl
  
  if disp_only == True:
    n_disp_2s    = integrate.quad_vec(integrand_for_l_s ,nu_l1, 20000*nu_l1,args=(nu,t0,8,4,Z),limit=200000,epsabs=1E-20,workers=16)[0]
    print("L-s done", n_disp_2s)
    n_disp_2p1_2 = integrate.quad_vec(integrand_for_l_p1,nu_l2, 20000*nu_l2,args=(nu,t0,8,4,Z),limit=200000,epsabs=1E-20,workers=16)[0]
    print("L-p1 done", n_disp_2p1_2)
    n_disp_2p3_2 = integrate.quad_vec(integrand_for_l_p3,nu_l3, 20000*nu_l3,args=(nu,t0,8,4,Z),limit=200000,epsabs=1E-20,workers=16)[0]
    print("L-p3 done", n_disp_2p3_2)

    return 2*n_disp_2s+2*n_disp_2p1_2+4*n_disp_2p3_2

  #REAL PART
  #When telling hte quad alorithm that nu is a "position to be carefull with" the results don't require epsilon
  integral_1 = perform_r_integration(nu,nu_l1, x_j1, t0,l_max,p_max,Z, integrand_for_l_s , M_matrizes)
  integral_2 = perform_r_integration(nu,nu_l2, x_j2, t0,l_max,p_max,Z, integrand_for_l_p1, M_matrizes)
  integral_3 = perform_r_integration(nu,nu_l3, x_j3, t0,l_max,p_max,Z, integrand_for_l_p3, M_matrizes)

  #IMAG PART
  imag_integral_1 = perform_i_integration(nu,nu_l1, x_j1, t0,l_max,p_max,Z, integrand_for_l_s , M_matrizes)
  imag_integral_2 = perform_i_integration(nu,nu_l2, x_j2, t0,l_max,p_max,Z, integrand_for_l_p1, M_matrizes)
  imag_integral_3 = perform_i_integration(nu,nu_l3, x_j3, t0,l_max,p_max,Z, integrand_for_l_p3, M_matrizes)
  
  alpha_L_sugiura_damp = complex(integral_1+integral_2+integral_3*2,imag_integral_1+imag_integral_2+2*imag_integral_3)
  ausdruck = (alpha_L_sugiura_damp)/(1-4.0*math.pi/3*alpha_L_sugiura_damp*prefactor)
  factor = pow(nu,2)
  return complex(-ausdruck.real*factor, -ausdruck.imag*factor)

def NF_M(Z: int=-1, nu:float=8047.8, t0:float = 0.0, l_max:int = 8, p_max:int = 4, disp_only:bool =False, M_matrizes=[]):
  nu_m1 = get_ionization_energy_3s(Z)    / h
  nu_m2 = get_ionization_energy_3p_1_2(Z)/ h
  nu_m3 = get_ionization_energy_3p_3_2(Z)/ h
  nu_m4 = get_ionization_energy_3d_3_2(Z)/ h
  nu_m5 = get_ionization_energy_3d_5_2(Z)/ h
  #omdelta_m1 = one_minus_delta_edge(Z,3,0,0)
  #omdelta_m2 = one_minus_delta_edge(Z,3,1,0.5)
  #omdelta_m3 = one_minus_delta_edge(Z,3,1,1.5)
  #omdelta_m4 = one_minus_delta_edge(Z,3,2,1.5)
  #omdelta_m5 = one_minus_delta_edge(Z,3,2,2.5)
  #n_0_1 = nu_m1/omdelta_m1 #22 in Hoenl
  #n_0_2 = nu_m2/omdelta_m2 #22 in Hoenl
  #n_0_3 = nu_m3/omdelta_m3 #22 in Hoenl
  #n_0_4 = nu_m4/omdelta_m4 #22 in Hoenl
  #n_0_5 = nu_m5/omdelta_m5 #22 in Hoenl

  if disp_only == True:
    n_disp_3s    = integrate.quad_vec(integrand_for_m_s ,nu_m1,10000*nu_m1,args=(nu,t0,8,4,Z),limit=200000,epsabs=1E-60,workers=16)[0]
    n_disp_3p1_2 = integrate.quad_vec(integrand_for_m_p1,nu_m2,10000*nu_m2,args=(nu,t0,8,4,Z),limit=200000,epsabs=1E-60,workers=16)[0]
    n_disp_3p3_2 = integrate.quad_vec(integrand_for_m_p3,nu_m3,10000*nu_m3,args=(nu,t0,8,4,Z),limit=200000,epsabs=1E-60,workers=16)[0]
    n_disp_3d3_2 = integrate.quad_vec(integrand_for_m_d3,nu_m4,10000*nu_m4,args=(nu,t0,8,4,Z),limit=200000,epsabs=1E-60,workers=16)[0]
    n_disp_3d5_2 = integrate.quad_vec(integrand_for_m_d5,nu_m5,10000*nu_m5,args=(nu,t0,8,4,Z),limit=200000,epsabs=1E-60,workers=16)[0]

    return 2*n_disp_3s+2*n_disp_3p1_2+4*n_disp_3p3_2+4*n_disp_3d3_2+6*n_disp_3d5_2

  x_j1 = get_line_width_Ls(Z)*(get_ionization_energy_3s(Z)/get_ionization_energy_2s(Z))/h
  x_j2 = get_line_width_Lp_1_2(Z)*(get_ionization_energy_3p_1_2(Z)/get_ionization_energy_2p1_2(Z))/h
  x_j3 = get_line_width_Lp_3_2(Z)*(get_ionization_energy_3p_3_2(Z)/get_ionization_energy_2p3_2(Z))/h
  x_j4 = get_line_width_Lp_3_2(Z)*(get_ionization_energy_3d_3_2(Z)/get_ionization_energy_2p3_2(Z))/h
  x_j5 = x_j4
  
  #REAL PART
  integral_1 = perform_r_integration(nu,nu_m1, x_j1, t0,l_max,p_max,Z, integrand_for_m_s , M_matrizes)
  integral_2 = perform_r_integration(nu,nu_m2, x_j2, t0,l_max,p_max,Z, integrand_for_m_p1, M_matrizes)
  integral_3 = perform_r_integration(nu,nu_m3, x_j3, t0,l_max,p_max,Z, integrand_for_m_p3, M_matrizes)
  integral_4 = perform_r_integration(nu,nu_m4, x_j4, t0,l_max,p_max,Z, integrand_for_m_d3, M_matrizes)
  integral_5 = perform_r_integration(nu,nu_m5, x_j5, t0,l_max,p_max,Z, integrand_for_m_d5, M_matrizes)

  #IMAG PART
  imag_integral_1 = perform_i_integration(nu,nu_m1, x_j1, t0,l_max,p_max,Z, integrand_for_m_s , M_matrizes)
  imag_integral_2 = perform_i_integration(nu,nu_m2, x_j2, t0,l_max,p_max,Z, integrand_for_m_p1, M_matrizes)
  imag_integral_3 = perform_i_integration(nu,nu_m3, x_j3, t0,l_max,p_max,Z, integrand_for_m_p3, M_matrizes)
  imag_integral_4 = perform_i_integration(nu,nu_m4, x_j4, t0,l_max,p_max,Z, integrand_for_m_d3, M_matrizes)
  imag_integral_5 = perform_i_integration(nu,nu_m5, x_j5, t0,l_max,p_max,Z, integrand_for_m_d5, M_matrizes)
  
  alpha_M_sugiura_damp = complex(integral_1     +integral_2+     2*integral_3+     3*integral_4+     2*integral_5,\
                                 imag_integral_1+imag_integral_2+2*imag_integral_3+3*imag_integral_4+2*imag_integral_5)
  ausdruck = (alpha_M_sugiura_damp)/(1-4.0*math.pi/3*alpha_M_sugiura_damp*prefactor)
  factor = pow(nu,2)
  return complex(-ausdruck.real*factor, -ausdruck.imag*factor)

def disp_el_int(func, nu_edge, nu,t0,l_max,p_max,Z, M_matrizes):
  t = integrate.quad(func ,nu_edge ,20000*nu_edge ,args=(nu,t0,l_max,p_max,Z, M_matrizes),limit=200000,epsabs=1E-20)[0]
  print(func," yields ",t)
  return t

def get_z_temp(el_nr, n0, orbital_type, nu_in):
  z_temp = -100
  if orbital_type == 1:
    if n0 == 1:
      z_temp = nu_in / (get_ionization_energy_1s(el_nr) / h)
      de = one_minus_delta_edge(el_nr,1,0,0)
    elif n0 == 2:
      z_temp = nu_in / (get_ionization_energy_2s(el_nr) / h)
      de = one_minus_delta_edge(el_nr,2,0,0)
    elif n0 == 3:
      z_temp = nu_in / (get_ionization_energy_3s(el_nr) / h)
      de = one_minus_delta_edge(el_nr,3,0,0)
    else:
      z_temp = 0
      de = 0
    if z_temp < 1: return 0
    z_temp *= de
  elif orbital_type == 2 or orbital_type == 3:
    if n0 == 1:
      return -100
    elif n0 == 2:
      if orbital_type == 2:
        z_temp = nu_in / (get_ionization_energy_2p1_2(el_nr) / h)
        de = one_minus_delta_edge(el_nr,2,1,0.5)
      else:
        z_temp = nu_in / (get_ionization_energy_2p3_2(el_nr) / h)
        de = one_minus_delta_edge(el_nr,2,1,1.5)
    elif n0 == 3:
      if orbital_type == 2:
        z_temp = nu_in / (get_ionization_energy_3p_1_2(el_nr) / h)
        de = one_minus_delta_edge(el_nr,3,1,0.5)
      else:
        z_temp = nu_in / (get_ionization_energy_3p_3_2(el_nr) / h)
        de = one_minus_delta_edge(el_nr,3,1,1.5)
    else:
      z_temp = 0
      de = 0
    if z_temp < 1: return 0
    z_temp *= de
  elif orbital_type == 4 or orbital_type == 5:
    if n0 == 1:
      return -100
    elif n0 == 2:
      return -100
    elif n0 == 3:
      if orbital_type == 4:
        z_temp = nu_in / (get_ionization_energy_3d_3_2(el_nr) / h)
        de = one_minus_delta_edge(el_nr,3,2,1.5)
      else:
        z_temp = nu_in / (get_ionization_energy_3d_5_2(el_nr) / h)
        de = one_minus_delta_edge(el_nr,3,2,2.5)
    else:
      z_temp = 0
      de = 0
    if z_temp < 1: return 0
    z_temp *= de
  return z_temp

def NF_DL(Z: int=-1, nu:float=8047.8, t0:float = 0.0, l_max:int = 5, p_max:int = 4, M_matrizes=[]):
  import multiprocessing
  from itertools import repeat

  funcs = [
    integrand_for_k_s ,
    integrand_for_l_s ,
    integrand_for_l_p1,
    integrand_for_l_p3,
    integrand_for_m_s ,
    integrand_for_m_p1,
    integrand_for_m_p3,
    integrand_for_m_d3,
    integrand_for_m_d5
  ]
  nus = [
    get_ionization_energy_3s(Z)    / h,
    get_ionization_energy_3p_1_2(Z)/ h,
    get_ionization_energy_3p_3_2(Z)/ h,
    get_ionization_energy_3d_3_2(Z)/ h,
    get_ionization_energy_3d_5_2(Z)/ h,
    get_ionization_energy_2s(Z)    / h,
    get_ionization_energy_2p1_2(Z) / h,
    get_ionization_energy_2p3_2(Z) / h,
    get_ionization_energy_1s(Z)    / h
  ]
  
  with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
    res = pool.starmap(disp_el_int, zip(funcs,nus, repeat(nu), repeat(t0), repeat(l_max), repeat(p_max), repeat(Z), repeat(M_matrizes)))
    return 2*res[0] + 2*res[1]+2*res[2]+4*res[3] + 2*res[4]+2*res[5]+4*res[6]+4*res[7]+6*res[8]

def NF(Z=38,nu=1E18,t0=0.0,l_max=8,p_max=4,disp_only=False,M_matrizes=[]):
    # this one will not correct for the Disp electron count!
    if disp_only == True:
      #K = NF_K(Z=Z,nu=nu,t0=t0,disp_only=disp_only)
      #print("K done!",K)
      #L = NF_L(Z=Z,nu=nu,t0=t0,disp_only=disp_only)
      #print("L done!",L)
      #M = NF_M(Z=Z,nu=nu,t0=t0,disp_only=disp_only)
      #print("M done!",M)
      #exit()
      #return K+L+M
      return NF_DL(Z,nu,t0,l_max,p_max,M_matrizes)
    else:
      K = NF_K(Z=Z,nu=nu,t0=t0,l_max=l_max,p_max=p_max,disp_only=disp_only,M_matrizes=M_matrizes)
      print("K:",K)
      L = NF_L(Z=Z,nu=nu,t0=t0,l_max=l_max,p_max=p_max,disp_only=disp_only,M_matrizes=M_matrizes)
      print("L:",L)
      M = NF_M(Z=Z,nu=nu,t0=t0,l_max=l_max,p_max=p_max,disp_only=disp_only,M_matrizes=M_matrizes)
      print("M:",M)
      return K+L+M