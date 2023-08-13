import numpy as np
import numpy.typing as npt
import math, cmath
import scipy.special as special
from typing import Any, Union, overload, Callable
from typing_extensions import Literal
import matplotlib
import matplotlib.pyplot as plt
from matrix_coefficients_v2 import K_recursive_from_z
from constants_and_atomic_properties import *

def sugiura_exps(z: Union[float,npt.NDArray], n_0: int) -> Union[float,npt.NDArray]:
  z_1 = z-1
  if z_1 < 0:
    return 0
  sqrt_var = np.sqrt(z_1)
  part2 = np.exp(-4*n_0/sqrt_var*np.arctan(sqrt_var))
  part3 = 1-np.exp(-2*n_0*math.pi/sqrt_var)
  return part2/part3

def sugiura_exps_fix(z: Union[float,npt.NDArray], n_0: int) -> Union[float,npt.NDArray]:
  z_1 = z-1
  if(z<=1):
    sqrt_var = np.sqrt(abs(z_1))
    temp = np.arctan(sqrt_var)/sqrt_var - 2/3.0*z_1*special.hyp2f1(0.75,1,1.75,z_1*z_1)
    part2 = np.exp(-4*n_0*temp)
    part3 = 1-np.exp(-2*n_0*math.pi/sqrt_var)
  else:
    sqrt_var = np.sqrt(z_1)
    part2 = np.exp(-4*n_0/sqrt_var*np.arctan(sqrt_var))
    part3 = 1-np.exp(-2*n_0*math.pi/sqrt_var)
  return part2/part3
  #The one above is without the correction for the part when z is <1 where we continue the function using the hyp2f1
  return np.exp(-4*n_0/np.sqrt(z-1)*np.arctan(np.sqrt(z-1)))/(1-np.exp(-2*n_0*math.pi/np.sqrt(z-1)))

def exp_from_z(z: float,n_0: int) -> float:
  if z<1:
    return 0
  return math.exp(-2*n_0/math.sqrt(z-1)*math.atan(math.sqrt(z-1)))

def exp_from_z_fix(z: float,n_0: int) -> float:
  z_1 = z-1
  temp = 0.0
  sqrt_var = np.sqrt(abs(z_1))
  if(z<=1):
    temp -= 2/3.0*z_1*special.hyp2f1(0.75,1,1.75,z_1*z_1)
  temp += np.arctan(sqrt_var)/sqrt_var
  return math.exp(-2*n_0*temp)

def N_square_from_z(l: int, m: int, b_: float, n_0: int, z: float) -> complex:
  if (m > l):
    return 0
  result = (2*l+1)*math.factorial(l-m)/math.factorial(l+m) * constant_factor / math.pi * n_0 * b_ * product_n_prime_from_z(n_0,z,l)
  if m >= 1:
    result *= 2
  return result

def n_prime_from_z(z:float,n_0:int) -> complex:
  return n_0/cmath.sqrt(z-1)

def product_n_prime_from_z(n_0: int,z: float,l: int) -> complex:
  n_p = n_prime_from_z(z,n_0)
  fact = complex(1.0)
  for nu in range(1,l+1):
    fact *= n_p * n_p + nu * nu
  denom = 1-cmath.exp(-2*math.pi*n_p)
  return fact/denom

def N_square_from_z_fix(l: int, m: int, b_: float, n_0: int, z: float) -> complex:
  if (m > l):
    return 0
  result = (2*l+1)*math.factorial(l-m)/math.factorial(l+m) * constant_factor / math.pi * n_0 * b_ * product_n_prime_from_z_fix(n_0,z,l)
  if m >= 1:
    result *= 2
  return result

def n_prime_from_z_fix(z:float,n_0:int) -> complex:
  return n_0/cmath.sqrt(z-1)

def product_n_prime_from_z_fix(n_0: int,z: float,l: int) -> complex:
  n_p = n_prime_from_z_fix(z,n_0)
  fact = complex(1.0)
  for nu in range(1,l+1):
    fact *= n_p * n_p + nu * nu
  denom = 1-cmath.exp(-2*math.pi*abs(n_p))
  return fact/denom

def linear_fit(z):
  return math.exp(-4)/3.0 * (4*z-1)

b_ = 30
p = 1
l = 1
z = 0.999999999999
n_0 = 1

print(K_recursive_from_z(p,l,b_,z,n_0))

#exit()
x = np.linspace(0, 4, 2000)

old = np.zeros_like(x)
new = np.zeros_like(x)
lin = np.zeros_like(x)
for i,z in enumerate(x):
  a = N_square_from_z(l,1,b_,1,x[i])
  old[i] = math.sqrt(a.imag * a.imag + a.real * a.real)
  a = N_square_from_z_fix(l,1,b_,1,x[i])
  new[i] = math.sqrt(a.imag * a.imag + a.real * a.real)
  #new[i] = N_square_from_z_fix(l,1,b_,1,x[i])
  #lin[i] = linear_fit(x[i])
plt.plot(x,old,label="old")
plt.plot(x,new,'r:',label="new")
#plt.plot(x,lin,'g:',label="linear_fit")
plt.legend()
plt.show()
