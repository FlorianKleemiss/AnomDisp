from matplotlib import ticker
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np
import numpy.typing as npt
import multiprocessing
from itertools import repeat
import tkinter as tk
from tkinter import ttk
from brennan import brennan
from time import time
import math, cmath
from typing import overload, Union
import pyximport
pyximport.install(language_level=3)


from legendre_plynomials import alpha_coef, beta_bar_coef, beta_coef, alpha_bar_coef
from matrix_coefficients_v2 import f_p_el_for_p,f_a_for_p,f_d_el_for_p
#from matrix_coefficients_v2 import f_a,f_c_0, f_c_2, integrand_matrix_p, integrand_matrix_s
from constants_and_atomic_properties import *
from constants_and_atomic_properties import speed_of_light, constant_factor, h, one_minus_delta_edge, a0, xn, b, elements, sugiura_exps, n_prime_from_z
from hoenl_like import calc_hoenllike, sugiura_k_purely_imag
from norbert_florian_integration import *
#from mpl_toolkits import mplot3d

Sr = [17.566,1.556,9.818,14.099,5.422,0.166,2.669,132.376,2.506,-1.586,3.325]
atom_x = 0.
atom_y = 0.
atom_z = 0.
atom_u = 0.04
nu = speed_of_light*1E10/0.71078

def DW(stol_sq, uiso) -> float:
  biso = uiso * 8*math.pi*math.pi
  return math.exp(-biso*stol_sq)
def stol_from_hkl(h,k,l) -> float: #This is hardcoded for the 10^3 Angs 90-90-90 cell
  lower = (pow(h, 2) + pow(k, 2) + pow(l, 2)) / 100.0 
  res = math.sqrt(1 / lower)
  return 1.0 / (2 * res)

def F_calc(indices,ndisp,precalc_array) -> "list[float]":
  h,k,l = indices
  if h==0 and k==0 and l == 0:
    two_theta = 0
    f_0 = 0
    for i in range(5):
      f_0 += Sr[2*i]
    fp_fdp = precalc_array[1][0]
    fp_fdp_0 = precalc_array[1][0]
    f  = (f_0 + fp_fdp)
    f0 = (f_0 + fp_fdp_0)
    phase = cmath.polar(f)[1]/math.pi*180
    print("Done with",h,k,l, "t0 ", two_theta, " fp_fdp ",fp_fdp, " f ",f, "f0: ",f0)
    return [h,k,l,abs(f*f), 0.01, abs(f0), phase]
  stol = stol_from_hkl(h,k,l)
  stol_sq = pow(stol,2)
  two_theta = math.asin(stol*0.71078)*2.0
  TS = np.sqrt((1.0+np.cos(two_theta)**2)/2.0)
  f_0 = complex(0)
  for i in range(4):
    f_0 += Sr[2*i] * math.exp(-Sr[2*i+1]*stol_sq)
    
  f_0 += Sr[8]
  fp_fdp = np.interp(two_theta,precalc_array[0],precalc_array[1])/TS
  #phi = 2*math.pi * (atom_x*h/10.0+atom_y*k/10.0+atom_z*10.0)
  phi = cmath.rect(1.0,0.0)
  T = DW(stol_sq,atom_u)
  #f = (f_0 + fp_fdp - ndisp) * phi * T
  f = (f_0 + fp_fdp) * phi * T
  phase = cmath.polar(f)[1]/math.pi*180
  fp_fdp_0 = precalc_array[1][0]
  f0 = (f_0 + fp_fdp_0) * T * phi
  print("Done with",h,k,l, "t0 ", two_theta, " fp_fdp ",fp_fdp, " f ",f, "f0: ",f0)
  return [h, k, l, abs(f*f), abs(f*f)*0.01, abs(f0), phase]

def inches(cm):
  return cm / 2.54

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

def calc_stuff_only_sums(nu_in, t0, l_max, p_max, el_nr) -> float:
  t1,t2,t3,t4,t5,t6,t7,t8,t9 = calc_stuff(nu_in, t0, l_max, p_max, el_nr)

  return sum([t1,t2,t3,t4+2*t5,t6+2*t7,3*t8+2*t9])

def calc_stuff_only_K(nu_in, t0, l_max, p_max, el_nr) -> float:
  t1,t2,t3,t4,t5,t6,t7,t8,t9 = calc_stuff(nu_in, t0, l_max, p_max, el_nr)

  return t1

def fdp_from_z(z:float, eins_minus_delta:float, osz_func) -> float:
  if z<1:  return 0
  return osz_func(eins_minus_delta * z)*2*math.pi
  #return 2**7/3.0 *sugiura_exps((1 - delta) * z,n0)/((1 - delta) * z)**3*math.pi

def eta_K(z: float,n0:int,nu_0:float,eins_minus_delta_K:float) -> float:
  if z<1: return 0
  z_eff = eins_minus_delta_K * z
  nu = z * nu_0
  #return 2**7/3/z_eff**4*sugiura_exps(z_eff,n0)/nu_0/nu*math.pi/2*emd
  return f_s_1_hoenl(z_eff)/nu**2*math.pi

def eta_K_unsers(z: float, nu_0:float, el_nr:int) -> float:
  nu_in = z * nu_0
  f = calc_stuff_only_K(nu_in,0,7,1,el_nr)
  #return 2**7/3/z_eff**4*sugiura_exps(z_eff,n0)/nu_0/nu*math.pi/2*emd
  return f/nu_in**2*math.pi

def eta_all_unser(nu_in, t0, l_max, p_max, el_nr):
  f = calc_stuff_only_sums(nu_in, t0, l_max, p_max, el_nr)
  return f/nu_in**2*math.pi

if __name__ == "__main__":
  test_1s = False
  higher_p_test = False
  M_shell_test = False
  M_shell_test_orthogonal = False
  form_factor_test = False
  numplot = False
  theta_plot = False
  fcf_test = False
  cmaps = ['blue','orange','green','red','black','blue','orange','green','red','black']
  t0:float = 0
  alp:float = 0
  l_max:int = 5
  p_max:int = 1
  el_nr:int = 52
  lam_min:float = 3.00
  lam_max:float = 0.15
  steps: int = 200
  theta_steps:int = 3
  overwrite_theta: bool = False
  x = None

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
  v8 = tk.StringVar()
  v8.set(str(theta_steps))
  v9 = tk.IntVar()
  options = [
        ('1s-test', 0),
        ('p test', 2),
        ('M miracle', 3),
        ('M miracle orthogonal', 5),
        ('Dispersion Coefficient Plot', 4),
        ('Numerical Plot', 6),
        ('Theta Plot', 7),
        ('fcf&cif of angular_effect', 8)
        ]
  result = 0
  def set_result():
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
  tk.Label(root,text="Choose number of theta steps",state="active",justify=tk.CENTER,padx=20).pack()
  tk.Entry(root, textvariable=v8).pack()
  tk.Checkbutton(root,text = "Choose theta0 = 0?", variable=v9).pack()
  tk.Label(root,text="Choose what to run",state="active",justify=tk.CENTER,padx=20).pack()
  for l,val in options:
    tk.Radiobutton(root,text=l,padx=20,indicatoron=False, width=30, height=2, variable=v,command=set_result,value=val).pack()
  root.mainloop()
  result = int(v.get())
  el_nr = int(v2.get())
  l_max = int(v3.get())
  p_max = int(v4.get())
  lam_max = float(v5.get())
  lam_min = float(v6.get())
  steps = int(v7.get())
  theta_steps = int(v8.get())
  overwrite_theta = bool(v9.get())
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
  elif result == 7:
    theta_plot = True
  elif result == 8:
    fcf_test = True

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
    np.random.seed(int(np.ceil(time())))
    t0 = np.random.random(1)[0] * math.pi
    alp = np.random.random(1)[0] * math.pi * 2
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
      s_0_result = np.zeros_like(x)
      s_1_1_result = np.zeros_like(x)#
      s_1_2_result = np.zeros_like(x)
      s_2_1_result = np.zeros_like(x)
      s_2_2_result = np.zeros_like(x)
      s_2_3_result = np.zeros_like(x)
    
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
        s_0_result[i] += apply_angle_part_s_parallel(res_0.real, 0, 0, l)
        res_0 = apply_angle_part_s_parallel(res_0.real, t0, alp, l)
        res_2 = f_a_for_p(el_nr, l, k_, z, nu_in, 1, 2)
        s_1_1_result[i] += apply_angle_part_s_parallel(res_2[0].real, 0, 0, l)
        s_1_1_result[i] += apply_angle_part_s_parallel(res_2[-1].real, 0, 0, l)
        s_1_2_result[i] += apply_angle_part_s_parallel(res_2[1].real, 0, 0, l)
        for runny in range(len(res_2)):
          res_2[runny] = apply_angle_part_s_parallel(res_2[runny].real,t0,alp, l)
        res_4 = f_a_for_p(el_nr,l,k_,z,nu_in,1,4)
        s_2_1_result[i] += apply_angle_part_s_parallel(res_4[0].real, 0, 0, l)
        s_2_1_result[i] += apply_angle_part_s_parallel(res_4[-1].real, 0, 0, l)
        s_2_2_result[i] += apply_angle_part_s_parallel(res_4[1].real, 0, 0, l)
        s_2_2_result[i] += apply_angle_part_s_parallel(res_4[-2].real, 0, 0, l)
        s_2_3_result[i] += apply_angle_part_s_parallel(res_4[2].real, 0, 0, l)
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

          b_l   = f_p_el_for_p(el_nr,l,1,_k,z,nu_in,2,6)
          c_0_l = f_p_el_for_p(el_nr,l,0,_k,z,nu_in,2,6)
          c_2_l = f_p_el_for_p(el_nr,l,2,_k,z,nu_in,2,6)
          d_l   = f_p_el_for_p(el_nr,l,2,_k,z,nu_in,2,6)
          for runny in range(len(b_l)):
            b_l[runny],c_0_l[runny],c_2_l[runny],d_l[runny] = apply_angle_part_p_parallel(b_l[runny].real,
                                                                                 c_0_l[runny].real,
                                                                                 c_2_l[runny].real,
                                                                                 d_l[runny].real,
                                                                                 _k, t0, alp, l)
          p6_0_result[i] += third * (b_l[0] + c_0_l[0] + c_2_l[0] + d_l[0] + b_l[-1] + c_0_l[-1] + c_2_l[-1] + d_l[-1])
          p6_1_result[i] += third * (b_l[1] + c_0_l[1] + c_2_l[1] + d_l[1] + b_l[-2] + c_0_l[-2] + c_2_l[-2] + d_l[-2])
          p6_2_result[i] += third * (b_l[2] + c_0_l[2] + c_2_l[2] + d_l[2] + b_l[-3] + c_0_l[-3] + c_2_l[-3] + d_l[-3])
          p6_3_result[i] += third * (b_l[3] + c_0_l[3] + c_2_l[3] + d_l[3])
          
        s_0_result[i] = apply_angle_part_s_parallel(s_0_result[i],t0,alp,1)
        s_1_1_result[i] = apply_angle_part_s_parallel(s_1_1_result[i],t0,alp,1)
        s_1_2_result[i] = apply_angle_part_s_parallel(s_1_2_result[i],t0,alp,2)
        s_2_1_result[i] = apply_angle_part_s_parallel(s_2_1_result[i],t0,alp,1)
        s_2_2_result[i] = apply_angle_part_s_parallel(s_2_2_result[i],t0,alp,2)
        s_2_3_result[i] = apply_angle_part_s_parallel(s_2_3_result[i],t0,alp,3)


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
      
        #for p = 6 using the 0 angle applying a common angle function later
        b_l   = f_p_el_for_p(el_nr,l,1,1,z,nu_in,2,6) 
        c_0_l = f_p_el_for_p(el_nr,l,0,0,z,nu_in,2,6) 
        c_2_l = f_p_el_for_p(el_nr,l,2,2,z,nu_in,2,6) 
        fac = [0,1,2,3,2,1,0]
        for runny in range(len(b_l)):
          res_6[runny] = apply_angle_part_s_parallel(third * (  b_l[runny].real * alpha_coef(l,1,1,0,0) 
                                                     + c_0_l[runny].real * alpha_coef(l,0,0,0,0)
                                                     + c_2_l[runny].real * (alpha_coef(l,2,2,0,0) + beta_bar_coef(l,2,2,0,0)))
                                            , t0, alp, fac[runny])
        p6_0_result_M[i] += ((res_6[0]+res_6[-1]))
        p6_1_result_M[i] += ((res_6[1]+res_6[-2]))
        p6_2_result_M[i] += ((res_6[2]+res_6[-3]))
        p6_3_result_M[i] += ((res_6[3]))
        
  
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
    axes[0,0].scatter(x,g_result,s=20,facecolors='none',edgecolors='b',marker='^',label="Hönl fs_1^(0)") # type: ignore
    axes[0,0].scatter(x,h_result,s=20,facecolors='none',edgecolors='g',marker='^',label="Hönl fs_1^(2)") # type: ignore
    axes[0,0].scatter(x,i_result,s=20,facecolors='none',edgecolors='r',marker='^',label="Hönl fs_2^(2)") # type: ignore
    axes[0,0].plot(x,a_result,color='b')# type: ignore
    axes[0,0].plot(x,b_result,color='g')# type: ignore
    axes[0,0].plot(x,c_result,color='r')# type: ignore
    axes[0,0].plot(x,b0_result,linestyle='dotted')# type: ignore
    axes[0,0].plot(x,b1_result,linestyle='dotted')# type: ignore
    axes[0,0].plot(x,b2_result,linestyle='dotted')# type: ignore
    axes[0,0].plot(x,k_s,color='black')# type: ignore
    axes[0,0].legend()# type: ignore
    axes[0,0].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
  
    axes[1,0].scatter(x,j_result,s=20,facecolors='none',edgecolors='b',marker='^',label="EM fs_1^(0)")# type: ignore
    axes[1,0].scatter(x,k_result,s=20,facecolors='none',edgecolors='g',marker='^',label="EM fs_1^(2)")# type: ignore
    axes[1,0].scatter(x,l_result,s=20,facecolors='none',edgecolors='r',marker='^',label="EM fs_2^(2)")# type: ignore
    axes[1,0].plot(x,d_result,color='b')# type: ignore
    axes[1,0].plot(x,e_result,color='g')# type: ignore
    axes[1,0].plot(x,f_result,color='r')# type: ignore
    axes[1,0].plot(x,e0_result,linestyle='dotted')# type: ignore
    axes[1,0].plot(x,e1_result,linestyle='dotted')# type: ignore
    axes[1,0].plot(x,e2_result,linestyle='dotted')# type: ignore
    axes[1,0].plot(x,l_s,color='black')# type: ignore
    axes[1,0].legend()# type: ignore
    axes[1,0].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    
    axes[0,1].scatter(x,m_result,s=20,facecolors='none',edgecolors='b',marker='^',label="EM fp_1^(0)").legend_elements()# type: ignore
    axes[0,1].scatter(x,n_result,s=20,facecolors='none',edgecolors='g',marker='^',label="EM fp_1^(2)").legend_elements()# type: ignore
    axes[0,1].scatter(x,o_result,s=20,facecolors='none',edgecolors='r',marker='^',label="EM fp_2^(2)").legend_elements()# type: ignore
    ours1, = axes[0,1].plot(x,p0_result,color='b')# type: ignore
    ours2, = axes[0,1].plot(x,p2_0_result,color='g')# type: ignore
    ours3, = axes[0,1].plot(x,p2_1_result,color='r')# type: ignore
    ours4, = axes[0,1].plot(x,p4_0_result,linestyle='dotted')# type: ignore
    ours5, = axes[0,1].plot(x,p4_1_result,linestyle='dotted')# type: ignore
    ours6, = axes[0,1].plot(x,p4_2_result,linestyle='dotted')# type: ignore
    ours7, = axes[0,1].plot(x,p6_0_result,linestyle='dashed')# type: ignore
    ours8, = axes[0,1].plot(x,p6_1_result,linestyle='dashed')# type: ignore
    ours9, = axes[0,1].plot(x,p6_2_result,linestyle='dashed')# type: ignore
    ours10, = axes[0,1].plot(x,p6_3_result,linestyle='dashed')# type: ignore
    ours11, = axes[0,1].plot(x,l_p,color='black')# type: ignore
    axes[0,1].legend()# type: ignore
    axes[0,1].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    legend1 = axes[1,1].legend(# type: ignore
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
                    bbox_to_anchor=(-0.15,1.0))# type: ignore
    
    axes[1,1].plot(x,p0_result  ,color='b')# type: ignore
    axes[1,1].plot(x,p2_0_result,color='g')# type: ignore
    axes[1,1].plot(x,p2_1_result,color='r')# type: ignore
    axes[1,1].plot(x,p4_0_result,linestyle='dotted')# type: ignore
    axes[1,1].plot(x,p4_1_result,linestyle='dotted')# type: ignore
    axes[1,1].plot(x,p4_2_result,linestyle='dotted')# type: ignore
    axes[1,1].plot(x,p6_0_result,linestyle='dashed')# type: ignore
    axes[1,1].plot(x,p6_1_result,linestyle='dashed')# type: ignore
    axes[1,1].plot(x,p6_2_result,linestyle='dashed')# type: ignore
    axes[1,1].plot(x,p6_3_result,linestyle='dashed')# type: ignore
    axes[1,1].plot(x,l_p,color='k')# type: ignore
    axes[1,1].scatter(x,p0_result_M  ,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 with M")# type: ignore
    axes[1,1].scatter(x,p2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")# type: ignore
    axes[1,1].scatter(x,p2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")# type: ignore
    axes[1,1].scatter(x,p4_0_result_M,s=20,facecolors='none',edgecolor='b',marker='o',label="p=4 0/4 with M")# type: ignore
    axes[1,1].scatter(x,p4_1_result_M,s=20,facecolors='none',edgecolor='orange',marker='o',label="p=4 1/3 with M")# type: ignore
    axes[1,1].scatter(x,p4_2_result_M,s=20,facecolors='none',edgecolor='g',marker='o',label="p=4 2/2 with M")# type: ignore
    axes[1,1].scatter(x,p6_0_result_M,s=20,facecolors='none',edgecolor='b',marker='v',label="p=6 0/6 with M")# type: ignore
    axes[1,1].scatter(x,p6_1_result_M,s=20,facecolors='none',edgecolor='orange',marker='v',label="p=6 1/5 with M")# type: ignore
    axes[1,1].scatter(x,p6_2_result_M,s=20,facecolors='none',edgecolor='g',marker='v',label="p=6 2/4 with M")# type: ignore
    axes[1,1].scatter(x,p6_3_result_M,s=20,facecolors='none',edgecolor='k',marker='v',label="p=6 3/3 with M")# type: ignore
    axes[1,1].scatter(x,l_p_M,s=20,facecolors='none',edgecolor='k',marker='*',label="sum")# type: ignore
    axes[1,1].legend(loc='upper right')# type: ignore
    axes[1,1].add_artist(legend1)# type: ignore
  
    plt.subplots_adjust(left=0.025, bottom=0.04, right=1.0, top=1.0, wspace=0.15, hspace=0.05)
    fig.suptitle("alpha = {:4.2f}, theta_0 = {:4.2f}".format(alp,t0))

    print("S ORBITAL FUNCITONS")
    print("------------tau = 0-------------")
    print(a_result/s_0_result)
    print("------------tau = 1-------------")
    print(b_result/s_1_1_result)
    print(c_result/s_1_2_result)
    print("------------tau = 2-------------")
    print(b0_result/s_2_1_result)
    print(b1_result/s_2_2_result)
    print(b2_result/s_2_3_result)

    print("P ORBITALS BELOW")
    print(p4_2_result/p4_2_result_M)
    print(p4_1_result/p4_1_result_M)
    print(p4_0_result/p4_0_result_M)

    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!P=6!!!!!!!!!!!!!!!!!!!!!!!!")
    print(p6_3_result/p6_3_result_M)
    print(p6_2_result/p6_2_result_M)
    print(p6_1_result/p6_1_result_M)
    print(p6_0_result/p6_0_result_M)

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
    np.random.seed(int(np.ceil(time())))
    t0 = np.random.random(1)[0] * math.pi
    alp = np.random.random(1)[0] * math.pi * 2
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
        p0_result_M[i] += apply_angle_part_s_parallel(third * (  b_l   * al11[l]
                                                      + c_0_l * al00[l]
                                                      + c_2_l * (al22[l] + bbl22[l])
                                                     )
                                             , t0, alp, 1) /constant_factor
        
        #for p = 2 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l, 1,1, z, nu_in, 2, 2)
        c_0_l   = f_p_el_for_p(el_nr, l, 0,0, z, nu_in, 2, 2)
        c_2_l   = f_p_el_for_p(el_nr, l, 2,2, z, nu_in, 2, 2)
        fac = [1,2,1]
        for num in range(3):
          res_2[num] = apply_angle_part_s_parallel(third * (  b_l[num].real * al11[l]
                                                 + c_0_l[num].real * al00[l]
                                                 + c_2_l[num].real * (al22[l]+bbl22[l])
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
        M_p0_result_M[i] += apply_angle_part_s_parallel(third * (  b_l   * al11[l]
                                                        + c_0_l * al00[l]
                                                        + c_2_l * (al22[l] + bbl22[l])
                                                        )
                                             , t0, alp, 1) /constant_factor
        
        #for p = 2 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l, 1,1, z, nu_in, 3, 2)
        c_0_l   = f_p_el_for_p(el_nr, l, 0,0, z, nu_in, 3, 2)
        c_2_l   = f_p_el_for_p(el_nr, l, 2,2, z, nu_in, 3, 2)
        fac = [1,2,1]
        for num in range(3):
          res_2[num] = apply_angle_part_s_parallel(third * (  b_l[num].real * al11[l]
                                                 + c_0_l[num].real * al00[l]
                                                 + c_2_l[num].real * (al22[l]+bbl22[l])
                                                  )
                                                 , t0, alp, fac[num])
        M_p2_0_result_M[i] += ((res_2[0]+res_2[-1])) /constant_factor
        M_p2_1_result_M[i] += ((res_2[1])) /constant_factor

        #d-orbital
        #Calculate all angle dependant parts
        ms = [0,2,2,1,3,4,1,3]
        for _k in range(5):
          res = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
          for index,m in enumerate(ms):
            res[index] = f_d_el_for_p(el_nr, l, m, _k, z, nu_in, 3, 0)[0].real
          res = apply_angle_part_d_parallel(res, _k, t0, alp, l)
          M_d0_result[i] += 0.2 * sum(res) /constant_factor

          res = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                 [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                 [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]]
          for index,m in enumerate(ms):
            temp = f_d_el_for_p(el_nr, l, m, _k, z, nu_in, 3, 2)
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
        mults = [al00[l],bbl11[l]+al11[l],al22[l]+bbl22[l],bbl33[l]+al33[l],al11[l]]
        res = 0
        for nummy in range(len(ms)):
          r = f_d_el_for_p(el_nr, l, ms[nummy], ms[nummy], z, nu_in, 3, 0)
          res += r[0].real * mults[nummy]
        M_d0_result_M[i] += apply_angle_part_s_parallel(0.2 * res, t0, alp, 1)  /constant_factor
        
        #for p = 2 using the 0 angle applying a common angle function later
        res_2 = [0.0,0.0,0.0]
        fac = [1,2,1]
        for nummy,m in enumerate(ms):
          r = f_d_el_for_p(el_nr, l, m, m, z, nu_in, 3, 2)
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
    axes[0,0].scatter(x,g_result,s=20,facecolors='none',edgecolors='b',marker='^',label="Hönl fs_1^(0)")# type: ignore
    axes[0,0].scatter(x,h_result,s=20,facecolors='none',edgecolors='g',marker='^',label="Hönl fs_1^(2)")# type: ignore
    axes[0,0].scatter(x,i_result,s=20,facecolors='none',edgecolors='r',marker='^',label="Hönl fs_2^(2)")# type: ignore
    axes[0,0].plot(x,a_result,color='b')# type: ignore
    axes[0,0].plot(x,b_result,color='g')# type: ignore
    axes[0,0].plot(x,c_result,color='r')# type: ignore
    axes[0,0].plot(x,k_s,color='black')# type: ignore
    axes[0,0].legend()# type: ignore
    axes[0,0].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[0,0].set_title("K-shell s-electrons", y=1.0, pad=-14)# type: ignore
    #L-shell s-orbital
    axes[1,0].scatter(x,j_result,s=20,facecolors='none',edgecolors='b',marker='^',label="EM fs_1^(0)")# type: ignore
    axes[1,0].scatter(x,k_result,s=20,facecolors='none',edgecolors='g',marker='^',label="EM fs_1^(2)")# type: ignore
    axes[1,0].scatter(x,l_result,s=20,facecolors='none',edgecolors='r',marker='^',label="EM fs_2^(2)")# type: ignore
    axes[1,0].plot(x,d_result,color='b')# type: ignore
    axes[1,0].plot(x,e_result,color='g')# type: ignore
    axes[1,0].plot(x,f_result,color='r')# type: ignore
    axes[1,0].plot(x,l_s,color='black')# type: ignore
    axes[1,0].legend()# type: ignore
    axes[1,0].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[1,0].set_title("L-shell s-electrons", y=1.0, pad=-14)# type: ignore
    #M-shell s-orbital
    axes[2,0].scatter(x,x_result,s=20,facecolors='none',edgecolors='b',marker='^',label="WA fs_1^(0)")# type: ignore
    axes[2,0].scatter(x,y_result,s=20,facecolors='none',edgecolors='g',marker='^',label="WA fs_1^(2)")# type: ignore
    axes[2,0].scatter(x,z_result,s=20,facecolors='none',edgecolors='r',marker='^',label="WA fs_2^(2)")# type: ignore
    axes[2,0].plot(x,Ms0_result,color='b')# type: ignore
    axes[2,0].plot(x,Ms1_result,color='g')# type: ignore
    axes[2,0].plot(x,Ms2_result,color='r')# type: ignore
    axes[2,0].plot(x,m_s,color='black')# type: ignore
    axes[2,0].legend()# type: ignore
    axes[2,0].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[2,0].set_title("M-shell s-electrons", y=1.0, pad=-14)# type: ignore
    
    ours1, = axes[1,1].plot(x,p0_result,color='b')# type: ignore
    ours2, = axes[1,1].plot(x,p2_0_result,color='g')# type: ignore
    ours3, = axes[1,1].plot(x,p2_1_result,color='r')# type: ignore
    ours11, = axes[1,1].plot(x,l_p,color='black')# type: ignore
    legend1 = axes[1,1].legend(# type: ignore
                    handles=[ours1,
                             ours2,ours3,
                             ours11],
                    labels=["0/0",
                            "0/2 + 2/0","1/1",
                            "sum"],
                    loc = 'upper left',
                    bbox_to_anchor=(0.5,1.5))# type: ignore
    
    axes[1,1].plot(x,p0_result,  color='b')# type: ignore
    axes[1,1].plot(x,p2_0_result,color='g')# type: ignore
    axes[1,1].plot(x,p2_1_result,color='r')# type: ignore
    axes[1,1].scatter(x,m_result,s=40,facecolors='none',edgecolors='b',marker='o',label="EM fp_1^(0)")# type: ignore
    axes[1,1].scatter(x,n_result,s=40,facecolors='none',edgecolors='g',marker='o',label="EM fp_1^(2)")# type: ignore
    axes[1,1].scatter(x,o_result,s=40,facecolors='none',edgecolors='r',marker='o',label="EM fp_2^(2)")# type: ignore
    axes[1,1].plot(x,l_p,color='black')# type: ignore
    axes[1,1].scatter(x,p0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")# type: ignore
    axes[1,1].scatter(x,p2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")# type: ignore
    axes[1,1].scatter(x,p2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")# type: ignore
    axes[1,1].scatter(x,l_p_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")# type: ignore
    axes[1,1].legend()# type: ignore
    axes[1,1].add_artist(legend1)# type: ignore
    axes[1,1].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[1,1].set_title("L-shell p-electrons", y=1.0, pad=-14)# type: ignore
    
    axes[2,1].plot(x,M_p0_result,  color='b')# type: ignore
    axes[2,1].plot(x,M_p2_0_result,color='g')# type: ignore
    axes[2,1].plot(x,M_p2_1_result,color='r')# type: ignore
    axes[2,1].plot(x,m_p,color='black')# type: ignore
    axes[2,1].scatter(x,M_p0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")# type: ignore
    axes[2,1].scatter(x,M_p2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")# type: ignore
    axes[2,1].scatter(x,M_p2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")# type: ignore
    axes[2,1].scatter(x,WA_p0_result,s=40,facecolors='none',edgecolor='b',marker='o',label="Wagenfeld p-els D")# type: ignore
    axes[2,1].scatter(x,m_p_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")# type: ignore
    axes[2,1].legend()# type: ignore
    axes[2,1].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[2,1].set_title("M-shell p-electrons", y=1.0, pad=-14)# type: ignore

    axes[0,1].set_axis_off()# type: ignore
    axes[0,2].set_axis_off()# type: ignore
    axes[1,2].set_axis_off()# type: ignore

    axes[2,2].plot(x,M_d0_result,  color='b')# type: ignore
    axes[2,2].plot(x,M_d2_0_result,color='g')# type: ignore
    axes[2,2].plot(x,M_d2_1_result,color='r')# type: ignore
    axes[2,2].plot(x,m_d,color='black')# type: ignore
    axes[2,2].scatter(x,M_d0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")# type: ignore
    axes[2,2].scatter(x,M_d2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")# type: ignore
    axes[2,2].scatter(x,M_d2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")# type: ignore
    axes[2,2].scatter(x,WA_d0_result,s=40,facecolors='none',edgecolor='b',marker='o',label="Wagenfeld d-els D")# type: ignore
    axes[2,2].scatter(x,m_d_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")# type: ignore
    axes[2,2].legend()# type: ignore
    axes[2,2].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[2,2].set_title("M-shell d-electrons", y=1.0, pad=-14)# type: ignore

    print(M_d0_result - M_d0_result_M)
  
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
    np.random.seed(int(np.ceil(time())))
    t0 = np.random.random(1)[0] * math.pi
    #t0 = math.pi / 2
    alp = np.random.random(1)[0] * math.pi * 2
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
        p0_result_M[i] += apply_angle_part_s_orthogonal(third * (  b_l   * al11[l]
                                                      + c_0_l * al00[l]
                                                      + c_2_l * (al22[l] + bbl22[l])
                                                     )
                                             , t0, alp, 1)
        
        #for p = 2 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l,1, 1, z, nu_in, 2, 2)
        c_0_l   = f_p_el_for_p(el_nr, l,0, 0, z, nu_in, 2, 2)
        c_2_l   = f_p_el_for_p(el_nr, l,2, 2, z, nu_in, 2, 2)
        fac = [1,2,1]
        for num in range(3):
          res_2[num] = apply_angle_part_s_orthogonal(third * (  b_l[num].real * al11[l]
                                                 + c_0_l[num].real * al00[l]
                                                 + c_2_l[num].real * (al22[l]+bbl22[l])
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
        M_p0_result_M[i] += apply_angle_part_s_orthogonal(third * (  b_l   * al11[l]
                                                        + c_0_l * al00[l]
                                                        + c_2_l * (al22[l] + bbl22[l])
                                                        )
                                             , t0, alp, 1)
        
        #for p = 2 using the 0 angle applying a common angle function later
        b_l     = f_p_el_for_p(el_nr, l,1, 1, z, nu_in, 3, 2)
        c_0_l   = f_p_el_for_p(el_nr, l,0, 0, z, nu_in, 3, 2)
        c_2_l   = f_p_el_for_p(el_nr, l,2, 2, z, nu_in, 3, 2)
        fac = [1,2,1]
        for num in range(3):
          res_2[num] = apply_angle_part_s_orthogonal(third * (  b_l[num].real * al11[l]
                                                 + c_0_l[num].real * al00[l]
                                                 + c_2_l[num].real * (al22[l] + bbl22[l])
                                                  )
                                                 , t0, alp, fac[num])
        M_p2_0_result_M[i] += ((res_2[0]+res_2[-1]))
        M_p2_1_result_M[i] += ((res_2[1]))

        #d-orbital
        #Calculate all angle dependant parts
        ms = [0,2,2,1,3,4,1,3]
        for _k in range(5):
          res = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
          for index,m in enumerate(ms):
            res[index] = f_d_el_for_p(el_nr, l, m, _k, z, nu_in, 3, 0)[0].real
          res = apply_angle_part_d_orthogonal(res, _k, t0, alp, l)
          M_d0_result[i] += 0.2 * sum(res)

          res = [[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                 [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                 [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]]
          for index,m in enumerate(ms):
            temp = f_d_el_for_p(el_nr, l, m, _k, z, nu_in, 3, 2)
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
        mults = [al00[l],bbl11[l]+al11[l],al22[l]+bbl22[l],bbl33[l]+al33[l],al11[l]]
        res = 0
        for nummy in range(len(ms)):
          r = f_d_el_for_p(el_nr, l, ms[nummy], ms[nummy], z, nu_in, 3, 0)
          res += r[0].real * mults[nummy]
        M_d0_result_M[i] += apply_angle_part_s_orthogonal(0.2 * res, t0, alp, 1)
        
        #for p = 2 using the 0 angle applying a common angle function later
        res_2 = [0.0,0.0,0.0]
        fac = [1,2,1]
        for nummy,m in enumerate(ms):
          r = f_d_el_for_p(el_nr, l, m, m, z, nu_in, 3, 2)
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
    axes[0,0].scatter(x,g_result,s=20,facecolors='none',edgecolors='b',marker='^',label="Hönl fs_1^(0)")# type: ignore
    axes[0,0].scatter(x,h_result,s=20,facecolors='none',edgecolors='g',marker='^',label="Hönl fs_1^(2)")# type: ignore
    axes[0,0].scatter(x,i_result,s=20,facecolors='none',edgecolors='r',marker='^',label="Hönl fs_2^(2)")# type: ignore
    axes[0,0].plot(x,a_result,color='b')# type: ignore
    axes[0,0].plot(x,b_result,color='g')# type: ignore
    axes[0,0].plot(x,c_result,color='r')# type: ignore
    axes[0,0].plot(x,k_s,color='black')# type: ignore
    axes[0,0].legend()# type: ignore
    axes[0,0].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[0,0].set_title("K-shell s-electrons", y=1.0, pad=-14)# type: ignore
    #L-shell s-orbital
    axes[1,0].scatter(x,j_result,s=20,facecolors='none',edgecolors='b',marker='^',label="EM fs_1^(0)")# type: ignore
    axes[1,0].scatter(x,k_result,s=20,facecolors='none',edgecolors='g',marker='^',label="EM fs_1^(2)")# type: ignore
    axes[1,0].scatter(x,l_result,s=20,facecolors='none',edgecolors='r',marker='^',label="EM fs_2^(2)")# type: ignore
    axes[1,0].plot(x,d_result,color='b')# type: ignore
    axes[1,0].plot(x,e_result,color='g')# type: ignore
    axes[1,0].plot(x,f_result,color='r')# type: ignore
    axes[1,0].plot(x,l_s,color='black')# type: ignore
    axes[1,0].legend()# type: ignore
    axes[1,0].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[1,0].set_title("L-shell s-electrons", y=1.0, pad=-14)# type: ignore
    #M-shell s-orbital
    axes[2,0].scatter(x,x_result,s=20,facecolors='none',edgecolors='b',marker='^',label="WA fs_1^(0)")# type: ignore
    axes[2,0].scatter(x,y_result,s=20,facecolors='none',edgecolors='g',marker='^',label="WA fs_1^(2)")# type: ignore
    axes[2,0].scatter(x,z_result,s=20,facecolors='none',edgecolors='r',marker='^',label="WA fs_2^(2)")# type: ignore
    axes[2,0].plot(x,Ms0_result,color='b')# type: ignore
    axes[2,0].plot(x,Ms1_result,color='g')# type: ignore
    axes[2,0].plot(x,Ms2_result,color='r')# type: ignore
    axes[2,0].plot(x,m_s,color='black')# type: ignore
    axes[2,0].legend()# type: ignore
    axes[2,0].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[2,0].set_title("M-shell s-electrons", y=1.0, pad=-14)# type: ignore
    
    ours1, = axes[1,1].plot(x,p0_result,color='b')# type: ignore
    ours2, = axes[1,1].plot(x,p2_0_result,color='g')# type: ignore
    ours3, = axes[1,1].plot(x,p2_1_result,color='r')# type: ignore
    ours11, = axes[1,1].plot(x,l_p,color='black')# type: ignore
    legend1 = axes[1,1].legend(# type: ignore
                    handles=[ours1,
                             ours2,ours3,
                             ours11],
                    labels=["0/0",
                            "0/2 + 2/0","1/1",
                            "sum"],
                    loc = 'upper left',
                    bbox_to_anchor=(0.5,1.5))# type: ignore
    
    axes[1,1].plot(x,p0_result,  color='b')# type: ignore
    axes[1,1].plot(x,p2_0_result,color='g')# type: ignore
    axes[1,1].plot(x,p2_1_result,color='r')# type: ignore
    axes[1,1].scatter(x,m_result,s=40,facecolors='none',edgecolors='b',marker='o',label="EM fp_1^(0)")# type: ignore
    axes[1,1].scatter(x,n_result,s=40,facecolors='none',edgecolors='g',marker='o',label="EM fp_1^(2)")# type: ignore
    axes[1,1].scatter(x,o_result,s=40,facecolors='none',edgecolors='r',marker='o',label="EM fp_2^(2)")# type: ignore
    axes[1,1].plot(x,l_p,color='black')# type: ignore
    axes[1,1].scatter(x,p0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")# type: ignore
    axes[1,1].scatter(x,p2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")# type: ignore
    axes[1,1].scatter(x,p2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")# type: ignore
    axes[1,1].scatter(x,l_p_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")# type: ignore
    axes[1,1].legend()# type: ignore
    axes[1,1].add_artist(legend1)# type: ignore
    axes[1,1].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[1,1].set_title("L-shell p-electrons", y=1.0, pad=-14)# type: ignore
    
    axes[2,1].plot(x,M_p0_result,  color='b')# type: ignore
    axes[2,1].plot(x,M_p2_0_result,color='g')# type: ignore
    axes[2,1].plot(x,M_p2_1_result,color='r')# type: ignore
    axes[2,1].plot(x,m_p,color='black')# type: ignore
    axes[2,1].scatter(x,M_p0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")# type: ignore
    axes[2,1].scatter(x,M_p2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")# type: ignore
    axes[2,1].scatter(x,M_p2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")# type: ignore
    axes[2,1].scatter(x,WA_p0_result,s=40,facecolors='none',edgecolor='b',marker='o',label="Wagenfeld p-els D")# type: ignore
    axes[2,1].scatter(x,m_p_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")# type: ignore
    axes[2,1].legend()# type: ignore
    axes[2,1].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[2,1].set_title("M-shell p-electrons", y=1.0, pad=-14)# type: ignore

    axes[0,1].set_axis_off()# type: ignore
    axes[0,2].set_axis_off()# type: ignore
    axes[1,2].set_axis_off()# type: ignore

    axes[2,2].plot(x,M_d0_result,  color='b')# type: ignore
    axes[2,2].plot(x,M_d2_0_result,color='g')# type: ignore
    axes[2,2].plot(x,M_d2_1_result,color='r')# type: ignore
    axes[2,2].plot(x,m_d,color='black')# type: ignore
    axes[2,2].scatter(x,M_d0_result_M,s=20,facecolors='none',edgecolor='b',marker='^',label="p=0 0/0 with M")# type: ignore
    axes[2,2].scatter(x,M_d2_0_result_M,s=20,facecolors='none',edgecolor='g',marker='^',label="p=2 0/2 with M")# type: ignore
    axes[2,2].scatter(x,M_d2_1_result_M,s=20,facecolors='none',edgecolor='r',marker='^',label="p=2 1/1 with M")# type: ignore
    axes[2,2].scatter(x,WA_d0_result,s=40,facecolors='none',edgecolor='b',marker='o',label="Wagenfeld d-els D")# type: ignore
    axes[2,2].scatter(x,m_d_M,s=20,facecolors='none',edgecolor='black',marker='*',label="sum")# type: ignore
    axes[2,2].legend()# type: ignore
    axes[2,2].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
    axes[2,2].set_title("M-shell d-electrons", y=1.0, pad=-14)# type: ignore
  
    plt.subplots_adjust(left=0.04, bottom=0.04, right=1.0, top=0.95, wspace=0.15, hspace=0.1)
    fig.suptitle("ORTHO alpha = {:4.2f}, theta_0 = {:4.2f}".format(alp,t0))
    plt.show()
    exit()
    
  if form_factor_test:
    #import cProfile
    #from pstats import Stats, SortKey
    #with cProfile.Profile() as pr:
      print("Performing Disp Correction Calculation")
      #axis x is in nu
      x = np.logspace(math.log10(speed_of_light*1E10/lam_min), 
                      math.log10(speed_of_light*1E10/lam_max), 
                      steps)
  
      np.random.seed(int(np.ceil(time())))
      t0 = np.random.random(1)[0] * math.pi
      if overwrite_theta == True:
        t0 = 0
      TS = math.sqrt((1.0+np.cos(t0)**2)/2.0)
      e = elements[el_nr]
      
      factor_x = 2*math.pi
  
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
        
        res = pool.starmap(calc_stuff_only_sums, zip(x, repeat(0), repeat(l_max), repeat(p_max), repeat(el_nr)))
        for i,r in enumerate(res):
          zero_angle[i] = r * factor_x
      
      second_x = speed_of_light / x * 1E10
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
      delta_K  = one_minus_delta_edge(el_nr,1,0,0)
      delta_l1 = one_minus_delta_edge(el_nr,2,0,0)
      delta_l2 = one_minus_delta_edge(el_nr,2,1,0.5)
      delta_l3 = one_minus_delta_edge(el_nr,2,1,1.5)
      delta_m1 = one_minus_delta_edge(el_nr,3,0,0)
      delta_m2 = one_minus_delta_edge(el_nr,3,1,0.5)
      delta_m3 = one_minus_delta_edge(el_nr,3,1,1.5)
      delta_m4 = one_minus_delta_edge(el_nr,3,2,1.5)
      delta_m5 = one_minus_delta_edge(el_nr,3,2,2.5)
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
        ax.scatter(second_x,o,s=20,facecolors='none',edgecolors='b',marker='^',label="orthogonal")
        ax.scatter(second_x,p,s=20,facecolors='none',edgecolors='g',marker='^',label="parallel")
        ax.plot(second_x,t,color='r',label="total")
        ax.axhline(y=0,linestyle='dashed',color='gray')
        ax.set_title(name, y=1.0, pad=-14)
  
      def plot(ax,t,p1=None,p2=None,name=""):
        if p1 is not None:
          ax.scatter(second_x,p1,s=20,facecolors='none',edgecolors='b',marker='^',label="subshell1")
        if p2 is not None:
          ax.scatter(second_x,p2,s=20,facecolors='none',edgecolors='g',marker='^',label="subshell2")
        ax.plot(second_x,t,color='r',label="total")
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
  
      res_K = k_s*factor_x/TS
      res_L = (l_s+l_p1+2*l_p3)*factor_x/TS
      res_M = (m_s+m_p1+2*m_p3+3*m_d3+2*m_d5)*factor_x/TS
      resy = res_K+res_L+res_M
  
      axes[0,2].plot(second_x,resy,color='k',label="Norbert & Florian")# type: ignore
      axes[0,2].axhline(y=0,linestyle='dashed',color='gray')# type: ignore
      axes[0,2].plot(second_x,brennan_fdp,color='r',label="brennan")# type: ignore
      axes[0,2].plot(second_x,hönl,color='b',label="hönl")# type: ignore
  
      axes[0,1].scatter(second_x,fdp_nach_hoenl,facecolors='none',edgecolors='g',marker='^',label="K purely paper")# type: ignore
      axes[0,1].scatter(second_x,fdp_nach_EM,facecolors='none',edgecolors='b',marker='^',label="L purely paper")# type: ignore
      axes[0,1].scatter(second_x,fdp_nach_Wa,facecolors='none',edgecolors='r',marker='^',label="M purely paper")# type: ignore
      
      #axes[0,2].scatter(second_x,fdp_nach_Wa+fdp_nach_EM+fdp_nach_hoenl,facecolors='none',edgecolors='k',marker='^',label="total purely paper")
  
      axes[1,2].plot(second_x,resy/brennan_fdp,color='k',label="ratio N&F/brennan")# type: ignore
      axes[1,2].plot(second_x,hönl/brennan_fdp,color='b',label="ratio hönl/brennan")# type: ignore
      axes[1,2].plot(second_x,resy/hönl,color='g',label="ratio N&F/hönl")# type: ignore
      axes[1,2].plot(second_x,(fdp_nach_Wa+fdp_nach_EM+fdp_nach_hoenl)/hönl,color='r',label="ratio hönl/hönl_integral")# type: ignore
      axes[1,2].axhline(y=1,linestyle='dashed',color='gray')# type: ignore
  
      axes[0,1].plot(second_x,res_K,label="K-shell electrons")# type: ignore
      axes[0,1].plot(second_x,res_L,label="L-shell electrons")# type: ignore
      axes[0,1].plot(second_x,res_M,label="M-shell electrons")# type: ignore
  
      axes[0,1].plot(second_x,hönl_K,linestyle='dashed',label="Hönl K-shell")# type: ignore
      axes[0,1].plot(second_x,hönl_L,linestyle='dashed',label="Hönl L-shell")# type: ignore
      axes[0,1].plot(second_x,hönl_M,linestyle='dashed',label="Hönl M-shell")# type: ignore
  
      axes[0,1].axhline(y=1,linestyle='dashed',color='gray')# type: ignore
      print("Starting angle calculations!")
      for i in range(1,theta_steps+1):
        theta = i*(1./theta_steps)*math.pi
        temp: npt.NDArray = np.zeros_like(zero_angle)
        print("{}/{} steps".format(i,theta_steps))
        with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
          res = pool.starmap(calc_stuff_only_sums, zip(x, repeat(theta), repeat(l_max), repeat(p_max), repeat(el_nr)))
          for i,r in enumerate(res):
            temp[i] = r
          temp *= factor_x
          axes[1,1].plot(second_x,zero_angle/temp * math.sqrt((1+np.cos(theta)**2)/2),label=r"$\theta_0 = ${:4.2f}$\pi$".format(theta/math.pi))# type: ignore
      print(factor_x)
      axes[1,1].plot(second_x,zero_angle/resy,color='k',label=r"$\theta_0 = ${:5.3f}$\pi$".format(t0/math.pi))# type: ignore
      axes[1,1].set_title(r"Ratio between $\theta_0=0$ and $\theta_0$", y=1.0, pad=-14)# type: ignore
      print("TS = " + str(TS))
  
      ticks = np.array([0.2,0.3,0.5,0.7,1.0,1.5,2.0,3.0])
  
      for ax in fig.get_axes():
        ax.set_xscale('log')
        ax.set_xlim(lam_min,lam_max)
        ax.set_xticks(ticks)
        ax.get_xaxis().get_major_formatter().labelOnlyBase = False
        ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
        ax.legend()
  
      ##This adds the frequency in ExaHz on top of the 1,2 axes
      #def tickfunc(x):
      #  N = speed_of_light / x * 1E-8 * h
      #  return  ["%.3f" %z for z in N]
      #ax2 = axes[1,2].twiny()# type: ignore
      #ax2.set_xscale('log')
      #ax2.set_xlim(lam_min,lam_max)
      #ax2.set_xticks(ticks)
      #ax2.set_xticklabels(tickfunc(ticks))
  
      
      fig.suptitle(r"$\theta_0$ = {:4.2f} $\pi$".format(t0/math.pi))
      plt.subplots_adjust(left=0.05, bottom=0.04, right=1.0, top=0.95, wspace=0.10, hspace=0.175)
      plt.show()
    #with open('factor_profiling_stats.txt', 'w') as stream:
    #    stats = Stats(pr, stream=stream)
    #    stats.strip_dirs()
    #    stats.sort_stats('time')
    #    stats.dump_stats('.prof_stats')
    #    stats.print_stats()
      exit()
  
  if numplot == True: 

    print("performing comparison of numerical and ealuation at a given energy")
    x = np.linspace(0.9,1.1,steps)
    eta_nach_hoenl = np.zeros_like(x)
    eta_nach_integral = np.zeros_like(x)
    eta_nach_unserem = np.zeros_like(x)
    nprime = np.zeros_like(x)
    inprime = np.zeros_like(x)
    K_recursive = np.zeros_like(x)
    K_recursive_imag = np.zeros_like(x)
    SE = np.zeros_like(x)
    
    nu_0 = (get_ionization_energy_1s(el_nr) / h)
    Z_s_sq = get_Zeff_1s(el_nr)**2
    e_ion = get_ionization_energy_1s(el_nr)
    #nu_k = e_ion / h
    omdelta_K = one_minus_delta_edge(el_nr, 1, 0, 0) #21a in Hoenl

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
      res = pool.starmap(eta_K, zip(x, repeat(1), repeat(nu_0),repeat(omdelta_K)))
      for i,r in enumerate(res):
        eta_nach_hoenl[i] = r
      res = pool.starmap(sugiura_k_purely_imag, zip(repeat(el_nr), x*get_ionization_energy_1s(el_nr)))
      for i,r in enumerate(res):
        eta_nach_integral[i] = r
      #res = pool.starmap(eta_K_unsers, zip(x,repeat(nu_0), repeat(el_nr)))
      #for i,r in enumerate(res):
      #  eta_nach_unserem[i] = r
        
      #from matrix_coefficients_v2 import K_recursive_from_z
      #res = pool.starmap(n_prime_from_z, zip(x,repeat(1)))
      #for i,r in enumerate(res):
      #  nprime[i] = r.real
      #  inprime[i] = r.imag
      #  
      #res = pool.starmap(K_recursive_from_z, zip(repeat(0),repeat(1), repeat(b(1,0,el_nr)), x,repeat(1)))
      #for i,r in enumerate(res):
      #  K_recursive[i] = r.real
      #  K_recursive_imag[i] = r.imag
      #res = pool.starmap(sugiura_exps, zip(x,repeat(1)))
      #for i,r in enumerate(res):
      #  SE[i] = r

    eta_nach_integral /= max(eta_nach_hoenl)  
    eta_nach_hoenl /= max(eta_nach_hoenl)
    

    fig, axes = plt.subplots(1,1)
    axes.plot(x,eta_nach_hoenl,label="without damping")# type: ignore
    #axes[0,0].plot(x,eta_nach_unserem,label="unsers")# type: ignore
    axes.plot(x,eta_nach_integral,"--",label="with damping")# type: ignore
    axes.set_ylabel(r"Normalized $f''$")
    axes.set_xlabel(r"$\frac{E_{X-ray}}{E_{Edge}}$")
    axes.legend()# type: ignore

    #ratio1 = np.zeros_like(eta_nach_hoenl)
    #ratio2 = np.zeros_like(eta_nach_hoenl)
    #for i,e in enumerate(eta_nach_hoenl):
    #  ratio1[i] = eta_nach_hoenl[i]  /eta_nach_integral[i]
    #  ratio2[i] = eta_nach_unserem[i]/eta_nach_integral[i]
    #  if x[i] <= 1:
    #    ratio1[i] += 1
    #    ratio2[i] += 1
#
    #axes[1,0].plot(x,ratio1,label="ratio")# type: ignore
    #axes[1,0].plot(x,ratio2,label="ratio unsers")# type: ignore
    #axes[1,0].legend()# type: ignore
#
    ##axes[0,0].axvline(x=1/omdelta_K,linestyle='dashed',color='gray')# type: ignore
    #axes[1,0].axvline(x=1/omdelta_K,linestyle='dashed',color='gray')# type: ignore
    #axes[1,1].axvline(x=1,linestyle='dashed',color='gray')# type: ignore
    #axes[0,1].axvline(x=1,linestyle='dashed',color='gray')# type: ignore
    #
    #axes[0,1].plot(x,nprime,label="Re prime(z)")# type: ignore
    #axes[0,1].plot(x,inprime,label="Im nprime(z)")# type: ignore
    #axes[0,1].legend()# type: ignore
    #axes[1,1].plot(x,SE*1E-10,label="sugiura(z)")# type: ignore
    #axes[1,1].plot(x,K_recursive,label="Re K(z)")# type: ignore
    #axes[1,1].plot(x,K_recursive_imag,label="Im K(z)")# type: ignore
    #axes[1,1].legend()# type: ignore
    #print(eta_nach_hoenl[-1]/eta_nach_integral[-1], eta_nach_unserem[-1]/eta_nach_integral[-1])
    
    fig.set_figheight(inches(17))
    fig.set_figwidth(inches(22))
    
    fig.savefig('integral_effect.png',dpi=300,bbox_inches='tight')
    exit()
  
  if theta_plot == True: 

    print("performing plot against theta angle at molybdnum wavelength")
    x = np.linspace(0,math.pi,steps)
    p0 = np.zeros_like(x)
    p2 = np.zeros_like(x)
    p4 = np.zeros_like(x)
    p6 = np.zeros_like(x)

    nu = speed_of_light*1E10/0.71078

    TS = np.sqrt((1.0+np.cos(x)**2)/2.0)

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
      res  = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(0), repeat(40)))
      res2 = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(2), repeat(40)))
      res4 = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(4), repeat(40)))
      res6 = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(6), repeat(40)))
      for i,r in enumerate(res):
        p0[i] = r
        p2[i] = res2[i]
        p4[i] = res4[i]
        p6[i] = res6[i]

    p0 /= TS
    p2 /= TS
    p4 /= TS
    p6 /= TS
    
    p0 *= 2*math.pi
    p2 *= 2*math.pi
    p4 *= 2*math.pi
    p6 *= 2*math.pi
    

    angle = x / math.pi * 180

    fig, axes = plt.subplots(1,1)
    axes.plot(angle,p0,"k",label="f'' p=0")# type: ignore
    axes.plot(angle,p2,"r-.",label="f'' p=2")# type: ignore
    axes.plot(angle,p4,"b--",label="f'' p=4")# type: ignore
    axes.plot(angle,p6,"g:",label="f'' p=6")# type: ignore
    axes.set_ylabel(r"$f''$ /$e$")
    axes.set_xlabel(r"$2\theta$")
    axes.axvline(x=50,linestyle='dashed',color='gray',label='Min. IUCr')# type: ignore
    axes.legend()# type: ignore
    
    fig.set_figheight(inches(17))
    fig.set_figwidth(inches(22))
    
    fig.savefig('Zr_atMo.png',dpi=300,bbox_inches='tight')

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
      res  = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(0), repeat(38)))
      res2 = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(2), repeat(38)))
      res4 = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(4), repeat(38)))
      res6 = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(6), repeat(38)))
      for i,r in enumerate(res):
        p0[i] = r
        p2[i] = res2[i]
        p4[i] = res4[i]
        p6[i] = res6[i]

    p0 /= TS
    p2 /= TS
    p4 /= TS
    p6 /= TS
    p0 *= 2*math.pi
    p2 *= 2*math.pi
    p4 *= 2*math.pi
    p6 *= 2*math.pi
    fig, axes = plt.subplots(1,1)
    axes.plot(angle,p0,"k",label="f'' p=0")# type: ignore
    axes.plot(angle,p2,"r-.",label="f'' p=2")# type: ignore
    axes.plot(angle,p4,"b--",label="f'' p=4")# type: ignore
    axes.plot(angle,p6,"g:",label="f'' p=6")# type: ignore
    axes.set_ylabel(r"$f''$ /$e$")
    axes.set_xlabel(r"$2\theta$")
    axes.axvline(x=50,linestyle='dashed',color='gray',label='Min. IUCr')# type: ignore
    axes.legend()# type: ignore
    
    fig.set_figheight(inches(17))
    fig.set_figwidth(inches(22))
    
    fig.savefig('Sr_atMo.png',dpi=300,bbox_inches='tight')
    
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
      res  = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(0), repeat(92)))
      res2 = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(2), repeat(92)))
      res4 = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(4), repeat(92)))
      res6 = pool.starmap(calc_stuff_only_sums, zip(repeat(nu), x, repeat(10), repeat(6), repeat(92)))
      for i,r in enumerate(res):
        p0[i] = r
        p2[i] = res2[i]
        p4[i] = res4[i]
        p6[i] = res6[i]

    p0 /= TS
    p2 /= TS
    p4 /= TS
    p6 /= TS
    p0 *= 2*math.pi
    p2 *= 2*math.pi
    p4 *= 2*math.pi
    p6 *= 2*math.pi
    fig, axes = plt.subplots(1,1)
    axes.plot(angle,p0,"k",label="f'' p=0")# type: ignore
    axes.plot(angle,p2,"r-.",label="f'' p=2")# type: ignore
    axes.plot(angle,p4,"b--",label="f'' p=4")# type: ignore
    axes.plot(angle,p6,"g:",label="f'' p=6")# type: ignore
    axes.set_ylabel(r"$f''$ /$e$")
    axes.set_xlabel(r"$2\theta$")
    axes.axvline(x=50,linestyle='dashed',color='gray',label='Min. IUCr')# type: ignore
    axes.legend()# type: ignore
    
    fig.set_figheight(inches(17))
    fig.set_figwidth(inches(22))
    
    fig.savefig('U_atMo.png',dpi=300,bbox_inches='tight')
    exit()
  if fcf_test == True:
    
    cif = """
data_Sr
_chemical_formula_moiety           'Sr1'
_chemical_formula_sum              'Sr1'
_chemical_formula_weight           87.62
loop_
  _atom_type_symbol
  _atom_type_scat_dispersion_real
  _atom_type_scat_dispersion_imag
  _atom_type_scat_source
  _atom_type_scat_dispersion_source
 Sr 0.35 0.522 'NoSpherA2: Chem.Sci. 2021, DOI:10.1039/D0SC05526C'
 'Brennan, Cowan, Rev. Sci. Instrum., 1992, 63, 850'

_space_group_crystal_system        'triclinic'
_space_group_IT_number             1
_space_group_name_H-M_alt          'P 1'
_space_group_name_Hall             'P 1'
loop_
  _space_group_symop_id
  _space_group_symop_operation_xyz
 1 x,y,z

_symmetry_Int_Tables_number        1
_cell_length_a                     10.000(2)
_cell_length_b                     10.000(3)
_cell_length_c                     10.000(2)
_cell_angle_alpha                  90.000(2)
_cell_angle_beta                   90.000(2)
_cell_angle_gamma                  90.000(2)
_cell_volume                       1000.0(5)
_cell_formula_units_Z              1
_cell_measurement_temperature      293(2)
_exptl_crystal_F_000               38.000
_diffrn_reflns_av_R_equivalents    0.0000
_diffrn_reflns_av_unetI/netI       0.0311
_diffrn_reflns_Laue_measured_fraction_full  0.9981
_diffrn_reflns_Laue_measured_fraction_max  0.9988
_diffrn_reflns_limit_h_max         15
_diffrn_reflns_limit_h_min         -15
_diffrn_reflns_limit_k_max         15
_diffrn_reflns_limit_k_min         -15
_diffrn_reflns_limit_l_max         15
_diffrn_reflns_limit_l_min         -15
_diffrn_reflns_number              7595
_diffrn_reflns_point_group_measured_fraction_full  0.9981
_diffrn_reflns_point_group_measured_fraction_max  0.9988
_diffrn_reflns_theta_full          25.2436
_diffrn_reflns_theta_max           30.05
_diffrn_reflns_theta_min           3.05
_diffrn_measured_fraction_theta_full  0.9981
_diffrn_measured_fraction_theta_max  0.9988
_diffrn_radiation_type             'Mo K\a'
_diffrn_radiation_wavelength       0.71078
_reflns_Friedel_coverage           0.0
_reflns_limit_h_max                15
_reflns_limit_h_min                -15
_reflns_limit_k_max                15
_reflns_limit_k_min                -15
_reflns_limit_l_max                15
_reflns_limit_l_min                -15
_reflns_number_gt                  5024
_reflns_number_total               7595
_reflns_threshold_expression       I>=2u(I)
_refine_diff_density_max           1.6950
_refine_diff_density_min           -0.8025
_refine_diff_density_rms           0.0870
_refine_ls_d_res_high              0.7097
_refine_ls_d_res_low               6.6752
_refine_ls_goodness_of_fit_ref     0.9034
_refine_ls_hydrogen_treatment      constr
_refine_ls_matrix_type             full
_refine_ls_number_constraints      0
_refine_ls_number_parameters       495
_refine_ls_number_reflns           7595
_refine_ls_number_restraints       29
_refine_ls_R_factor_all            0.0548
_refine_ls_R_factor_gt             0.0361
_refine_ls_restrained_S_all        0.9015
_refine_ls_shift/su_max            -0.4325
_refine_ls_shift/su_mean           0.0077
_refine_ls_structure_factor_coef   Fsqd
_refine_ls_weighting_details      
 'w=1/[\s^2^(Fo^2^)]'
_refine_ls_weighting_scheme        calc
_refine_ls_wR_factor_gt            0.0871
_refine_ls_wR_factor_ref           0.0924

loop_
  _atom_site_label
  _atom_site_type_symbol
  _atom_site_fract_x
  _atom_site_fract_y
  _atom_site_fract_z
  _atom_site_U_iso_or_equiv
  _atom_site_adp_type
  _atom_site_occupancy
  _atom_site_refinement_flags_posn
  _atom_site_refinement_flags_adp
 Sr Sr 0.000(18) 0.000(14) 0.000(2) 0.04(8) Uiso 1.000000 . .

"""
    fcf_header = """
data_Sr
_computing_structure_refinement   'olex2.refine 1.5-dev (Bourhis et al., 2015)'
_shelx_refln_list_code            6
loop_
  _space_group_symop_id
  _space_group_symop_operation_xyz
  1  x,y,z

_space_group_crystal_system       triclinic
_space_group_IT_number            1
_space_group_name_H-M_alt         'P 1'
_space_group_name_Hall            'P 1'
_symmetry_space_group_name_H-M    'P 1'
_symmetry_space_group_name_Hall   'P 1'
_symmetry_Int_Tables_number       1
_cell_length_a                    10.000(2)
_cell_length_b                    10.000(3)
_cell_length_c                    10.000(2)
_cell_angle_alpha                 90.00(2)
_cell_angle_beta                  90.00(2)
_cell_angle_gamma                 90.00(2)
_cell_volume                      1000.0(5)
loop_
  _refln_index_h
  _refln_index_k
  _refln_index_l
  _refln_F_squared_meas
  _refln_F_squared_sigma
  _refln_F_calc
  _refln_phase_calc"""
    
    h = range(-30,31)
    k = range(-30,31)
    l = range(-30,31)
    #h = range(-1,2)
    #k = range(-1,2)
    #l = range(-1,2)
    listy = [h,k,l]
    import itertools
    from matrix_coefficients_v2 import prepare_M
    indices_temp = list(itertools.product(*listy))
    indices = []
    for i in indices_temp:
      if i[0] == i[1] == i[2] == 0: continue
      if stol_from_hkl(i[0],i[1],i[2]) > 1.0: continue
      indices.append(i)
    
    read = True
    if read == False:
      M_matrizes = []
      for shell in range(0,6):
        M_matrizes.append([])
      for shell in range(1,6):
        for n_0 in range(0,4):
          M_matrizes[shell].append([])
      for shell in range(1,6):
        for n_0 in range(1,4):
          for l in range(l_max+2):
            M_matrizes[shell][n_0].append([])
      for shell in range(1,6):
        for n_0 in range(1,4):
          for l in range(l_max+2):
            for p in range(p_max+6):
              z = get_z_temp(38,n_0,shell,nu)
              if z == -100: continue
              M_matrizes[shell][n_0][l].append(prepare_M(p,l,z,n_0))
      precalc_0  = NF(38,nu, 0,4,0,False,M_matrizes)
      x = np.linspace(0,math.pi,128)
      with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        precalc  = pool.starmap(NF, zip(repeat(38), repeat(nu), x, repeat(4), repeat(2), repeat(False), repeat(M_matrizes)))
      
      precalc *= 2
      preacalc_arrays = (x,precalc,precalc_0)
      with open("t0.npy","wb") as f:
        np.save(f,x)
      with open("fp_fdp.npy","wb") as f:
        np.save(f,precalc)
    else:
      x = np.load("t0.npy")
      angle = x / math.pi * 180
      precalc = np.load("fp_fdp.npy")
      precalc *= math.pi
      im = precalc[0:128].imag
      re = precalc[0:128].real
      fig, axes = plt.subplots(1,1)
      axes.plot(angle,re,"b",label="f'")# type: ignore
      axes.plot(angle,im,"k",label="f''")# type: ignore
      axes.set_ylabel(r"$f$ /$e$")
      axes.set_xlabel(r"$2\theta$")
      axes.axvline(x=50,linestyle='dashed',color='gray',label='Min. IUCr')# type: ignore
      axes.legend()# type: ignore
      
      fig.set_figheight(inches(17))
      fig.set_figwidth(inches(22))
      
      fig.savefig('Sr_integrated.png',dpi=300,bbox_inches='tight')
      precalc_0 = -0.5155468923122816+1.6399873146939559j
      preacalc_arrays = (x,precalc[0:128],precalc_0)
      
    print("We have f' and f''!")
    print("Calculating number of dispersion electrons...")
    #n_disp = NF_DL(38,nu,0.0,5,2,M_matrizes)
    n_disp = 0.0
    #print("Number of disp els: ", n_disp)
    #print(NF(38,nu,0,5,4,False,M_matrizes))
    with open("Sr.cif","wt") as ciffy:
      ciffy.write(cif)
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
      print("starting calc!")
      res  = pool.starmap(F_calc, zip(indices,repeat(n_disp),repeat(preacalc_arrays)))
      with open("Sr.fcf","wt") as fcfy:
        fcfy.write(fcf_header+'\n')
        fcfy.write('\n'.join(' '.join(map(str, r)) for r in res))
    
    print("Done with calculation!")
    exit()

