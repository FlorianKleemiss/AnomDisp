import math
import scipy.integrate as integrate
from constants_and_atomic_properties import *
# distutils: language=Py3

parallel = True

def calc_stuff(position: int, number: int, n_disp_K: float, n_disp_L:float, n_disp_M:float, axis: str, edge: float) -> "list[float]":
  energy = 0
  x = 0
  if axis == "relative":
    energy = 1/(position * 0.001) * edge
    rel_pos = position * 0.001
    x = (rel_pos) #For relative to K-edge
  elif axis == "keV":
    energy = position
    x = (energy/1000) #for absolute energy spectra in keV
  elif axis == "Angstrom" or axis == "fixed":
    energy = 1E10*h*speed_of_light/(position*0.0001)
    wavelength = (position*0.0001)
    x = (wavelength) #for wavelength in A
  if energy == edge:
    return [0,0,0,0,0,0,0,0,0]
  else:
    result_s, result_s_imag = sugiura_k(number,energy)
    y_s = 2*result_s - n_disp_K
    y_s_imag = 2*result_s_imag
    result_real, result_imag = l_with_imag(number,energy)
    y_l_real = 2*result_real-n_disp_L
    y_l_imag = 2*result_imag
    result_real_M, result_imag_M = m_with_imag(number,energy)
    y_m_real = 2*result_real_M-n_disp_M
    y_m_imag = 2*result_imag_M

    y_d = 2*result_s+2*result_real+2*result_real_M-(n_disp_K+n_disp_L+n_disp_M)
    y_d_imag = 2*result_s_imag+2*result_imag+2*result_imag_M
    if axis == "fixed":
      print("{:5.3f} ".format(x)+"{:8e} ".format(2*result_s-n_disp_K)+"{:8e} ".format(2*result_real-n_disp_L)+"{:8e} ".format(2*result_s_imag)+"{:8e} ".format(2*result_imag)+"{:8e} ".format(math.sqrt(pow(i+2*result_s-n_disp_K+2*result_real-n_disp_L,2)+pow(2*result_s_imag+2*result_imag,2))) )
    return [x,y_s,y_s_imag,y_l_real,y_l_imag,y_m_real,y_m_imag,y_d,y_d_imag]

def calc_hoenllike(energy_in_ev: float, Z: int) -> "list[float]":
  junk, result_s_imag = sugiura_k(Z,energy_in_ev)
  y_s_imag = 2*result_s_imag
  junk, result_imag = l_with_imag(Z,energy_in_ev)
  y_l_imag = 2*result_imag
  junk, result_imag_M = m_with_imag(Z,energy_in_ev)
  y_m_imag = 2*result_imag_M

  y_d_imag = 2*result_s_imag+2*result_imag+2*result_imag_M
  return [y_s_imag,y_l_imag,y_m_imag,y_d_imag]

def integrand_for_disp_els_K(n_j:float,n_0:float) -> float:
  z=n_j/n_0
  res = f_s_1_hoenl(z) * 2/n_0/z
  return res
  part1 = 128/n_0/(3*pow(z,4))
  return part1*sugiura_exps(z,1)

def integrand_for_disp_els_L1(n_j:float,n_0:float) -> float:
  z=n_j/n_0
  res = f_s_1_EM(z) * 2/n_0/z
  return res
  part1 = 1024/n_0*(z+3)/(3*pow(z,5))
  return part1*sugiura_exps(z,2)

def integrand_for_disp_els_L2_3(n_j:float,n_0:float) -> float:
  z=n_j/n_0
  res = f_p_1_EM(z) * 2/n_0/z
  return res
  part1 = 1024/n_0*(z+8.0/3.0)/(3*pow(z,6))
  return part1*sugiura_exps(z,2)

def integrand_for_disp_els_M1(n_j:float,n_0:float) -> float:
  z=n_j/n_0
  res = f_s_1_WA(z) * 2/n_0/z
  return res
  part1 = 64*2/n_0*pow(3*z+4,2)*(z+8)/(pow(z,7))
  return part1*sugiura_exps(z,3)

def integrand_for_disp_els_M_2_3(n_j:float,n_0:float) -> float:
  z=n_j/n_0
  res = f_p_1_WA(z) * 2 /n_0 /z
  return res
  part1 = 512*2/n_0*(3*z*z+26*z+28)/(pow(z,7))
  return part1*sugiura_exps(z,3)

def integrand_for_disp_els_M_4_5(n_j:float,n_0:float) -> float:
  z=n_j/n_0
  res = f_d_1_WA(z) * 2 /n_0 /z
  return res
  part1 = 1024/5.0*2/n_0*(5*z*z+46*z+48)/(pow(z,8))
  return part1*sugiura_exps(z,3)

def real_general_integration(n_j:float,chi_j:float,nu:float,n_0:float,oscillator_density_function) -> float:
  res = (pow(n_j,2)-pow(nu,2))/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*oscillator_density_function(n_j,n_0)

def imag_general_integration(n_j:float,chi_j:float,nu:float,n_0:float,oscillator_density_function) -> float:
  res = -nu*chi_j/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*oscillator_density_function(n_j,n_0)

def perform_imag_integration(nu: float, nu_edge: float, x_j:float, n_0:float, func: Callable) -> float:
  epsilon2 = 0.1
  lower_limit = nu_edge
  if (1+epsilon2)*nu < nu_edge:
    upper_limit = 1000*nu_edge
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral = integrate.quad(imag_general_integration,
                                 lower_limit,
                                 upper_limit,
                                 args=(x_j,nu,n_0,func),
                                 points=nu,
                                 limit=200000,
                                 epsabs=1E-55,
                                 epsrel=1E-10)
  if imag_integral[1] > 0.5*abs(imag_integral[0]):
    print("!!! inaccurate IMAG integral at nu_edge/nu = " + "{:6.4f}".format(nu_edge/nu) + " ratio:" + str(imag_integral[1]/abs(imag_integral[0])) + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))
  return imag_integral[0]

def perform_real_integration(nu: float, nu_edge:float, x_j:float, n_0: float,func: Callable) -> float:
  inte1 = integrate.quad(real_general_integration,
                         nu_edge,
                         nu,
                         args=(x_j,nu,n_0,func),
                         points=nu,
                         limit=200000,
                         epsabs=1E-50,
                         epsrel=1E-10)
  inte2 = integrate.quad(real_general_integration,
                         nu,
                         100*nu,
                         args=(x_j,nu,n_0,func),
                         points=nu,
                         limit=200000,
                         epsabs=1E-50,
                         epsrel=1E-10)
  if inte1[1] > 0.5*abs(inte1[0]):
    print("!!! inaccurate REAL1 integral at nu_edge/nu = " + "{:8.4f}".format(nu_edge/nu) + " " + str(inte1) + str(inte1[1]/abs(inte1[0])))
  if inte2[1] > 0.5*abs(inte2[0]):
    print("!!! inaccurate REAL2 integral at nu_edge/nu = " + "{:8.4f}".format(nu_edge/nu) + " " + str(inte2) + str(inte2[1]/abs(inte2[0])))
  integral = inte1[0] + inte2[0]
  return integral

@overload
def sugiura_k(Z: int = -1, ener:float =8047.8, disp_only: Literal[False] = ...) -> "list[float]": ...

@overload
def sugiura_k(Z: int = -1, ener:float =8047.8, disp_only: Literal[True] = ...) -> float: ...

def sugiura_k(Z: int = -1, ener:float =8047.8, disp_only: bool=False) -> Union[float,"list[float]"]:
  # Z is integer number of 
  if Z == -1:
    raise ValueError("Z MUST BE POSITIVE INTEGER!")
  e_ion = get_ionization_energy_1s(Z)
  nu_k = e_ion / h
  omdelta_K = one_minus_delta_edge(Z,1,0,0) #21a in Hoenl
  n_0 = nu_k/omdelta_K #22 in Hoenl
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  x_j = get_line_width_K(Z)/h
  #integral_test = integrate.quad(integrand_damped_abs,nu_k,50*nu,args=(x_j,nu,n_0),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  #integrating directly as one intergal with the singulary position as point does not work.
  
  if disp_only == True:
    n_disp_el = integrate.quad(integrand_for_disp_els_K,nu_k,20000*nu_k,args=(n_0),limit=200000,epsabs=1E-60)[0]

    return 2*n_disp_el

  #REAL PART
  integral = perform_real_integration(nu,nu_k, x_j, n_0, integrand_for_disp_els_K)
  #IMAG PART
  imag_integral = perform_imag_integration(nu,nu_k, x_j, n_0, integrand_for_disp_els_K)

  alpha_K_sugiura_damp = complex(integral,imag_integral)
  ausdruck = alpha_K_sugiura_damp/(1-4.0*math.pi/3*alpha_K_sugiura_damp*prefactor)
  factor = pow(nu,2)
  return [-ausdruck.real*factor, -ausdruck.imag*factor]

@overload
def l_with_imag(Z: int = -1, ener:float =8047.8, disp_only: Literal[False] = ...) -> "list[float]": ...

@overload
def l_with_imag(Z: int = -1, ener:float =8047.8, disp_only: Literal[True] = ...) -> float: ...

def l_with_imag(Z: int = -1, ener:float =8047.8, disp_only: bool=False) -> Union[float,"list[float]"]:
  # Z is integer number of 
  if Z == -1:
    raise ValueError("Z MUST BE POSITIVE INTEGER!")
  nu_l1 = get_ionization_energy_2s(Z)   /h
  nu_l2 = get_ionization_energy_2p1_2(Z)/h
  nu_l3 = get_ionization_energy_2p3_2(Z)/h
  omdelta_l1 = one_minus_delta_edge(Z,2,0,0)
  omdelta_l2 = one_minus_delta_edge(Z,2,1,0.5)
  omdelta_l3 = one_minus_delta_edge(Z,2,1,1.5)
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  x_j1 = get_line_width_Ls(Z)/h
  x_j2 = get_line_width_Lp_1_2(Z)/h
  x_j3 = get_line_width_Lp_3_2(Z)/h
  n_0_1 = nu_l1/omdelta_l1 #22 in Hoenl
  n_0_2 = nu_l2/omdelta_l2 #22 in Hoenl
  n_0_3 = nu_l3/omdelta_l3 #22 in Hoenl
  
  if disp_only == True:
    n_disp_2s = integrate.quad(integrand_for_disp_els_L1,nu_l1,80000*nu_l1,args=(n_0_1),limit=200000,epsabs=1E-12)[0]
    n_disp_2p1_2 = integrate.quad(integrand_for_disp_els_L2_3,nu_l2,100000*nu_l2,args=(n_0_2),limit=200000,epsabs=1E-12)[0]
    n_disp_2p3_2 = integrate.quad(integrand_for_disp_els_L2_3,nu_l3,100000*nu_l3,args=(n_0_3),limit=200000,epsabs=1E-12)[0]

    return 2*n_disp_2s+2*n_disp_2p1_2+4*n_disp_2p3_2

  #REAL PART
  #When telling hte quad alorithm that nu is a "position to be carefull with" the results don't require epsilon
  integral_1 = perform_real_integration(nu,nu_l1, x_j1, n_0_1, integrand_for_disp_els_L1)
  integral_2 = perform_real_integration(nu,nu_l2, x_j2, n_0_2, integrand_for_disp_els_L2_3)
  integral_3 = perform_real_integration(nu,nu_l3, x_j3, n_0_3, integrand_for_disp_els_L2_3)

  #IMAG PART
  imag_integral_1 = perform_imag_integration(nu,nu_l1, x_j1, n_0_1, integrand_for_disp_els_L1)
  imag_integral_2 = perform_imag_integration(nu,nu_l2, x_j2, n_0_2, integrand_for_disp_els_L2_3)
  imag_integral_3 = perform_imag_integration(nu,nu_l3, x_j3, n_0_3, integrand_for_disp_els_L2_3)
  
  alpha_L_sugiura_damp = complex(integral_1+integral_2+integral_3*2,imag_integral_1+imag_integral_2+2*imag_integral_3)
  ausdruck = (alpha_L_sugiura_damp)/(1-4.0*math.pi/3*alpha_L_sugiura_damp*prefactor)
  factor = pow(nu,2)
  return [-ausdruck.real*factor, -ausdruck.imag*factor]

@overload
def m_with_imag(Z: int = -1, ener:float =8047.8, disp_only: Literal[False] = ...) -> "list[float]": ...

@overload
def m_with_imag(Z: int = -1, ener:float =8047.8, disp_only: Literal[True] = ...) -> float: ...

def m_with_imag(Z: int=-1, ener:float=8047.8, disp_only:bool =False) -> Union[float,"list[float]"]:
  # Z is integer number of 
  if Z == -1:
    raise ValueError("Z MUST BE POSITIVE INTEGER!")
  nu_m1 = get_ionization_energy_3s(Z)    / h
  nu_m2 = get_ionization_energy_3p_1_2(Z)/ h
  nu_m3 = get_ionization_energy_3p_3_2(Z)/ h
  nu_m4 = get_ionization_energy_3d_3_2(Z)/ h
  nu_m5 = get_ionization_energy_3d_5_2(Z)/ h
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  omdelta_m1 = one_minus_delta_edge(Z,3,0,0)
  omdelta_m2 = one_minus_delta_edge(Z,3,1,0.5)
  omdelta_m3 = one_minus_delta_edge(Z,3,1,1.5)
  omdelta_m4 = one_minus_delta_edge(Z,3,2,1.5)
  omdelta_m5 = one_minus_delta_edge(Z,3,2,2.5)
  n_0_1 = nu_m1/omdelta_m1 #22 in Hoenl
  n_0_2 = nu_m2/omdelta_m2 #22 in Hoenl
  n_0_3 = nu_m3/omdelta_m3 #22 in Hoenl
  n_0_4 = nu_m4/omdelta_m4 #22 in Hoenl
  n_0_5 = nu_m5/omdelta_m5 #22 in Hoenl

  if disp_only == True:
    n_disp_3s    = integrate.quad(integrand_for_disp_els_M1   ,nu_m1,10000*nu_m1,args=(n_0_1),limit=200000,epsabs=1E-60)[0]
    n_disp_3p1_2 = integrate.quad(integrand_for_disp_els_M_2_3,nu_m2,10000*nu_m2,args=(n_0_2),limit=200000,epsabs=1E-60)[0]
    n_disp_3p3_2 = integrate.quad(integrand_for_disp_els_M_2_3,nu_m3,10000*nu_m3,args=(n_0_3),limit=200000,epsabs=1E-60)[0]
    n_disp_3d3_2 = integrate.quad(integrand_for_disp_els_M_4_5,nu_m4,10000*nu_m4,args=(n_0_4),limit=200000,epsabs=1E-60)[0]
    n_disp_3d5_2 = integrate.quad(integrand_for_disp_els_M_4_5,nu_m5,10000*nu_m5,args=(n_0_5),limit=200000,epsabs=1E-60)[0]

    return 2*n_disp_3s+2*n_disp_3p1_2+4*n_disp_3p3_2+4*n_disp_3d3_2+6*n_disp_3d5_2

  x_j1 = get_line_width_Ls(Z)*(get_ionization_energy_3s(Z)/get_ionization_energy_2s(Z))/h
  x_j2 = get_line_width_Lp_1_2(Z)*(get_ionization_energy_3p_1_2(Z)/get_ionization_energy_2p1_2(Z))/h
  x_j3 = get_line_width_Lp_3_2(Z)*(get_ionization_energy_3p_3_2(Z)/get_ionization_energy_2p3_2(Z))/h
  x_j4 = get_line_width_Lp_3_2(Z)*(get_ionization_energy_3d_3_2(Z)/get_ionization_energy_2p3_2(Z))/h
  x_j5 = x_j4
  
  #REAL PART
  integral_1 = perform_real_integration(nu,nu_m1, x_j1, n_0_1, integrand_for_disp_els_M1)
  integral_2 = perform_real_integration(nu,nu_m2, x_j2, n_0_2, integrand_for_disp_els_M_2_3)
  integral_3 = perform_real_integration(nu,nu_m3, x_j3, n_0_3, integrand_for_disp_els_M_2_3)
  integral_4 = perform_real_integration(nu,nu_m4, x_j4, n_0_4, integrand_for_disp_els_M_4_5)
  integral_5 = perform_real_integration(nu,nu_m5, x_j5, n_0_5, integrand_for_disp_els_M_4_5)

  #IMAG PART
  imag_integral_1 = perform_imag_integration(nu,nu_m1, x_j1, n_0_1, integrand_for_disp_els_M1)
  imag_integral_2 = perform_imag_integration(nu,nu_m2, x_j2, n_0_2, integrand_for_disp_els_M_2_3)
  imag_integral_3 = perform_imag_integration(nu,nu_m3, x_j3, n_0_3, integrand_for_disp_els_M_2_3)
  imag_integral_4 = perform_imag_integration(nu,nu_m4, x_j4, n_0_4, integrand_for_disp_els_M_4_5)
  imag_integral_5 = perform_imag_integration(nu,nu_m5, x_j5, n_0_5, integrand_for_disp_els_M_4_5)
  
  alpha_M_sugiura_damp = complex(integral_1     +integral_2+     2*integral_3+     3*integral_4+     2*integral_5,\
                                 imag_integral_1+imag_integral_2+2*imag_integral_3+3*imag_integral_4+2*imag_integral_5)
  ausdruck = (alpha_M_sugiura_damp)/(1-4.0*math.pi/3*alpha_M_sugiura_damp*prefactor)
  factor = pow(nu,2)
  return [-ausdruck.real*factor, -ausdruck.imag*factor]

if __name__ == '__main__':
  import matplotlib.pyplot as plt
  from matplotlib import ticker
  fig, axs = plt.subplots(4,2,sharex=True)
  stop = False
  Zs = [42,52,92]
  edges = []
  L_edges = []
  M_edges = []
  axis = "Angstrom"
  #while stop == False:
  #  new_Z = int(input("Please enter an element number (0 will terminate input and start calculation): "))
  #  if new_Z != 0:
  #    Zs.append(new_Z)
  #  else:
  #    stop = True
  for i in Zs:
    edge = get_ionization_energy_1s(i)
    if axis == "Angstrom":
      edges.append(h*speed_of_light*1E10/edge)
      L_edges.append(h*speed_of_light*1E10/get_ionization_energy_2s(i))
      L_edges.append(h*speed_of_light*1E10/get_ionization_energy_2p1_2(i))
      L_edges.append(h*speed_of_light*1E10/get_ionization_energy_2p3_2(i))
      M_edges.append(h*speed_of_light*1E10/get_ionization_energy_3s(i))
      M_edges.append(h*speed_of_light*1E10/get_ionization_energy_3p_1_2(i))
      M_edges.append(h*speed_of_light*1E10/get_ionization_energy_3p_3_2(i))
      M_edges.append(h*speed_of_light*1E10/get_ionization_energy_3d_3_2(i))
      M_edges.append(h*speed_of_light*1E10/get_ionization_energy_3d_5_2(i))
    elif axis=="keV":
      edges.append(edge/1000)
      L_edges.append(get_ionization_energy_2s(i)/1000)
      L_edges.append(get_ionization_energy_2p1_2(i)/1000)
      L_edges.append(get_ionization_energy_2p3_2(i)/1000)
      M_edges.append(get_ionization_energy_3s(i)/1000)
      M_edges.append(get_ionization_energy_3p_1_2(i)/1000)
      M_edges.append(get_ionization_energy_3p_3_2(i)/1000)
      M_edges.append(get_ionization_energy_3d_3_2(i)/1000)
      M_edges.append(get_ionization_energy_3d_5_2(i)/1000)
    elif axis == "relative":
      edges.append(1.0)
      L_edges.append(get_ionization_energy_2s(i)/get_ionization_energy_1s(i))
      L_edges.append(get_ionization_energy_2p1_2(i)/get_ionization_energy_1s(i))
      L_edges.append(get_ionization_energy_2p3_2(i)/get_ionization_energy_1s(i))
    x = []
    y_s_real = []
    y_s_imag = []
    y_d_real = []
    y_d_imag = []
    y_l_real = []
    y_l_imag = []
    y_m_real = []
    y_m_imag = []
    steps = []
    if axis == "relative":
      minimal = 5
      maximal = 2000
      nr_steps = 5
      for step in range(int(minimal),int(maximal),int(nr_steps)):
        steps.append(step)
    elif axis == "keV":
      minimal = 100
      maximal = 200000
      nr_steps = 10
      for step in range(int(minimal),int(maximal),int(nr_steps)):
        steps.append(step)
    elif axis == "Angstrom":
      minimal = 1000
      maximal = 30000
      nr_steps = 50
      for step in range(int(minimal),int(maximal),int(nr_steps)):
        steps.append(step)
    elif axis == "fixed":
      if i == 52:
        steps = [78,156,233,272,311,350,369,381,388.9,389.1,397,408,428,467,545,558,622,700,778,1167,1200,1556,1700,1945,2100,2400,2450,2480,2504.9,2505.1,2520,2530,2550,2600,2650,2681.9,2682.1,2700,2720,2723,2750,2770,2800,2849.9,2850.1,2900,3000,3500,4000,100000]
      elif i == 74:
        steps = [80,120,150,160,170,175,177.9,178.1,180,185,200,250,400,600,800,900,980,1010,1021.9,1022.1,1030,1040,1050,1065,1071.9,1072.1,1080,1120,1160,1205,1212.9,1213.1,1220,1300,1500,2000]
    n_disp_K = sugiura_k(i,disp_only=True)
    n_disp_L = l_with_imag(i,disp_only=True)
    n_disp_M = m_with_imag(i,disp_only=True)
    print(" For Z="+str(i)+" el disps: "+str(n_disp_K) + " " + str(n_disp_L) + " " + str(n_disp_M))
    if parallel == True:
      import multiprocessing
      import time
      from itertools import repeat
      time_1 = time.perf_counter()
      with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        res = pool.starmap(calc_stuff, zip(steps, repeat(i), repeat(n_disp_K), repeat(n_disp_L), repeat(n_disp_M), repeat(axis), repeat(edge)))
        for r in res:
          x.append(r[0])
          y_s_real.append(r[1])
          y_s_imag.append(r[2])
          y_l_real.append(r[3])
          y_l_imag.append(r[4])
          y_m_real.append(r[5])
          y_m_imag.append(r[6])
          y_d_real.append(r[7])
          y_d_imag.append(r[8])
      print("Timing: " + str(time.perf_counter()-time_1))
    
    else:
      import time
      time_1 = time.perf_counter()
      for step in steps:
        result = calc_stuff(step, i, n_disp_K, n_disp_L, n_disp_M, axis, edge)
        x.append(result[0])
        y_s_real.append(result[1])
        y_s_imag.append(result[2])
        y_l_real.append(result[3])
        y_l_imag.append(result[4])
        y_m_real.append(result[5])
        y_m_imag.append(result[6])
        y_d_real.append(result[7])
        y_d_imag.append(result[8])
      print("Timing: " + str(time.perf_counter()-time_1))
      
    axs[0,0].plot(x,y_s_real,':',label="%s K_edge"%elements[i]) # type: ignore
    axs[0,1].plot(x,y_s_imag,':',label="%s K_edge"%elements[i]) # type: ignore

    axs[1,0].plot(x,y_l_real,':',label="%s L_edges"%elements[i]) # type: ignore
    axs[1,1].plot(x,y_l_imag,':',label="%s L_edges"%elements[i]) # type: ignore

    axs[2,0].plot(x,y_m_real,':',label="%s M_edges"%elements[i]) # type: ignore
    axs[2,1].plot(x,y_m_imag,':',label="%s M_edges"%elements[i]) # type: ignore

    axs[3,0].plot(x,y_d_real,':',label="%s"%elements[i]) # type: ignore
    axs[3,1].plot(x,y_d_imag,':',label="%s"%elements[i]) # type: ignore
    axs[3,0].axhline(y=0, color="gray", linestyle="--") # type: ignore

  axs[0, 0].set_title('K Real') # type: ignore
  axs[1, 0].set_title('L Real') # type: ignore
  axs[2, 0].set_title('M Real') # type: ignore
  axs[3, 0].set_title('Complete Real') # type: ignore
  axs[0, 1].set_title('K Imag') # type: ignore
  axs[1, 1].set_title('L Imag') # type: ignore
  axs[2, 1].set_title('M Imag') # type: ignore
  axs[3, 1].set_title('Complete Imag') # type: ignore
  axs[0, 1].legend(loc='upper right') # type: ignore 
  axs[1, 1].legend(loc='upper right')# type: ignore
  axs[2, 1].legend(loc='upper right')# type: ignore
  axs[3, 1].legend(loc='upper right')# type: ignore
  for j in range(len(edges)):
      axs[0,0].axvline(x=edges[j], color="gray", linestyle=":")# type: ignore
      axs[0,1].axvline(x=edges[j], color="gray", linestyle=":")# type: ignore
      axs[3,0].axvline(x=edges[j], color="gray", linestyle=":")# type: ignore
      axs[3,1].axvline(x=edges[j], color="gray", linestyle=":")# type: ignore
  for j in range(len(L_edges)):
      axs[1,0].axvline(x=L_edges[j], color="gray", linestyle=":")# type: ignore
      axs[1,1].axvline(x=L_edges[j], color="gray", linestyle=":")# type: ignore
      axs[3,0].axvline(x=L_edges[j], color="gray", linestyle=":")# type: ignore
      axs[3,1].axvline(x=L_edges[j], color="gray", linestyle=":")# type: ignore
  for j in range(len(M_edges)):
      axs[2,0].axvline(x=M_edges[j], color="gray", linestyle=":")# type: ignore
      axs[2,1].axvline(x=M_edges[j], color="gray", linestyle=":")# type: ignore
      axs[3,0].axvline(x=M_edges[j], color="gray", linestyle=":")# type: ignore
      axs[3,1].axvline(x=M_edges[j], color="gray", linestyle=":")# type: ignore
  if axis == "Angstrom" or axis == "fixed":
    axs[3, 0].set(xlabel=r'$\lambda\ /\AA$')# type: ignore
    axs[3, 1].set(xlabel=r'$\lambda\ /\AA$')# type: ignore
  elif axis == "keV":
    axs[3, 0].set(xlabel=r'$E /keV$')# type: ignore
    axs[3, 1].set(xlabel=r'$E /keV$')# type: ignore
  elif axis == "relative":
    axs[3, 0].set(xlabel=r'$\frac{\lambda}{\lambda_K}$')# type: ignore
    axs[3, 1].set(xlabel=r'$\frac{\lambda}{\lambda_K}$')# type: ignore
  axs[0, 0].set(ylabel=r'$\frac{1-n_0}{\lambda^2} * \frac{2 \pi m c^2}{N e^2}$')# type: ignore
  axs[1, 0].set(ylabel=r'$\frac{1-n_0}{\lambda^2} * \frac{2 \pi m c^2}{N e^2}$')# type: ignore
  axs[2, 0].set(ylabel=r'$\frac{1-n_0}{\lambda^2} * \frac{2 \pi m c^2}{N e^2}$')# type: ignore
  axs[3, 0].set(ylabel=r'$\Delta f\ /e$')# type: ignore
  for ax in fig.get_axes():
    if axis == "Angstrom" or axis == "fixed":
      #ax.set_xscale('log')
      ax.set_xlim(3.0,0.05)
      ax.set_xticks([0.1,0.2,0.3,0.5,0.7,1.0,1.5,2.0,3.0])
      ax.get_xaxis().get_major_formatter().labelOnlyBase = False
      ax.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    elif axis == "relative":
      ax.set_xlim(0.0,2.0) #for relative spectra
    elif axis == "keV":
      ax.set_xlim(0.1,200) #for energy spectra in keV
  plt.show()

  print("Done")

def sugiura_k_purely_imag(Z:int=-1, ener:float=8047.8) -> float:
  # Z is integer number of 
  if Z <= 0:
    raise ValueError("Z MUST BE Positive INTEGER!")
  e_ion = get_ionization_energy_1s(Z)
  nu_k = e_ion / h
  omdelta_K = one_minus_delta_edge(Z,1,0,0) #21a in Hoenl
  #delta_K = 0
  n_0 = nu_k/omdelta_K #22 in Hoenl
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  x_j = get_line_width_K(Z)/h

  #IMAG PART
  imag_integral = perform_imag_integration(nu,nu_k,x_j,n_0,integrand_for_disp_els_K)
  
  return -imag_integral