import math
#from matplotlib.pyplot import get
import numpy
import scipy.integrate as integrate
import scipy.special as special
import matrix_coeficients
from constants_and_atomic_properties import *
#from scipy.optimize._lsq.least_squares import arctan

parallel = True

def calc_stuff(position, number, n_disp_K, n_disp_L, n_disp_M, axis, edge):
  if axis == "relative":
    energy = 1/(position * 0.001) * edge
    rel_pos = position * 0.001
  elif axis == "keV":
    energy = position
  elif axis == "Angstrom" or axis == "fixed":
    energy = 1E10*h*speed_of_light/(position*0.0001)
    wavelength = (position*0.0001)
  if energy == edge:
    return
  else:
    result_s, result_s_imag = sugiura_k(number,energy)
    if axis == "relative":
      x = (rel_pos) #For relative to K-edge
    elif axis == "keV":
      x = (energy/1000) #for absolute energy spectra in keV
    elif axis == "Angstrom" or axis == "fixed":
      x = (wavelength) #for wavelength in A
    y_s = 2*result_s - n_disp_K
    y_s_imag = 2*result_s_imag
    result_real, result_imag = l_with_imag(number,energy)
    if result_imag == math.nan:
      print("problem!")
    y_l_real = 2*result_real-n_disp_L
    y_l_imag = 2*result_imag
    result_real_M, result_imag_M = m_with_imag(number,energy)
    y_m_real = 2*result_real_M-n_disp_M
    y_m_imag = 2*result_imag_M

    y_d = 2*result_s+2*result_real+2*result_real_M-(n_disp_K+n_disp_L+n_disp_M)
    y_d_imag = 2*result_s_imag+2*result_imag+2*result_imag_M
    if axis == "fixed":
      print("{:5.3f} ".format(x[len(x)-1])+"{:8e} ".format(2*result_s-n_disp_K)+"{:8e} ".format(2*result_real-n_disp_L)+"{:8e} ".format(2*result_s_imag)+"{:8e} ".format(2*result_imag)+"{:8e} ".format(math.sqrt(pow(i+2*result_s-n_disp_K+2*result_real-n_disp_L,2)+pow(2*result_s_imag+2*result_imag,2))) )
    return x,y_s,y_s_imag,y_l_real,y_l_imag,y_m_real,y_m_imag,y_d,y_d_imag

def integrand_for_disp_els_K(n_j,n_0):
  z=n_j/n_0
  part1 = 128/n_0/(3*pow(z,4))
  z_1=z-1
  if(n_j < n_0):
    sqrt_var = numpy.sqrt(abs(z_1))
    temp = numpy.arctan(sqrt_var)/sqrt_var - 2*1.0/3.0*z_1*special.hyp2f1(0.75,1,1.75,z_1*z_1)
    part2 = numpy.exp(-4*temp)
    part3 = 1-numpy.exp(-2*math.pi/sqrt_var)
  else:
    sqrt_var = numpy.sqrt(z_1)
    part2 = numpy.exp(-4/sqrt_var*numpy.arctan(sqrt_var))
    part3 = 1-numpy.exp(-2*math.pi/sqrt_var)
  return part1*part2/part3

def integrand_for_disp_els_L1(n_j,n_0):
  z=n_j/n_0
  part1 = 1024/n_0*(z+3)/(3*pow(z,5))
  z_1=z-1
  if(n_j < n_0):
    sqrt_var = numpy.sqrt(abs(z_1))
    temp = numpy.arctan(sqrt_var)/sqrt_var - 2*1.0/3.0*z_1*special.hyp2f1(0.75,1,1.75,z_1*z_1)
    part2 = numpy.exp(-8*temp)
    part3 = 1-numpy.exp(-4*math.pi/sqrt_var)
  else:
    sqrt_var = numpy.sqrt(z_1)
    part2 = numpy.exp(-8/sqrt_var*numpy.arctan(sqrt_var))
    part3 = 1-numpy.exp(-4*math.pi/sqrt_var)
  return part1*part2/part3

def integrand_for_disp_els_L2_3(n_j,n_0):
  z=n_j/n_0
  part1 = 1024/n_0*(z+8.0/3.0)/(3*pow(z,6))
  z_1=z-1
  if(n_j < n_0):
    sqrt_var = numpy.sqrt(abs(z_1))
    temp = numpy.arctan(sqrt_var)/sqrt_var - 2*1.0/3.0*z_1*special.hyp2f1(0.75,1,1.75,z_1*z_1)
    part2 = numpy.exp(-8*temp)
    part3 = 1-numpy.exp(-4*math.pi/sqrt_var)
  else:
    sqrt_var = numpy.sqrt(z_1)
    part2 = numpy.exp(-8/sqrt_var*numpy.arctan(sqrt_var))
    part3 = 1-numpy.exp(-4*math.pi/sqrt_var)
  return part1*part2/part3

def integrand_for_disp_els_M1(n_j,n_0):
  z=n_j/n_0
  part1 = 64*2/n_0*pow(3*z+4,2)*(z+8)/(pow(z,7))
  z_1=z-1
  if(n_j < n_0):
    sqrt_var = numpy.sqrt(abs(z_1))
    temp = numpy.arctan(sqrt_var)/sqrt_var - 2*1/3.0*z_1*special.hyp2f1(0.75,1,1.75,z_1*z_1)
    part2 = numpy.exp(-12*temp)
    part3 = 1-numpy.exp(-6*math.pi/sqrt_var)
  else:
    sqrt_var = numpy.sqrt(z_1)
    part2 = numpy.exp(-16/sqrt_var*numpy.arctan(sqrt_var))
    part3 = 1-numpy.exp(-6*math.pi/sqrt_var)
  return part1*part2/part3

def integrand_for_disp_els_M_2_3(n_j,n_0):
  z=n_j/n_0
  part1 = 512*2/n_0*(3*z*z+26*z+28)/(pow(z,7))
  z_1=z-1
  if(n_j < n_0):
    sqrt_var = numpy.sqrt(abs(z_1))
    temp = numpy.arctan(sqrt_var)/sqrt_var - 2*1/3.0*z_1*special.hyp2f1(0.75,1,1.75,z_1*z_1)
    part2 = numpy.exp(-12*temp)
    part3 = 1-numpy.exp(-6*math.pi/sqrt_var)
  else:
    sqrt_var = numpy.sqrt(z_1)
    part2 = numpy.exp(-12/sqrt_var*numpy.arctan(sqrt_var))
    part3 = 1-numpy.exp(-6*math.pi/sqrt_var)
  return part1*part2/part3

def integrand_for_disp_els_M_4_5(n_j,n_0):
  z=n_j/n_0
  part1 = 1024/5.0*2/n_0*(5*z*z+46*z+48)/(pow(z,8))
  z_1=z-1
  if(n_j < n_0):
    sqrt_var = numpy.sqrt(abs(z_1))
    temp = numpy.arctan(sqrt_var)/sqrt_var - 2*1/3.0*z_1*special.hyp2f1(0.75,1,1.75,z_1*z_1)
    part2 = numpy.exp(-12*temp)
    part3 = 1-numpy.exp(-6*math.pi/sqrt_var)
  else:
    sqrt_var = numpy.sqrt(z_1)
    part2 = numpy.exp(-12/sqrt_var*numpy.arctan(sqrt_var))
    part3 = 1-numpy.exp(-6*math.pi/sqrt_var)
  return part1*part2/part3

def real_general_integration(n_j,chi_j,nu,n_0,oscillator_density_function):
  res = (pow(n_j,2)-pow(nu,2))/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*oscillator_density_function(n_j,n_0)

def imag_general_integration(n_j,chi_j,nu,n_0,oscillator_density_function):
  res = -nu*chi_j/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*oscillator_density_function(n_j,n_0)

def k_disp_els(Z=None):
    # Z is integer number of 
  if type(Z) != type(int(20)):
    print("Z MUST BE INTEGER!")
    return
  if Z == None:
    Z = 6
  Z_s_sq = pow(get_Zeff_1s(Z),2)
  e_ion = get_ionization_energy_1s(Z)
  nu_k = e_ion / h
  delta_K = 1 + alpha_sq * Z_s_sq / 4 - e_ion/(Ryd_ener * Z_s_sq) #21a in Hoenl
  n_0 = nu_k/(1-delta_K) #22 in Hoenl

  #this yields the number of dispersion electrons, which is independant of the wavelenght used
  n_disp_el = integrate.quad(integrand_for_disp_els_K,nu_k,2000*nu_k,args=(n_0),limit=200000,epsabs=1E-60)[0]

  return 2*n_disp_el

def l_disp_els(Z=None):
  # Z is integer number of 
  if type(Z) != type(int(20)):
    print("Z MUST BE INTEGER!")
    return
  if Z == None:
    Z = 6
  Z_s_sq1 = pow(get_Zeff_2s(Z),2)
  Z_s_sq2 = pow(get_Zeff_2p_1_2(Z),2)
  Z_s_sq3 = pow(get_Zeff_2p_3_2(Z),2)
  e_ion_1 = get_ionization_energy_2s(Z)
  e_ion_2 = get_ionization_energy_2p1_2(Z)
  e_ion_3 = get_ionization_energy_2p3_2(Z)
  nu_l1 = e_ion_1 / h
  nu_l2 = e_ion_2 / h
  nu_l3 = e_ion_3 / h
  delta_l1 = 1 + alpha_sq * Z_s_sq1 * 0.3125 - 4*e_ion_1/(Ryd_ener * Z_s_sq1) #33 in Hoenl
  delta_l2 = 1 + alpha_sq * Z_s_sq2 * 0.3125 - 4*e_ion_2/(Ryd_ener * Z_s_sq2) #33 in Hoenl
  delta_l3 = 1 + alpha_sq * Z_s_sq3 * 0.0625 - 4*e_ion_3/(Ryd_ener * Z_s_sq3) #33 in Hoenl
  n_0_1 = nu_l1/(1-delta_l1) #22 in Hoenl
  n_0_2 = nu_l2/(1-delta_l2) #22 in Hoenl
  n_0_3 = nu_l3/(1-delta_l3) #22 in Hoenl

  n_disp_2s = integrate.quad(integrand_for_disp_els_L1,nu_l1,20000*nu_l1,args=(n_0_1),limit=200000,epsabs=1E-60)[0]
  n_disp_2p1_2 = integrate.quad(integrand_for_disp_els_L2_3,nu_l2,20000*nu_l2,args=(n_0_2),limit=200000,epsabs=1E-60)[0]
  n_disp_2p3_2 = integrate.quad(integrand_for_disp_els_L2_3,nu_l3,20000*nu_l3,args=(n_0_3),limit=200000,epsabs=1E-60)[0]

  return 2*n_disp_2s+2*n_disp_2p1_2+4*n_disp_2p3_2

def m_disp_els(Z=None):
  # Z is integer number of 
  if type(Z) != type(int(20)):
    print("Z MUST BE INTEGER!")
    return
  if Z == None:
    Z = 52
  Z_s_sq1 = pow(get_Zeff_3s(Z),2)
  Z_s_sq2 = pow(get_Zeff_3p_1_2(Z),2)
  Z_s_sq3 = pow(get_Zeff_3p_3_2(Z),2)
  Z_s_sq4 = pow(get_Zeff_3d_3_2(Z),2)
  Z_s_sq5 = pow(get_Zeff_3d_5_2(Z),2)
  e_ion_1 = get_ionization_energy_3s(Z)
  e_ion_2 = get_ionization_energy_3p_1_2(Z)
  e_ion_3 = get_ionization_energy_3p_3_2(Z)
  e_ion_4 = get_ionization_energy_3d_3_2(Z)
  e_ion_5 = get_ionization_energy_3d_5_2(Z)
  nu_m1 = e_ion_1 / h
  nu_m2 = e_ion_2 / h
  nu_m3 = e_ion_3 / h
  nu_m4 = e_ion_4 / h
  nu_m5 = e_ion_5 / h
  delta_m1 = 1 + alpha_sq * Z_s_sq1 * 0.25 - 9*e_ion_1/(Ryd_ener * Z_s_sq1) #33 in Hoenl
  delta_m2 = 1 + alpha_sq * Z_s_sq2 * 0.25 - 9*e_ion_2/(Ryd_ener * Z_s_sq2) #33 in Hoenl
  delta_m3 = 1 + alpha_sq * Z_s_sq3 / 12 - 9*e_ion_3/(Ryd_ener * Z_s_sq3) #33 in Hoenl
  delta_m4 = 1 + alpha_sq * Z_s_sq4 / 12 - 9*e_ion_4/(Ryd_ener * Z_s_sq4) #33 in Hoenl
  delta_m5 = 1 + alpha_sq * Z_s_sq5 / 36 - 9*e_ion_5/(Ryd_ener * Z_s_sq5) #33 in Hoenl
  n_0_1 = nu_m1/(1-delta_m1) #22 in Hoenl
  n_0_2 = nu_m2/(1-delta_m2) #22 in Hoenl
  n_0_3 = nu_m3/(1-delta_m3) #22 in Hoenl
  n_0_4 = nu_m4/(1-delta_m4) #22 in Hoenl
  n_0_5 = nu_m5/(1-delta_m5) #22 in Hoenl

  n_disp_3s = integrate.quad(integrand_for_disp_els_M1,nu_m1,20000*nu_m1,args=(n_0_1),limit=200000,epsabs=1E-60)[0]
  n_disp_3p1_2 = integrate.quad(integrand_for_disp_els_M_2_3,nu_m2,20000*nu_m2,args=(n_0_2),limit=200000,epsabs=1E-60)[0]
  n_disp_3p3_2 = integrate.quad(integrand_for_disp_els_M_2_3,nu_m3,20000*nu_m3,args=(n_0_3),limit=200000,epsabs=1E-60)[0]
  n_disp_3d3_2 = integrate.quad(integrand_for_disp_els_M_4_5,nu_m4,20000*nu_m4,args=(n_0_4),limit=200000,epsabs=1E-60)[0]
  n_disp_3d5_2 = integrate.quad(integrand_for_disp_els_M_4_5,nu_m5,20000*nu_m5,args=(n_0_5),limit=200000,epsabs=1E-60)[0]

  return 2*n_disp_3s+2*n_disp_3p1_2+4*n_disp_3p3_2+4*n_disp_3d3_2+6*n_disp_3d5_2

def sugiura_k(Z=None, ener=8047.8):
    # Z is integer number of 
  if type(Z) != type(int(20)):
    print("Z MUST BE INTEGER!")
    return
  if Z == None:
    Z = 6
  Z_s_sq = pow(get_Zeff_1s(Z),2)
  e_ion = get_ionization_energy_1s(Z)
  nu_k = e_ion / h
  delta_K = 1 + alpha_sq * Z_s_sq / 4 - e_ion/(Ryd_ener * Z_s_sq) #21a in Hoenl
  n_0 = nu_k/(1-delta_K) #22 in Hoenl
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  epsilon2 = 0.005
  x_j = nu_k*1E-6
  #integral_test = integrate.quad(integrand_damped_abs,nu_k,50*nu,args=(x_j,nu,n_0),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  #integrating directly as one intergal with the singulary position as point does not work.

  #REAL PART
  #When telling the quad alorithm that nu is a "position to be carefull with" the results don't require epsilon, one can integrate from left and right
  inte1 = integrate.quad(real_general_integration,nu_k,nu,args=(x_j,nu,n_0,integrand_for_disp_els_K),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2 = integrate.quad(real_general_integration,nu,40*nu,args=(x_j,nu,n_0,integrand_for_disp_els_K),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  if inte1[1] > 0.5*abs(inte1[0]):
    print("!!! inaccurate REAL1 integral at nu_k/nu = " + "{:8.4f}".format(nu_k/nu) + " " + str(inte1) + str(inte1[1]/abs(inte1[0])))
  if inte2[1] > 0.5*abs(inte2[0]):
    print("!!! inaccurate REAL2 integral at nu_k/nu = " + "{:8.4f}".format(nu_k/nu) + " " + str(inte2) + str(inte2[1]/abs(inte2[0])))
  integral = inte1[0] + inte2[0]
  #IMAG PART
  lower_limit = nu_k
  if (1-epsilon2)*nu > nu_k:
    lower_limit = (1-epsilon2)*nu
  if (1+epsilon2)*nu < nu_k:
    upper_limit = 10*nu_k
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral = integrate.quad(imag_general_integration,lower_limit,upper_limit,args=(x_j,nu,n_0,integrand_for_disp_els_K),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10,full_output=2)
  if imag_integral[1] > 0.5*abs(imag_integral[0]):
    print("!!! inaccurate IMAG integral at nu_k/nu = " + "{:6.4f}".format(nu_k/nu) + " " + str(imag_integral) + " " + str(imag_integral[1]/abs(imag_integral[0])) + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))

  #this yields the number of dispersion electrons, which is independant of the wavelenght used
  #n_disp_el = integrate.quad(integrand_for_disp_els,nu_k,2000*nu_k,args=(n_0),limit=200000,epsabs=1E-60)[0]
  #print(n_disp_el)

  alpha_K_sugiura_damp = complex(integral,imag_integral[0])
  ausdruck = alpha_K_sugiura_damp/(1-4.0*math.pi/3*alpha_K_sugiura_damp*prefactor)
  factor = pow(nu,2)
  return -ausdruck.real*factor, -ausdruck.imag*factor

def l_with_imag(Z=None, ener=8047.8):
  # Z is integer number of 
  if type(Z) != type(int(20)):
    print("Z MUST BE INTEGER!")
    return
  if Z == None:
    Z = 6
  Z_s_sq1 = pow(get_Zeff_2s(Z),2)
  Z_s_sq2 = pow(get_Zeff_2p_1_2(Z),2)
  Z_s_sq3 = pow(get_Zeff_2p_3_2(Z),2)
  e_ion_1 = get_ionization_energy_2s(Z)
  e_ion_2 = get_ionization_energy_2p1_2(Z)
  e_ion_3 = get_ionization_energy_2p3_2(Z)
  nu_l1 = e_ion_1 / h
  nu_l2 = e_ion_2 / h
  nu_l3 = e_ion_3 / h
  delta_l1 = 1 + alpha_sq * Z_s_sq1 * 0.3125 - 4*e_ion_1/(Ryd_ener * Z_s_sq1) #33 in Hoenl
  delta_l2 = 1 + alpha_sq * Z_s_sq2 * 0.3125 - 4*e_ion_2/(Ryd_ener * Z_s_sq2) #33 in Hoenl
  delta_l3 = 1 + alpha_sq * Z_s_sq3 * 0.0625 - 4*e_ion_3/(Ryd_ener * Z_s_sq3) #33 in Hoenl
  n_0_1 = nu_l1/(1-delta_l1) #22 in Hoenl
  n_0_2 = nu_l2/(1-delta_l2) #22 in Hoenl
  n_0_3 = nu_l3/(1-delta_l3) #22 in Hoenl
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  epsilon2 = 0.0025
  x_j1 = nu_l1*1E-6
  x_j2 = nu_l2*1E-6
  x_j3 = nu_l3*1E-6

  #REAL PART
  #When telling hte quad alorithm that nu is a "position to be carefull with" the results don't require epsilon
  inte1_1 = integrate.quad(real_general_integration,nu_l1,nu,args=(x_j1,nu,n_0_1,integrand_for_disp_els_L1),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2_1 = integrate.quad(real_general_integration,nu,40*nu,args=(x_j1,nu,n_0_1,integrand_for_disp_els_L1),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  if inte1_1[1] > 0.5*abs(inte1_1[0]):
    print("!!! inaccurate REAL1 integral at nu_l/nu = " + "{:8.4f}".format(nu_l1/nu) + " " + str(inte1_1) + str(inte1_1[1]/abs(inte1_1[0])))
  if inte2_1[1] > 0.5*abs(inte2_1[0]):
    print("!!! inaccurate REAL2 integral at nu_l/nu = " + "{:8.4f}".format(nu_l1/nu) + " " + str(inte2_1) + str(inte2_1[1]/abs(inte2_1[0])))
  integral_1 = inte1_1[0] + inte2_1[0]

  inte1_2 = integrate.quad(real_general_integration,nu_l2,nu,args=(x_j2,nu,n_0_2,integrand_for_disp_els_L2_3),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2_2 = integrate.quad(real_general_integration,nu,40*nu,args=(x_j2,nu,n_0_2,integrand_for_disp_els_L2_3),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  if inte1_2[1] > 0.5*abs(inte1_2[0]):
    print("!!! inaccurate REAL1 integral at nu_l/nu = " + "{:8.4f}".format(nu_l2/nu) + " " + str(inte1_2) + str(inte1_2[1]/abs(inte1_2[0])))
  if inte2_2[1] > 0.5*abs(inte2_2[0]):
    print("!!! inaccurate REAL2 integral at nu_l/nu = " + "{:8.4f}".format(nu_l2/nu) + " " + str(inte2_2) + str(inte2_2[1]/abs(inte2_2[0])))
  integral_2 = inte1_2[0] + inte2_2[0]

  inte1_3 = integrate.quad(real_general_integration,nu_l3,nu,args=(x_j3,nu,n_0_3,integrand_for_disp_els_L2_3),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2_3 = integrate.quad(real_general_integration,nu,40*nu,args=(x_j3,nu,n_0_3,integrand_for_disp_els_L2_3),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  if inte1_3[1] > 0.5*abs(inte1_3[0]):
    print("!!! inaccurate REAL1 integral at nu_l/nu = " + "{:8.4f}".format(nu_l3/nu) + " " + str(inte1_3) + str(inte1_3[1]/abs(inte1_3[0])))
  if inte2_3[1] > 0.5*abs(inte2_3[0]):
    print("!!! inaccurate REAL2 integral at nu_l/nu = " + "{:8.4f}".format(nu_l3/nu) + " " + str(inte2_3) + str(inte2_3[1]/abs(inte2_3[0])))
  integral_3 = inte1_3[0] + inte2_3[0]

  #IMAG PART
  lower_limit = nu_l1
  if (1-epsilon2)*nu > nu_l1:
    lower_limit = (1-epsilon2)*nu
  if (1+epsilon2)*nu < nu_l1:
    upper_limit = 10*nu_l1
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral_1 = integrate.quad(imag_general_integration,lower_limit,upper_limit,args=(x_j1,nu,n_0_1,integrand_for_disp_els_L1),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10,full_output=2)
  if imag_integral_1[1] > 0.5*abs(imag_integral_1[0]):
    print("!!! inaccurate IMAG integral at nu_l/nu = " + "{:6.4f}".format(nu_l1/nu) + " " + str(imag_integral_1) + " " + str(imag_integral_1[1]/abs(imag_integral_1[0])) 
    + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))

  lower_limit = nu_l2
  if (1-epsilon2)*nu > nu_l2:
    lower_limit = (1-epsilon2)*nu
  if (1+epsilon2)*nu < nu_l2:
    upper_limit = 10*nu_l2
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral_2 = integrate.quad(imag_general_integration,lower_limit,upper_limit,args=(x_j2,nu,n_0_2,integrand_for_disp_els_L2_3),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10,full_output=2)
  if imag_integral_2[1] > 0.5*abs(imag_integral_2[0]):
    print("!!! inaccurate IMAG integral at nu_l/nu = " + "{:6.4f}".format(nu_l2/nu) + " " + str(imag_integral_2) + " " + str(imag_integral_2[1]/abs(imag_integral_2[0])) 
    + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))
  
  lower_limit = nu_l3
  if (1-epsilon2)*nu > nu_l3:
    lower_limit = (1-epsilon2)*nu
  if (1+epsilon2)*nu < nu_l3:
    upper_limit = 10*nu_l3
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral_3 = integrate.quad(imag_general_integration,lower_limit,upper_limit,args=(x_j3,nu,n_0_3,integrand_for_disp_els_L2_3),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10,full_output=2)
  if imag_integral_3[1] > 0.5*abs(imag_integral_3[0]):
    print("!!! inaccurate IMAG integral at nu_l/nu = " + "{:6.4f}".format(nu_l3/nu) + " " + str(imag_integral_3) + " " + str(imag_integral_3[1]/abs(imag_integral_3[0])) 
    + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))

  #n_disp_2s = integrate.quad(integrand_for_disp_els_L1,nu_l1,2000*nu_l1,args=(n_0_1),limit=200000,epsabs=1E-60)[0]
  #n_disp_2p1_2 = integrate.quad(integrand_for_disp_els_L2_3,nu_l2,2000*nu_l2,args=(n_0_2),limit=200000,epsabs=1E-60)[0]
  #n_disp_2p3_2 = integrate.quad(integrand_for_disp_els_L2_3,nu_l3,2000*nu_l3,args=(n_0_3),limit=200000,epsabs=1E-60)[0]

  alpha_L_sugiura_damp = complex(integral_1+integral_2+integral_3*2,imag_integral_1[0]+imag_integral_2[0]+2*imag_integral_3[0])
  ausdruck = (alpha_L_sugiura_damp)/(1-4.0*math.pi/3*alpha_L_sugiura_damp*prefactor)
  factor = pow(nu,2)
  return -ausdruck.real*factor, -ausdruck.imag*factor

def m_with_imag(Z=None, ener=8047.8):
  # Z is integer number of 
  if type(Z) != type(int(20)):
    print("Z MUST BE INTEGER!")
    return
  if Z == None:
    Z = 6
  Z_s_sq1 = pow(get_Zeff_3s(Z),2)
  Z_s_sq2 = pow(get_Zeff_3p_1_2(Z),2)
  Z_s_sq3 = pow(get_Zeff_3p_3_2(Z),2)
  Z_s_sq4 = pow(get_Zeff_3d_3_2(Z),2)
  Z_s_sq5 = pow(get_Zeff_3d_5_2(Z),2)
  e_ion_1 = get_ionization_energy_3s(Z)
  e_ion_2 = get_ionization_energy_3p_1_2(Z)
  e_ion_3 = get_ionization_energy_3p_3_2(Z)
  e_ion_4 = get_ionization_energy_3d_3_2(Z)
  e_ion_5 = get_ionization_energy_3d_5_2(Z)
  nu_m1 = e_ion_1 / h
  nu_m2 = e_ion_2 / h
  nu_m3 = e_ion_3 / h
  nu_m4 = e_ion_4 / h
  nu_m5 = e_ion_5 / h
  delta_m1 = 1 + alpha_sq * Z_s_sq1 * 0.25 - 9*e_ion_1/(Ryd_ener * Z_s_sq1) #33 in Hoenl
  delta_m2 = 1 + alpha_sq * Z_s_sq2 * 0.25 - 9*e_ion_2/(Ryd_ener * Z_s_sq2) #33 in Hoenl
  delta_m3 = 1 + alpha_sq * Z_s_sq3 / 12 - 9*e_ion_3/(Ryd_ener * Z_s_sq3) #33 in Hoenl
  delta_m4 = 1 + alpha_sq * Z_s_sq4 / 12 - 9*e_ion_2/(Ryd_ener * Z_s_sq4) #33 in Hoenl
  delta_m5 = 1 + alpha_sq * Z_s_sq5 / 36 - 9*e_ion_3/(Ryd_ener * Z_s_sq5) #33 in Hoenl
  n_0_1 = nu_m1/(1-delta_m1) #22 in Hoenl
  n_0_2 = nu_m2/(1-delta_m2) #22 in Hoenl
  n_0_3 = nu_m3/(1-delta_m3) #22 in Hoenl
  n_0_4 = nu_m4/(1-delta_m4) #22 in Hoenl
  n_0_5 = nu_m5/(1-delta_m5) #22 in Hoenl
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  epsilon2 = 0.0025
  x_j1 = nu_m1*1E-6
  x_j2 = nu_m2*1E-6
  x_j3 = nu_m3*1E-6
  x_j4 = nu_m4*1E-6
  x_j5 = nu_m5*1E-6

  #REAL PART
  #When telling hte quad alorithm that nu is a "position to be carefull with" the results don't require epsilon
  inte1_1 = integrate.quad(real_general_integration,nu_m1,nu,args=(x_j1,nu,n_0_1,integrand_for_disp_els_M1),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2_1 = integrate.quad(real_general_integration,nu,40*nu,args=(x_j1,nu,n_0_1,integrand_for_disp_els_M1),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  if inte1_1[1] > 0.5*abs(inte1_1[0]):
    print("!!! inaccurate REAL1 integral at nu_m1/nu = " + "{:8.4f}".format(nu_m1/nu) + " " + str(inte1_1) + str(inte1_1[1]/abs(inte1_1[0])))
  if inte2_1[1] > 0.5*abs(inte2_1[0]):
    print("!!! inaccurate REAL2 integral at nu_m1/nu = " + "{:8.4f}".format(nu_m1/nu) + " " + str(inte2_1) + str(inte2_1[1]/abs(inte2_1[0])))
  integral_1 = inte1_1[0] + inte2_1[0]

  inte1_2 = integrate.quad(real_general_integration,nu_m2,nu,args=(x_j2,nu,n_0_2,integrand_for_disp_els_M_2_3),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2_2 = integrate.quad(real_general_integration,nu,40*nu,args=(x_j2,nu,n_0_2,integrand_for_disp_els_M_2_3),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  if inte1_2[1] > 0.5*abs(inte1_2[0]):
    print("!!! inaccurate REAL1 integral at nu_m2/nu = " + "{:8.4f}".format(nu_m2/nu) + " " + str(inte1_2) + str(inte1_2[1]/abs(inte1_2[0])))
  if inte2_2[1] > 0.5*abs(inte2_2[0]):
    print("!!! inaccurate REAL2 integral at nu_m2/nu = " + "{:8.4f}".format(nu_m2/nu) + " " + str(inte2_2) + str(inte2_2[1]/abs(inte2_2[0])))
  integral_2 = inte1_2[0] + inte2_2[0]

  inte1_3 = integrate.quad(real_general_integration,nu_m3,nu,args=(x_j3,nu,n_0_3,integrand_for_disp_els_M_2_3),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2_3 = integrate.quad(real_general_integration,nu,40*nu,args=(x_j3,nu,n_0_3,integrand_for_disp_els_M_2_3),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  if inte1_3[1] > 0.5*abs(inte1_3[0]):
    print("!!! inaccurate REAL1 integral at nu_m3/nu = " + "{:8.4f}".format(nu_m3/nu) + " " + str(inte1_3) + str(inte1_3[1]/abs(inte1_3[0])))
  if inte2_3[1] > 0.5*abs(inte2_3[0]):
    print("!!! inaccurate REAL2 integral at nu_m3/nu = " + "{:8.4f}".format(nu_m3/nu) + " " + str(inte2_3) + str(inte2_3[1]/abs(inte2_3[0])))
  integral_3 = inte1_3[0] + inte2_3[0]

  inte1_4 = integrate.quad(real_general_integration,nu_m3,nu,args=(x_j4,nu,n_0_4,integrand_for_disp_els_M_4_5),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2_4 = integrate.quad(real_general_integration,nu,40*nu,args=(x_j4,nu,n_0_4,integrand_for_disp_els_M_4_5),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  if inte1_4[1] > 0.5*abs(inte1_4[0]):
    print("!!! inaccurate REAL1 integral at nu_m3/nu = " + "{:8.4f}".format(nu_m3/nu) + " " + str(inte1_4) + str(inte1_4[1]/abs(inte1_4[0])))
  if inte2_4[1] > 0.5*abs(inte2_4[0]):
    print("!!! inaccurate REAL2 integral at nu_m3/nu = " + "{:8.4f}".format(nu_m3/nu) + " " + str(inte2_4) + str(inte2_4[1]/abs(inte2_4[0])))
  integral_4 = inte1_4[0] + inte2_4[0]

  inte1_5 = integrate.quad(real_general_integration,nu_m3,nu,args=(x_j5,nu,n_0_5,integrand_for_disp_els_M_4_5),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2_5 = integrate.quad(real_general_integration,nu,40*nu,args=(x_j5,nu,n_0_5,integrand_for_disp_els_M_4_5),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  if inte1_5[1] > 0.5*abs(inte1_5[0]):
    print("!!! inaccurate REAL1 integral at nu_m3/nu = " + "{:8.4f}".format(nu_m3/nu) + " " + str(inte1_5) + str(inte1_5[1]/abs(inte1_5[0])))
  if inte2_5[1] > 0.5*abs(inte2_5[0]):
    print("!!! inaccurate REAL2 integral at nu_m3/nu = " + "{:8.4f}".format(nu_m3/nu) + " " + str(inte2_5) + str(inte2_5[1]/abs(inte2_5[0])))
  integral_5 = inte1_5[0] + inte2_5[0]

  #IMAG PART
  lower_limit = nu_m1
  if (1-epsilon2)*nu > nu_m1:
    lower_limit = (1-epsilon2)*nu
  if (1+epsilon2)*nu < nu_m1:
    upper_limit = 10*nu_m1
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral_1 = integrate.quad(imag_general_integration,lower_limit,upper_limit,args=(x_j1,nu,n_0_1,integrand_for_disp_els_M1),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10,full_output=2)
  if imag_integral_1[1] > 0.5*abs(imag_integral_1[0]):
    print("!!! inaccurate IMAG integral at nu_m1/nu = " + "{:6.4f}".format(nu_m1/nu) + " " + str(imag_integral_1) + " " + str(imag_integral_1[1]/abs(imag_integral_1[0])) 
    + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))

  lower_limit = nu_m2
  if (1-epsilon2)*nu > nu_m2:
    lower_limit = (1-epsilon2)*nu
  if (1+epsilon2)*nu < nu_m2:
    upper_limit = 10*nu_m2
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral_2 = integrate.quad(imag_general_integration,lower_limit,upper_limit,args=(x_j2,nu,n_0_2,integrand_for_disp_els_M_2_3),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10,full_output=2)
  if imag_integral_2[1] > 0.5*abs(imag_integral_2[0]):
    print("!!! inaccurate IMAG integral at nu_m2/nu = " + "{:6.4f}".format(nu_m2/nu) + " " + str(imag_integral_2) + " " + str(imag_integral_2[1]/abs(imag_integral_2[0])) 
    + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))
  
  lower_limit = nu_m3
  if (1-epsilon2)*nu > nu_m3:
    lower_limit = (1-epsilon2)*nu
  if (1+epsilon2)*nu < nu_m3:
    upper_limit = 10*nu_m3
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral_3 = integrate.quad(imag_general_integration,lower_limit,upper_limit,args=(x_j3,nu,n_0_3,integrand_for_disp_els_M_2_3),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10,full_output=2)
  if imag_integral_3[1] > 0.5*abs(imag_integral_3[0]):
    print("!!! inaccurate IMAG integral at nu_m3/nu = " + "{:6.4f}".format(nu_m3/nu) + " " + str(imag_integral_3) + " " + str(imag_integral_3[1]/abs(imag_integral_3[0])) 
    + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))

  lower_limit = nu_m4
  if (1-epsilon2)*nu > nu_m4:
    lower_limit = (1-epsilon2)*nu
  if (1+epsilon2)*nu < nu_m4:
    upper_limit = 10*nu_m4
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral_4 = integrate.quad(imag_general_integration,lower_limit,upper_limit,args=(x_j4,nu,n_0_4,integrand_for_disp_els_M_4_5),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10,full_output=2)
  if imag_integral_4[1] > 0.5*abs(imag_integral_4[0]):
    print("!!! inaccurate IMAG integral at nu_m3/nu = " + "{:6.4f}".format(nu_m4/nu) + " " + str(imag_integral_4) + " " + str(imag_integral_4[1]/abs(imag_integral_4[0])) 
    + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))

  lower_limit = nu_m5
  if (1-epsilon2)*nu > nu_m5:
    lower_limit = (1-epsilon2)*nu
  if (1+epsilon2)*nu < nu_m5:
    upper_limit = 10*nu_m5
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral_5 = integrate.quad(imag_general_integration,lower_limit,upper_limit,args=(x_j5,nu,n_0_5,integrand_for_disp_els_M_4_5),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10,full_output=2)
  if imag_integral_5[1] > 0.5*abs(imag_integral_5[0]):
    print("!!! inaccurate IMAG integral at nu_m3/nu = " + "{:6.4f}".format(nu_m5/nu) + " " + str(imag_integral_5) + " " + str(imag_integral_5[1]/abs(imag_integral_5[0])) 
    + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))

  alpha_M_sugiura_damp = complex(integral_1+integral_2+integral_3*2+integral_4*3+integral_5*2,imag_integral_1[0]+imag_integral_2[0]+2*imag_integral_3[0]+3*imag_integral_4[0]+2*imag_integral_5[0])
  ausdruck = (alpha_M_sugiura_damp)/(1-4.0*math.pi/3*alpha_M_sugiura_damp*prefactor)
  factor = pow(nu,2)
  return -ausdruck.real*factor, -ausdruck.imag*factor

def test_florian():
  a=52
  #
  #x = numpy.linspace(3,113,110)
  #y1 = []
  #y2 = []
  #y3 = []
  #y4 = []
  #for i in x:
  #  y1.append(get_Zeff_1s(int(i)))
  #  y2.append(get_Zeff_2s(int(i)))
  #  y3.append(get_Zeff_2p_1_2(int(i)))
  #  y4.append(get_Zeff_2p_3_2(int(i)))
  #
  #plt.plot(x,y1,label='1s')
  #plt.plot(x,y2,label='2s')
  #plt.plot(x,y3,label='2p')
  #plt.plot(x,y3,label='2p2')
  #plt.legend()
  #plt.show()
  #
  Z_s_sq = pow(get_Zeff_1s(a),2)
  ZK = get_Zeff_1s(a)
  e_ion = get_ionization_energy_1s(a)
  nu_k = e_ion / h
  delta_K = 1 + alpha_sq * Z_s_sq / 4 - e_ion/(Ryd_ener * Z_s_sq) #21a in Hoenl
  n_0 = nu_k/(1-delta_K) #22 in Hoenl
  #this yields the number of dispersion electrons, which is independant of the wavelenght used
  n_disp_el = integrate.quad(integrand_for_disp_els_K,nu_k,20000*nu_k,args=(n_0),limit=200000,epsabs=1E-60)
  print("K_electrons: " + str(2*n_disp_el[0]) +  " delta_K: " + str(delta_K))
  print("Missed part right: " + str(integrate.quad(integrand_for_disp_els_K,20000*nu_k,50000*nu_k,args=(n_0),limit=2000000,epsabs=1E-30,epsrel=1E-10)) )
  #
  #
  ##en = get_ionization_energy_2s(a)
  ZL1 = get_Zeff_2s(a)
  Z_s_sq1 = pow(ZL1,2)
  ZL2 = get_Zeff_2p_1_2(a)
  Z_s_sq2 = pow(ZL2,2) #according to the tables there is no difference between the effective core charge for 2p electrons
  ZL3 = get_Zeff_2p_3_2(a)
  Z_s_sq3 = pow(ZL3,2)
  e_ion_1 = get_ionization_energy_2s(a)
  e_ion_2 = get_ionization_energy_2p1_2(a)
  e_ion_3 = get_ionization_energy_2p3_2(a)
  nu_l1 = e_ion_1 / h
  nu_l2 = e_ion_2 / h
  nu_l3 = e_ion_3 / h
  delta_l1 = 1 + alpha_sq * Z_s_sq1 * 0.3125 - 4*e_ion_1/(Ryd_ener * Z_s_sq1) #33 in Hoenl
  delta_l2 = 1 + alpha_sq * Z_s_sq2 * 0.3125 - 4*e_ion_2/(Ryd_ener * Z_s_sq2) #33 in Hoenl
  delta_l3 = 1 + alpha_sq * Z_s_sq3 * 0.0625 - 4*e_ion_3/(Ryd_ener * Z_s_sq3) #33 in Hoenl
  n_0_1 = nu_l1/(1-delta_l1) #22 in Hoenl
  n_0_2 = nu_l2/(1-delta_l2) #22 in Hoenl
  n_0_3 = nu_l3/(1-delta_l3) #22 in Hoenl
  
  #plt.axvline(x=nu_l1, color="blue", linestyle="--")
  #plt.axvline(x=nu_l2, color="orange", linestyle="--")
  #plt.axvline(x=nu_l3, color="green", linestyle="--")
  #if delta_l1 <0:
  #  lower_limit=(1+delta_l1)*nu_l1/n_0_1
  #else:
  lower_limit=nu_l1
  upper_limit=20000*nu_l1
  #plt.axvline(x=upper_limit, color="gray", linestyle="--")
  #x = numpy.linspace(lower_limit,upper_limit,5000)
  #y = []
  inte1 = integrate.quad(integrand_for_disp_els_L1,lower_limit,upper_limit,args=(n_0_1),limit=200000,epsabs=1E-60,epsrel=1E-10)
  print("LI n_disp: " +str(inte1) + " delta_L: " + str(delta_l1))
  #inte1_2 = integrate.quad(integrand_for_disp_els_L1,nu_l1,upper_limit,args=(n_0_1),limit=200000,epsabs=1E-60,epsrel=1E-10)
  #print("nu_l1 Imag: " +str(inte1))
  print("Missed part right: " + str(integrate.quad(integrand_for_disp_els_L1,upper_limit,50000*nu_l1,args=(n_0_1),limit=2000000,epsabs=1E-30,epsrel=1E-10)) )
  #for i in x:
  #  y.append(integrand_for_disp_els_L1(i,n_0_1))
  #plt.plot(x,y,'+:',label="LI")
  
  lower_limit=nu_l2
  #x = numpy.linspace(lower_limit,upper_limit,5000)
  #y = []
  inte2 = integrate.quad(integrand_for_disp_els_L2_3,lower_limit,upper_limit,args=(n_0_2),limit=200000,epsabs=1E-60,epsrel=1E-10)
  print("LII n_disp: " +str(inte2) + " delta_L: " + str(delta_l2))
  #for i in x:
  #  y.append(integrand_for_disp_els_L2_3(i,n_0_2))
  #plt.plot(x,y,'+:',label="LII")
  
  lower_limit=nu_l3
  #x = numpy.linspace(lower_limit,upper_limit,5000)
  #y = []
  inte3 = integrate.quad(integrand_for_disp_els_L2_3,lower_limit,upper_limit,args=(n_0_3),limit=200000,epsabs=1E-60,epsrel=1E-10)
  print("LIII n_disp: " +str(inte3) + " delta_L: " + str(delta_l3))
  #for i in x:
  #  y.append(integrand_for_disp_els_L2_3(i,n_0_3))
  #plt.plot(x,y,'+:',label="LIII")
  print("Number: " + str(2*inte1[0]+2*inte2[0]+4*inte3[0]))
  #plt.axhline(y=0, color="gray", linestyle="--")
  #plt.legend(loc='upper right')
  #plt.show()
  #exit()
  print("done")

#test_florian()

if __name__ == '__main__':
  import matplotlib.pyplot as plt
  fig, axs = plt.subplots(4,2,sharex=True)
  stop = False
  Zs = [42]
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
      nr_steps = 30
      for step in range(int(minimal),int(maximal),int(nr_steps)):
        steps.append(step)
    elif axis == "fixed":
      if i == 52:
        steps = [78,156,233,272,311,350,369,381,388.9,389.1,397,408,428,467,545,558,622,700,778,1167,1200,1556,1700,1945,2100,2400,2450,2480,2504.9,2505.1,2520,2530,2550,2600,2650,2681.9,2682.1,2700,2720,2723,2750,2770,2800,2849.9,2850.1,2900,3000,3500,4000,100000]
      elif i == 74:
        steps = [80,120,150,160,170,175,177.9,178.1,180,185,200,250,400,600,800,900,980,1010,1021.9,1022.1,1030,1040,1050,1065,1071.9,1072.1,1080,1120,1160,1205,1212.9,1213.1,1220,1300,1500,2000]
    n_disp_K = k_disp_els(i)
    n_disp_L = l_disp_els(i)
    n_disp_M = m_disp_els(i)
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
      
    axs[0,0].plot(x,y_s_real,'+:',label="%s K_edge"%elements[i])
    axs[0,1].plot(x,y_s_imag,'+:',label="%s K_edge"%elements[i])

    axs[1,0].plot(x,y_l_real,'+:',label="%s L_edges"%elements[i])
    axs[1,1].plot(x,y_l_imag,'+:',label="%s L_edges"%elements[i])

    axs[2,0].plot(x,y_m_real,'+:',label="%s M_edges"%elements[i])
    axs[2,1].plot(x,y_m_imag,'+:',label="%s M_edges"%elements[i])

    axs[3,0].plot(x,y_d_real,'+:',label="%s"%elements[i])
    axs[3,1].plot(x,y_d_imag,'+:',label="%s"%elements[i])
    axs[3,0].axhline(y=-(n_disp_K+n_disp_L), color="gray", linestyle="--")

  axs[0, 0].set_title('K Real')
  axs[1, 0].set_title('L Real')
  axs[2, 0].set_title('M Real')
  axs[3, 0].set_title('Complete Real')
  axs[0, 1].set_title('K Imag')
  axs[1, 1].set_title('L Imag')
  axs[2, 1].set_title('M Imag')
  axs[3, 1].set_title('Complete Imag')
  axs[0, 1].legend(loc='upper right')
  axs[1, 1].legend(loc='upper right')
  axs[2, 1].legend(loc='upper right')
  axs[3, 1].legend(loc='upper right')
  for j in range(len(edges)):
      axs[0,0].axvline(x=edges[j], color="gray", linestyle=":")
      axs[0,1].axvline(x=edges[j], color="gray", linestyle=":")
      axs[3,0].axvline(x=edges[j], color="gray", linestyle=":")
      axs[3,1].axvline(x=edges[j], color="gray", linestyle=":")
  for j in range(len(L_edges)):
      axs[1,0].axvline(x=L_edges[j], color="gray", linestyle=":")
      axs[1,1].axvline(x=L_edges[j], color="gray", linestyle=":")
      axs[3,0].axvline(x=L_edges[j], color="gray", linestyle=":")
      axs[3,1].axvline(x=L_edges[j], color="gray", linestyle=":")
  for j in range(len(M_edges)):
      axs[2,0].axvline(x=M_edges[j], color="gray", linestyle=":")
      axs[2,1].axvline(x=M_edges[j], color="gray", linestyle=":")
      axs[3,0].axvline(x=M_edges[j], color="gray", linestyle=":")
      axs[3,1].axvline(x=M_edges[j], color="gray", linestyle=":")
  if axis == "Angstrom" or axis == "fixed":
    axs[3, 0].set(xlabel=r'$\lambda\ /\AA$')
    axs[3, 1].set(xlabel=r'$\lambda\ /\AA$')
  elif axis == "keV":
    axs[3, 0].set(xlabel=r'$E /keV$')
    axs[3, 1].set(xlabel=r'$E /keV$')
  elif axis == "relative":
    axs[3, 0].set(xlabel=r'$\frac{\lambda}{\lambda_K}$')
    axs[3, 1].set(xlabel=r'$\frac{\lambda}{\lambda_K}$')
  axs[0, 0].set(ylabel=r'$\frac{1-n_0}{\lambda^2} * \frac{2 \pi m c^2}{N e^2}$')
  axs[1, 0].set(ylabel=r'$\frac{1-n_0}{\lambda^2} * \frac{2 \pi m c^2}{N e^2}$')
  axs[2, 0].set(ylabel=r'$\frac{1-n_0}{\lambda^2} * \frac{2 \pi m c^2}{N e^2}$')
  axs[3, 0].set(ylabel=r'$\Delta f\ /e$')
  import matplotlib
  for ax in fig.get_axes():
    if axis == "Angstrom" or axis == "fixed":
      #ax.set_xscale('log')
      ax.set_xlim(3.0,0.05)
      ax.set_xticks([0.1,0.2,0.3,0.5,0.7,1.0,1.5,2.0,3.0])
      ax.get_xaxis().get_major_formatter().labelOnlyBase = False
      ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    elif axis == "relative":
      ax.set_xlim(0.0,2.0) #for relative spectra
    elif axis == "keV":
      ax.set_xlim(0.1,200) #for energy spectra in keV
  plt.show()

  print("Done")

# OLD STUFF BELOW
def integrand_a_la_hoenl(n_j,nu,n_0):
  var = (n_j/n_0)
  part1 = 128/(3.0*n_0)
  part2 = numpy.exp(-4)/3.0*(4*var-1)
  part3 = 1/pow(var,4)

  res = 1/(pow(n_j,2)-pow(nu,2))
  return res*part1*part2*part3

def integrand_a_la_hoenl_damp(n_j,chi_j,nu,n_0):
  var = (n_j/n_0)
  part1 = 128/(3.0*n_0)
  part2 = numpy.exp(-4)/3.0*(4*var-1)
  part3 = 1/pow(var,4)
  res = (pow(n_j,2)-pow(nu,2))/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*part1*part2*part3

def integrand(n_j,nu,n_0):
  var = (n_j/n_0)-1
  sqrt_var = numpy.sqrt(var)
  part1 = 128*pow(n_0,3)/(3*pow(n_j,4))
  part2 = numpy.exp(-4/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-2*math.pi/sqrt_var)
  res = 1/(pow(n_j,2)-pow(nu,2))
  return res*part1*part2/part3

def integrand_abs(n_j,nu,n_0):
  var = abs((n_j/n_0)-1)
  sqrt_var = numpy.sqrt(var)
  part1 = 128*pow(n_0,3)/(3*pow(n_j,4))
  part2 = numpy.exp(-4/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-2*math.pi/sqrt_var)
  res = 1/(pow(n_j,2)-pow(nu,2))
  return res*part1*part2/part3

def integrand_damped(n_j,chi_j,nu,n_0):
  var = (n_j/n_0)-1
  sqrt_var = numpy.sqrt(var)
  part1 = 128*pow(n_0,3)/(3*pow(n_j,4))
  part2 = numpy.exp(-4/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-2*math.pi/sqrt_var)
  res = (pow(n_j,2)-pow(nu,2))/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*part1*part2/part3

def imag_integrand_damped(n_j,chi_j,nu,n_0):
  var = (n_j/n_0)-1
  sqrt_var = numpy.sqrt(var)
  part1 = 128*pow(n_0,3)/(3*pow(n_j,4))
  part2 = numpy.exp(-4/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-2*math.pi/sqrt_var)
  res = -nu*chi_j/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*part1*part2/part3

def imag_integrand_damped_a_la_hoenl(n_j,chi_j,nu,n_0):
  var = (n_j/n_0)
  part1 = 128/(3.0*n_0)
  part2 = numpy.exp(-4)/3.0*(4*var-1)
  part3 = 1/pow(var,4)

  res = -nu*chi_j/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*part1*part2*part3

def oszillator_sugiura_K(z):
  if(z<1):
    var = numpy.sqrt(abs(z-1))
    zaehler = numpy.exp(-4/var*numpy.arctan(var)) - 2*1/3.0*(z-1)*special.hyp2f1(0.75,1,1.75,(z-1)*(z-1))
  else:
    var = numpy.sqrt(z-1)
    zaehler = numpy.exp(-4/var*numpy.arctan(var))
  nenner = 1-numpy.exp(-2*numpy.pi/var)
  return zaehler/nenner

def hoenl_26b(nu,n_0,n_k):
  if nu >= n_k: 
    oscillator = oszillator_sugiura_K(nu/n_0)
    #factor_26b = numpy.pi*2*pow(el_charge,2)/(4*pow(numpy.pi,2)*el_mass*2*nu)
    #factor_alpha = 2*numpy.pi
    #factor_plot = 2*math.pi*el_mass*pow(nu,2)/pow(el_charge,2)
    #all_factor = factor_26b*factor_alpha*factor_plot/n_0
    factor = math.pi*nu/n_0
    return oscillator*factor
  else:
    return 0.0

def hoenl_27a(nu_K,nu,delta_K):
  x = nu_K/nu
  term1 = -32*pow(el_charge,2)*numpy.exp(-4)/(9*pow(math.pi,2)*el_mass*pow(nu,2))
  term2 = 4/pow(1-delta_K,2)*(1+pow(x,2)*math.log(abs(1-1/pow(x,2))))
  term3 = -1/pow(1-delta_K,3)*(2.0/3.0 + 2*pow(x,2) + pow(x,3)*math.log(abs((1-x)/(1+x))))
  result = term1*(term2+term3)
  return result

def hoenl_27b(nu_K,nu,delta_K):
  x = nu_K/nu
  if x > 1:
    return 0.0
  term1 = -32*pow(el_charge,2)*numpy.exp(-4)/(9*math.pi*el_mass*pow(nu,2))
  term2 = 4/pow(1-delta_K,2)*pow(x,2)
  term3 = -pow(x,3)/pow(1-delta_K,3)
  return term1*(term2+term3)

def hoenl_like_k_print(steps=20, Z=None, ener=8047.8):
  #Z = int(input("Please enter element number: "))
  if Z == None:
    Z = 6
  #s = 5/16
  Z_s_sq = pow(get_Zeff_1s(Z),2)
  e_ion = get_ionization_energy_1s(Z)
  nu_k = e_ion / h
  delta_K = 1 + pow(alpha/2,2) * Z_s_sq - e_ion/(Ryd_ener * Z_s_sq)
  nu_0 = nu_k/(1-delta_K)
  n_0 = nu_0
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  
  print("nu_k:    "+"{:.4e}".format(nu_k)+" e_K:  "+"{:.4e}".format(nu_k*h))
  print("nu_0:    "+"{:.4e}".format(nu_0)+" e_0:  "+"{:.4e}".format(nu_0*h))
  print("n_0:     "+"{:.4e}".format(n_0)+" e_0:  "+"{:.4e}".format(n_0*h))
  print("nu:      "+"{:.4e}".format(nu)+" e_nu: "+"{:.4e}".format(nu*h))
  print("delta_k: "+"{:.4e}".format(delta_K))
  #energies = []
  #alpha_K = []
  epsilon = 1E-12
  #######print("UNDAMPED:")
  ########for i in range(steps):
  #######integral1 = 0.0
  #######integral2 = 0.0
  #######integral = 0.0
  #######if(nu >= n_0):
  #######  integral1 = integrate.quad(integrand_abs,nu_k,(1-epsilon)*nu,limit=20000,epsabs=1E-50)
  #######  integral2 = integrate.quad(integrand,(1+epsilon)*nu,1000*nu,limit=20000,epsabs=1E-50)
  #######  integral = integral1[0] + integral2[0]
  #######else:
  #######  integral = integrate.quad(integrand_abs,nu_k,1000*nu,limit=20000,epsabs=1E-50)[0]
  #######integral1_a_la_hoenl = 0.0
  #######integral2_a_la_hoenl = 0.0
  #######integral_a_la_hoenl = 0.0
  #######if(nu >= n_0):
  #######  integral1_a_la_hoenl = integrate.quad(integrand_a_la_hoenl,nu_k,(1-epsilon)*nu,limit=2000,epsabs=1E-50)
  #######  integral2_a_la_hoenl = integrate.quad(integrand_a_la_hoenl,(1+epsilon)*nu,1000*nu,limit=2000,epsabs=1E-50)
  #######  integral_a_la_hoenl = integral1_a_la_hoenl[0] + integral2_a_la_hoenl[0]
  #######else:
  #######  integral_a_la_hoenl = integrate.quad(integrand_a_la_hoenl,nu_k,1000*nu,limit=2000,epsabs=1E-50)[0]
  #######
  #######print("Epsilon: " + "{:.2e}".format(epsilon) 
  #######+ " nu_0: " + "{:.4e}".format(n_0) +" nu: " + "{:.4e}".format(nu) 
  #######+ " Integral: " + "{:.8e}".format(integral)
  #######+ " PrefInt: " + "{:.8e}".format(prefactor*integral)
  #######+ " IntHoenl: " + "{:.8e}".format(integral_a_la_hoenl)
  #######+ " PrefIntH: " + "{:.8e}".format(prefactor*integral_a_la_hoenl)
  #######)
  #x = numpy.arange(1.01*nu_0,(1-epsilon)*nu,1E16)
  #y = integrand(x)
  #plt.plot(x,y)
  #plt.show()
  epsilon2 = 0.001
  x_j = nu*1E-8
  #x = numpy.arange((1-0.0001)*nu,(1+0.0001)*nu,1E12)
  #y0 = imag_integrand_damped(x,1E12)
  #y = imag_integrand_damped(x,1E13)
  #y2 = imag_integrand_damped(x,1E14)
  #y3 = imag_integrand_damped(x,1E15)
  #plt.plot(x,y0,'r')
  #plt.plot(x,y, 'g')
  #plt.plot(x,y2,'b')
  #plt.plot(x,y3,'y')
  #plt.show()
  print("DAMPED: ")
  #for i in range(steps):
  #nu = 6.4077E18
  #x_j = 0.1 * x_j
  print("x_j:     " + "{:.2e}".format(x_j))
  integral1 = integrate.quad(integrand_damped_abs,nu_k,(1-epsilon)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
  integral2 = integrate.quad(integrand_damped_abs,(1+epsilon)*nu,1000*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
  integral = integral1[0] + integral2[0]
  print(" Int:          " + "{: .8e}".format(integral) +                 " PrefInt:          " + "{: .8e}".format(integral*prefactor))
  integral_a_la_hoenl_damp = integrate.quad(integrand_a_la_hoenl_damp,nu_k,(1-epsilon)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)[0] \
    + integrate.quad(integrand_a_la_hoenl_damp,(1+epsilon)*nu,1000*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)[0]
  print(" IntHoenl:     " + "{: .8e}".format(integral_a_la_hoenl_damp) + " PrefIntHoenl:     " + "{: .8e}".format(integral_a_la_hoenl_damp*prefactor))
  lower_limit = n_0
  if (1-epsilon2)*nu > n_0:
    lower_limit = (1-epsilon2)*nu
  imag_integral1 = integrate.quad(imag_integrand_damped_abs,lower_limit,(1+epsilon2)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
  imag_integral = imag_integral1[0]
  print(" ImagInt:      " + "{: .8e}".format(imag_integral) +            " ImagPrefInt:      " + "{: .8e}".format(imag_integral*prefactor))
  imag_integral_a_la_hoenl = integrate.quad(imag_integrand_damped_a_la_hoenl,lower_limit,(1+epsilon2)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)[0]
  print(" ImagIntHoenl: " + "{: .8e}".format(imag_integral_a_la_hoenl) + " ImagPrefIntHoenl: " + "{: .8e}".format(imag_integral_a_la_hoenl*prefactor))

  alpha_K_sugiura_damp = complex(integral*prefactor,imag_integral*prefactor)
  alpha_K_hoenl_damp = complex(integral_a_la_hoenl_damp*prefactor,imag_integral_a_la_hoenl*prefactor)
  alpha_K_Hoenl_27 = complex(hoenl_27a(nu_k,nu,delta_K),hoenl_27b(nu_k,nu,delta_K))
  print("Hoenl 27a: " +"{: .8e}".format(alpha_K_Hoenl_27.real))
  print("Hoenl 27b: " +"{: .8e}".format(alpha_K_Hoenl_27.imag))

  print(" Sugiura:      " + "{: .8e}".format(alpha_K_sugiura_damp.real)  + " " + "{:+.8e}".format(alpha_K_sugiura_damp.imag) + " i"
        +"\n Hoenl_damp:   " + "{: .8e}".format(alpha_K_hoenl_damp.real) + " " + "{:+.8e}".format(alpha_K_hoenl_damp.imag) + " i"
        +"\n Hoenl(27a,b): " + "{: .8e}".format(alpha_K_Hoenl_27.real)   + " " + "{:+.8e}".format(alpha_K_Hoenl_27.imag) + " i")
  x = [integral*prefactor,integral_a_la_hoenl_damp*prefactor,hoenl_27a(n_0,nu,0)]
  y = [imag_integral*prefactor,imag_integral_a_la_hoenl*prefactor,hoenl_27b(n_0,nu,0)]
  #plt.scatter(x,y,marker="*")
  #plt.show()
  #fourpithirdalpha = 4.0*math.pi/3.0*alpha_K_sugiura_damp-1
  ausdruck = (12.0*math.pi/3.0*alpha_K_sugiura_damp)/(1-4.0*math.pi/3*alpha_K_sugiura_damp)
  ausdruck_Hoenl_damp = (12.0*math.pi/3.0*alpha_K_hoenl_damp)/(1-4.0*math.pi/3*alpha_K_hoenl_damp)
  ausdruck_Hoenl = (12.0*math.pi/3.0*alpha_K_Hoenl_27)/(1-4.0*math.pi/3*alpha_K_Hoenl_27)
  print("Halber Ausdruck: " + "{:.8e}".format(ausdruck/2.0))
  print("2 Pi Xhi: " + "{:.8e}".format(2*math.pi*alpha_K_sugiura_damp.real))
  #lamb = speed_of_light/nu
  factor = 2*math.pi*el_mass*pow(nu,2)/pow(el_charge,2)
  print("Sugiura curve:    " + "{:.8e}".format(-ausdruck.real/2.0*factor))
  print("Hoenl_damp curve: " + "{:.8e}".format(-ausdruck_Hoenl_damp.real/2.0*factor))
  print("Hoenl_27 curve:   " + "{:.8e}".format(-ausdruck_Hoenl.real/2.0*factor))
  #print("n-1 (that is the complex refractive index):")
  #n_minus1 = (-2*fourpithirdalpha+numpy.sqrt(4*pow(fourpithirdalpha,2)-16*math.pi*alpha_K_sugiura_damp*fourpithirdalpha))/2*fourpithirdalpha
  #n_minus1_2 = (-2*fourpithirdalpha-numpy.sqrt(4*pow(fourpithirdalpha,2)-16*math.pi*alpha_K_sugiura_damp*fourpithirdalpha))/2*fourpithirdalpha
  #print(n_minus1,n_minus1_2)

def hoenl_like_k(Z=None, ener=8047.8):
  # Z is integer number of 
  if type(Z) != type(int(20)):
    print("Z MUST BE INTEGER!")
    return
  if Z == None:
    Z = 6
  Z_s_sq = pow(get_Zeff_1s(Z),2)
  e_ion = get_ionization_energy_1s(Z)
  nu_k = e_ion / h
  delta_K = 1 + pow(alpha/2,2) * Z_s_sq - e_ion/(Ryd_ener * Z_s_sq)
  n_0 = nu_k/(1-delta_K)
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  epsilon = 1E-12
  epsilon2 = 0.001
  x_j = nu*1E-7
  integral1 = integrate.quad(integrand_damped_abs,nu_k,(1-epsilon)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
  integral2 = integrate.quad(integrand_damped_abs,(1+epsilon)*nu,1000*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
  integral = integral1[0] + integral2[0]
  integral_a_la_hoenl_damp = integrate.quad(integrand_a_la_hoenl_damp,nu_k,(1-epsilon)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)[0] \
    + integrate.quad(integrand_a_la_hoenl_damp,(1+epsilon)*nu,1000*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)[0]
  lower_limit = n_0
  if (1-epsilon2)*nu > n_0:
    lower_limit = (1-epsilon2)*nu
  imag_integral1 = integrate.quad(imag_integrand_damped_abs,lower_limit,(1+epsilon2)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
  imag_integral = imag_integral1[0]
  imag_integral_a_la_hoenl = integrate.quad(imag_integrand_damped_a_la_hoenl,lower_limit,(1+epsilon2)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)[0]

  alpha_K_sugiura_damp = complex(integral*prefactor,imag_integral*prefactor)
  alpha_K_hoenl_damp = complex(integral_a_la_hoenl_damp*prefactor,imag_integral_a_la_hoenl*prefactor)
  alpha_K_Hoenl_27 = complex(hoenl_27a(nu_k,nu,delta_K),hoenl_27b(nu_k,nu,delta_K))
  ausdruck = (12.0*math.pi/3.0*alpha_K_sugiura_damp)/(1-4.0*math.pi/3*alpha_K_sugiura_damp)
  ausdruck_Hoenl_damp = (12.0*math.pi/3.0*alpha_K_hoenl_damp)/(1-4.0*math.pi/3*alpha_K_hoenl_damp)
  ausdruck_Hoenl = (12.0*math.pi/3.0*alpha_K_Hoenl_27)/(1-4.0*math.pi/3*alpha_K_Hoenl_27)
  factor = 2*math.pi*el_mass*pow(nu,2)/pow(el_charge,2)
  return -ausdruck.real/2.0*factor, -ausdruck_Hoenl_damp.real/2.0*factor, -ausdruck_Hoenl.real/2.0*factor

def hoenl_like_k_with_imag(Z=None, ener=8047.8):
  # Z is integer number of 
  if type(Z) != type(int(20)):
    print("Z MUST BE INTEGER!")
    return
  if Z == None:
    Z = 6
  Z_s_sq = pow(get_Zeff_1s(Z),2)
  e_ion = get_ionization_energy_1s(Z)
  nu_k = e_ion / h
  delta_K = 1 + pow(alpha/2,2) * Z_s_sq - e_ion/(Ryd_ener * Z_s_sq) #21a in Hoenl
  Z_s_sq_Hoenl = pow(Z-5.0/16.0,2)
  delta_K_Hoenl = 1 + pow(alpha/2,2) * Z_s_sq_Hoenl - e_ion/(Ryd_ener * Z_s_sq_Hoenl)
  n_0 = nu_k/(1-delta_K) #22 in Hoenl
  n_0_Hoenl = nu_k/(1-delta_K_Hoenl)
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  epsilon2 = 0.005
  x_j = nu_k*1E-6
  #integral_test = integrate.quad(integrand_damped_abs,nu_k,50*nu,args=(x_j,nu,n_0),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  #integrating directly as one intergal with the singulary position as point does not work.

  #REAL PART
  #When telling hte quad alorithm that nu is a "position to be carefull with" the results don't require epsilon
  inte1 = integrate.quad(integrand_damped_abs,nu_k,nu,args=(x_j,nu,n_0),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2 = integrate.quad(integrand_damped_abs,nu,20*nu,args=(x_j,nu,n_0),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  if inte1[1] > 0.5*abs(inte1[0]):
    print("!!! inaccurate REAL1 integral at nu_k/nu = " + "{:8.4f}".format(nu_k/nu) + " " + str(inte1) + str(inte1[1]/abs(inte1[0])))
  if inte2[1] > 0.5*abs(inte2[0]):
    print("!!! inaccurate REAL2 integral at nu_k/nu = " + "{:8.4f}".format(nu_k/nu) + " " + str(inte2) + str(inte2[1]/abs(inte2[0])))
  integral = inte1[0] + inte2[0]
  integral_a_la_hoenl_damp = integrate.quad(integrand_a_la_hoenl_damp,nu_k,nu,args=(x_j,nu,n_0_Hoenl),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)[0] \
    + integrate.quad(integrand_a_la_hoenl_damp,nu,20*nu,args=(x_j,nu,n_0_Hoenl),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)[0]
  #IMAG PART
  lower_limit = nu_k
  if (1-epsilon2)*nu > nu_k:
    lower_limit = (1-epsilon2)*nu
  if (1+epsilon2)*nu < nu_k:
    upper_limit = 10*nu_k
  else:
    upper_limit = (1+epsilon2)*nu
  imag_integral = integrate.quad(imag_integrand_damped_abs,lower_limit,upper_limit,args=(x_j,nu,n_0),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10)
  imag_integral_a_la_hoenl = integrate.quad(imag_integrand_damped_a_la_hoenl,lower_limit,upper_limit,points=nu,args=(x_j,nu,n_0_Hoenl),limit=200000,epsabs=1E-55,epsrel=1E-10)
  if imag_integral[1] > 0.5*abs(imag_integral[0]):
    print("!!! inaccurate IMAG integral at nu_k/nu = " + "{:6.4f}".format(nu_k/nu) + " " + str(imag_integral) + " " + str(imag_integral[1]/abs(imag_integral[0])) + " lower_limit is with epsilon?: " + str((1-epsilon2)*nu == lower_limit))

  #this yields the number of dispersion electrons, which is independant of the wavelenght used
  #n_disp_el = integrate.quad(integrand_for_disp_els,nu_k,2000*nu_k,args=(n_0),limit=200000,epsabs=1E-60)[0]
  #print(n_disp_el)

  alpha_K_sugiura_damp = complex(integral,imag_integral[0])
  alpha_K_hoenl_damp = complex(integral_a_la_hoenl_damp,imag_integral_a_la_hoenl[0])
  ausdruck = alpha_K_sugiura_damp/(1-4.0*math.pi/3*alpha_K_sugiura_damp*prefactor)
  ausdruck_Hoenl_damp = alpha_K_hoenl_damp/(1-4.0*math.pi/3*alpha_K_hoenl_damp*prefactor)
  factor = pow(nu,2)
  return -ausdruck.real*factor, -ausdruck_Hoenl_damp.real*factor, -ausdruck.imag*factor, -ausdruck_Hoenl_damp.imag*factor#, n_disp_el

def hoenl_27_with_imag(Z=None, ener=8047.8):
  # Z is integer number of 
  if type(Z) != type(int(20)):
    print("Z MUST BE INTEGER!")
    return
  if Z == None:
    Z = 6
  Z_s_sq = pow(Z-5.0/16.0,2)
  e_ion = get_ionization_energy_1s(Z)
  nu_k = e_ion / h
  delta_K = 1 + pow(alpha/2,2) * Z_s_sq - e_ion/(Ryd_ener * Z_s_sq) #21a in Hoenl
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  alpha_K_Hoenl_27 = complex(hoenl_27a(nu_k,nu,delta_K),hoenl_27b(nu_k,nu,delta_K))
  ausdruck = (2.0*math.pi*alpha_K_Hoenl_27)/(1-4.0*math.pi/3*alpha_K_Hoenl_27)
  factor = 2*math.pi*el_mass*pow(nu,2)/pow(el_charge,2)
  return -ausdruck.real/2.0*factor, -ausdruck.imag/2.0*factor

def hoenl_27_with_imag_and_real_Z(Z=None, ener=8047.8):
  # Z is integer number of 
  if type(Z) != type(int(20)):
    print("Z MUST BE INTEGER!")
    return
  if Z == None:
    Z = 6
  Z_s_sq = pow(get_Zeff_1s(Z),2)
  e_ion = get_ionization_energy_1s(Z)
  nu_k = e_ion / h
  delta_K = 1 + pow(alpha/2,2) * Z_s_sq - e_ion/(Ryd_ener * Z_s_sq) #21a in Hoenl
  nu = ener / h #Cu alpha: 8047.8 eV, molly 17450 eV
  alpha_K_Hoenl_27 = complex(hoenl_27a(nu_k,nu,delta_K),hoenl_27b(nu_k,nu,delta_K))
  ausdruck = (2.0*math.pi*alpha_K_Hoenl_27)/(1-4.0*math.pi/3*alpha_K_Hoenl_27)
  factor = 2*math.pi*el_mass*pow(nu,2)/pow(el_charge,2)
  return -ausdruck.real/2.0*factor, -ausdruck.imag/2.0*factor
