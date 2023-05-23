import numpy as np
import math
from constants_and_atomic_properties import *
from matrix_coefficients_v2 import *
from hoenl_like import *
import scipy.integrate as integrate
import scipy.special as special

def integrand_a_la_hoenl(n_j,nu,n_0):
  var = (n_j/n_0)
  part1 = 128/(3.0*n_0)
  part2 = np.exp(-4)/3.0*(4*var-1)
  part3 = 1/pow(var,4)

  res = 1/(pow(n_j,2)-pow(nu,2))
  return res*part1*part2*part3

def integrand_a_la_hoenl_damp(n_j,chi_j,nu,n_0):
  var = (n_j/n_0)
  part1 = 128/(3.0*n_0)
  part2 = np.exp(-4)/3.0*(4*var-1)
  part3 = 1/pow(var,4)
  res = (pow(n_j,2)-pow(nu,2))/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*part1*part2*part3

def integrand(n_j,nu,n_0):
  var = (n_j/n_0)-1
  sqrt_var = np.sqrt(var)
  part1 = 128*pow(n_0,3)/(3*pow(n_j,4))
  part2 = np.exp(-4/sqrt_var*np.arctan(sqrt_var))
  part3 = 1-np.exp(-2*math.pi/sqrt_var)
  res = 1/(pow(n_j,2)-pow(nu,2))
  return res*part1*part2/part3

def integrand_abs(n_j,nu,n_0):
  var = abs((n_j/n_0)-1)
  sqrt_var = np.sqrt(var)
  part1 = 128*pow(n_0,3)/(3*pow(n_j,4))
  part2 = np.exp(-4/sqrt_var*np.arctan(sqrt_var))
  part3 = 1-np.exp(-2*math.pi/sqrt_var)
  res = 1/(pow(n_j,2)-pow(nu,2))
  return res*part1*part2/part3

def integrand_damped(n_j,chi_j,nu,n_0):
  var = (n_j/n_0)-1
  sqrt_var = np.sqrt(var)
  part1 = 128*pow(n_0,3)/(3*pow(n_j,4))
  part2 = np.exp(-4/sqrt_var*np.arctan(sqrt_var))
  part3 = 1-np.exp(-2*math.pi/sqrt_var)
  res = (pow(n_j,2)-pow(nu,2))/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*part1*part2/part3

def imag_integrand_damped(n_j,chi_j,nu,n_0):
  var = (n_j/n_0)-1
  sqrt_var = np.sqrt(var)
  part1 = 128*pow(n_0,3)/(3*pow(n_j,4))
  part2 = np.exp(-4/sqrt_var*np.arctan(sqrt_var))
  part3 = 1-np.exp(-2*math.pi/sqrt_var)
  res = -nu*chi_j/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*part1*part2/part3

def imag_integrand_damped_a_la_hoenl(n_j,chi_j,nu,n_0):
  var = (n_j/n_0)
  part1 = 128/(3.0*n_0)
  part2 = np.exp(-4)/3.0*(4*var-1)
  part3 = 1/pow(var,4)

  res = -nu*chi_j/(pow(pow(n_j,2)-pow(nu,2),2)+pow(nu*chi_j,2))
  return res*part1*part2*part3

def oszillator_sugiura_K(z):
  if(z<1):
    var = np.sqrt(abs(z-1))
    zaehler = np.exp(-4/var*np.arctan(var)) - 2*1/3.0*(z-1)*special.hyp2f1(0.75,1,1.75,(z-1)*(z-1))
  else:
    var = np.sqrt(z-1)
    zaehler = np.exp(-4/var*np.arctan(var))
  nenner = 1-np.exp(-2*np.pi/var)
  return zaehler/nenner

def hoenl_26b(nu,n_0,n_k):
  if nu >= n_k: 
    oscillator = oszillator_sugiura_K(nu/n_0)
    #factor_26b = np.pi*2*pow(el_charge,2)/(4*pow(np.pi,2)*el_mass*2*nu)
    #factor_alpha = 2*np.pi
    #factor_plot = 2*math.pi*el_mass*pow(nu,2)/pow(el_charge,2)
    #all_factor = factor_26b*factor_alpha*factor_plot/n_0
    factor = math.pi*nu/n_0
    return oscillator*factor
  else:
    return 0.0

def hoenl_27a(nu_K,nu,delta_K):
  x = nu_K/nu
  term1 = -32*pow(el_charge,2)*np.exp(-4)/(9*pow(math.pi,2)*el_mass*pow(nu,2))
  term2 = 4/pow(1-delta_K,2)*(1+pow(x,2)*math.log(abs(1-1/pow(x,2))))
  term3 = -1/pow(1-delta_K,3)*(2.0/3.0 + 2*pow(x,2) + pow(x,3)*math.log(abs((1-x)/(1+x))))
  result = term1*(term2+term3)
  return result

def hoenl_27b(nu_K,nu,delta_K):
  x = nu_K/nu
  if x > 1:
    return 0.0
  term1 = -32*pow(el_charge,2)*np.exp(-4)/(9*math.pi*el_mass*pow(nu,2))
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
  #x = np.arange(1.01*nu_0,(1-epsilon)*nu,1E16)
  #y = integrand(x)
  #plt.plot(x,y)
  #plt.show()
  epsilon2 = 0.001
  x_j = nu*1E-8
  #x = np.arange((1-0.0001)*nu,(1+0.0001)*nu,1E12)
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
  integral1 = integrate.quad(real_general_integration,nu_k,(1-epsilon)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
  integral2 = integrate.quad(real_general_integration,(1+epsilon)*nu,1000*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
  integral = integral1[0] + integral2[0]
  print(" Int:          " + "{: .8e}".format(integral) +                 " PrefInt:          " + "{: .8e}".format(integral*prefactor))
  integral_a_la_hoenl_damp = integrate.quad(integrand_a_la_hoenl_damp,nu_k,(1-epsilon)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)[0] \
    + integrate.quad(integrand_a_la_hoenl_damp,(1+epsilon)*nu,1000*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)[0]
  print(" IntHoenl:     " + "{: .8e}".format(integral_a_la_hoenl_damp) + " PrefIntHoenl:     " + "{: .8e}".format(integral_a_la_hoenl_damp*prefactor))
  lower_limit = n_0
  if (1-epsilon2)*nu > n_0:
    lower_limit = (1-epsilon2)*nu
  imag_integral1 = integrate.quad(imag_general_integration,lower_limit,(1+epsilon2)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
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
  #n_minus1 = (-2*fourpithirdalpha+np.sqrt(4*pow(fourpithirdalpha,2)-16*math.pi*alpha_K_sugiura_damp*fourpithirdalpha))/2*fourpithirdalpha
  #n_minus1_2 = (-2*fourpithirdalpha-np.sqrt(4*pow(fourpithirdalpha,2)-16*math.pi*alpha_K_sugiura_damp*fourpithirdalpha))/2*fourpithirdalpha
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
  integral1 = integrate.quad(real_general_integration,nu_k,(1-epsilon)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
  integral2 = integrate.quad(real_general_integration,(1+epsilon)*nu,1000*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
  integral = integral1[0] + integral2[0]
  integral_a_la_hoenl_damp = integrate.quad(integrand_a_la_hoenl_damp,nu_k,(1-epsilon)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)[0] \
    + integrate.quad(integrand_a_la_hoenl_damp,(1+epsilon)*nu,1000*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)[0]
  lower_limit = n_0
  if (1-epsilon2)*nu > n_0:
    lower_limit = (1-epsilon2)*nu
  imag_integral1 = integrate.quad(imag_general_integration,lower_limit,(1+epsilon2)*nu,args=(x_j,nu,n_0),limit=200000,epsabs=1E-60)
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
  inte1 = integrate.quad(real_general_integration,nu_k,nu,args=(x_j,nu,n_0),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
  inte2 = integrate.quad(real_general_integration,nu,20*nu,args=(x_j,nu,n_0),points=nu,limit=200000,epsabs=1E-50,epsrel=1E-10)
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
  imag_integral = integrate.quad(imag_general_integration,lower_limit,upper_limit,args=(x_j,nu,n_0),points=nu,limit=200000,epsabs=1E-55,epsrel=1E-10)
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
  conjugate_function = dummy_func
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
  conjugate_function = dummy_func
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
  conjugate_function = dummy_func
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
  conjugate_function = dummy_func
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
    
def A_l_from_z(b_, z, n_0, l, nu, p_limit):
  part1 = b_/2/pow(-2*b_*math.sqrt(z-1),l+1)
  sum = 0
  for p in range(p_limit):
    n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
    J1 = J(1,p,1,l,W11)
    if (J == 0):
      continue
    K1 = K_recursive_from_z(p,l,b_,z, n_0)
    sum += n1 * J1 * K1
  return part1 * sum

def C_l_from_z(b_, z, n_0, l, nu, p_limit):
  part1 = b_/pow(-2*b_*math.sqrt(z-1),l+1)
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
  
def B2_from_z(b_, z, n_0, l, nu, p_limit):
  part1 = b_*b_/(4*pow(-2*b_*math.sqrt(z-1),l+1))
  sum = 0
  for p in range(p_limit):
    n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
    J1 = J(2,p,2,l,W22)
    if (J1 == 0):
      continue
    K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
    sum += n1 * J1 * K1
  return part1 * sum

def B1_from_z(b_, z, n_0, l, nu, p_limit):
  part1 = b_*b_/(2*pow(-2*b_*math.sqrt(z-1),l+1))
  sum = 0
  for p in range(p_limit):
    n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
    J1 = J(1,p+1,1,l,W11)
    if J == 0:
      continue
    K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
    sum += n1 * J1 * K1
  return part1 * sum
  
def B0_from_z(b_, z, n_0, l, nu, p_limit):
  part1 = b_/pow(-2*b_*math.sqrt(z-1),l+1)
  sum = 0
  for p in range(p_limit):
    n1 = pow(complex(0,-q(nu)),p) / math.factorial(p)
    J1 = 0.25 * J(2,p,0,l, W20)
    J2 = J(0,p,0,l, W00)
    if J1 == 0 and J2 == 0:
      continue
    K1 = K_recursive_from_z(p+1,l,b_,z,n_0)
    K2 = K_recursive_from_z(p,l,b_,z,n_0)
    sum += n1 * (2*b_*J1 * K1 - J2 * K2)
  return part1 * sum

def J2(p,l):
  if p == 0:
    return 16/5 * kronecke_delta(l,2)
  elif p == 1:
    return 16/7 * kronecke_delta(l,3)
  elif p == 2:
    return 16/35 * kronecke_delta(l,2) + 32/21 * kronecke_delta(l,4)
  elif p == 3:
    return 16/21 * kronecke_delta(l,3) + 32/33 * kronecke_delta(l,5)
  else:
    return J(2,p,2,l,W22)
def J1(p,l):
  if p == 0:
    return 4/3 * kronecke_delta(l,1)
  elif p == 1:
    return 4/5 * kronecke_delta(l,2)
  elif p == 2:
    return 4/15 * kronecke_delta(l,1) + 16/35 * kronecke_delta(l,3)
  elif p == 3:
    return 12/35 * kronecke_delta(l,2) + 16/63 * kronecke_delta(l,4)
  elif p == 4:
    return 4/35 * kronecke_delta(l,1) + 32/105 * kronecke_delta(l,3) + 32/231 * kronecke_delta(l,5)
  else:
    return J(1,p,1,l,W11)
def J0_bar(p,l):
  if p == 0:
    return 2 * kronecke_delta(l,0)
  elif p == 1:
    return 2/3 * kronecke_delta(l,1)
  elif p == 2:
    return 2/3 * kronecke_delta(l,0) + 4/15 * kronecke_delta(l,2)
  elif p == 3:
    return 2/5 * kronecke_delta(l,1) + 4/35 * kronecke_delta(l,3)
  elif p == 4:
    return 2/5 * kronecke_delta(l,0) + 8/35 * kronecke_delta(l,2) + 16/315 * kronecke_delta(l,4)
  elif p == 5:
    return 2/7 * kronecke_delta(l,1) + 8/63 * kronecke_delta(l,3) + 16/693 * kronecke_delta(l,5)
  else:
    return J(0,p,0,l,W00)
def J0_hat(p,l):
  if p == 0:
    return 1/3 * kronecke_delta(l,0) - 1/15 * kronecke_delta(l,2)
  elif p == 1:
    return 1/15 * kronecke_delta(l,1) - 1/35 * kronecke_delta(l,3)
  elif p == 2:
    return 1/15 * kronecke_delta(l,0) + 1/105 * kronecke_delta(l,2) - 4/315 * kronecke_delta(l,4)
  elif p == 3:
    return 1/35 * kronecke_delta(l,1) - 1/315 * kronecke_delta(l,3) - 4/693 * kronecke_delta(l,5)
  else:
    return 0.25 * J(2,p,0,l, W20)
  
def test_florian() -> None:
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