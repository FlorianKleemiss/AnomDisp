from legendre_plynomials import *
from matrix_coefficients_v2 import *
import matplotlib.pyplot as plt
from brennan import brennan
import scipy.integrate as integrate
from scipy.special import factorial2
from mpl_toolkits import mplot3d

def getting_agreement():
  import random
  random.seed()
  alpha = math.pi/4
  theta0 = np.linspace(0,math.pi,100)
  y_l = []
  y_l2 = []
  y_r = []
  y_parts = [[] for x in range(15)]
  labels = [
  r'$\alpha_{10}^1$',
  r'$\alpha_{11}^1$',
  r'$\alpha_{12}^1$',
  r'$\alpha_{10}^0$',
  r'$\alpha_{10}^2$',
  r'$\alpha_{11}^0$',
  r'$\alpha_{11}^2$',
  r'$\alpha_{12}^0$',
  r'$\alpha_{12}^2$',
  r'$\beta_{12}^0$',
  r'$\beta_{12}^2$',
  r'$\bar \alpha_{10}^2$',
  r'$\bar \alpha_{11}^2$',
  r'$\bar \alpha_{12}^2$',
  r'$\bar \beta_{12}^2$'
  ]
  sa = math.sin(alpha)
  ca = math.cos(alpha)
  for theta in theta0:
    ct0 = math.cos(theta)
    st0 = math.sin(theta)
    #nomen clature a_l_m'_m = alpha_coef(l,m,m')
    # m'=0
    a_1_0_1 = alpha_coef(1,1,0,theta,alpha)
    a_1_0_0 = alpha_coef(1,0,0,theta,alpha)
    a_1_0_2 = alpha_coef(1,2,0,theta,alpha)
    # m'=1
    a_1_1_0 = alpha_coef(1,0,1,theta,alpha)
    a_1_1_1 = alpha_coef(1,1,1,theta,alpha)
    a_1_1_2 = alpha_coef(1,2,1,theta,alpha)
    # m'=2
    a_1_2_0 = alpha_coef(1,0,2,theta,alpha)
    a_1_2_1 = alpha_coef(1,1,2,theta,alpha)
    a_1_2_2 = alpha_coef(1,2,2,theta,alpha)
    #beta
    b_1_2_0 = beta_coef (1,0,2,theta,alpha)
    b_1_2_2 = beta_coef (1,2,2,theta,alpha)
    ab_1_0_2 = alpha_bar_coef(1,2,0,theta,alpha)
    ab_1_1_2 = alpha_bar_coef(1,2,1,theta,alpha)
    ab_1_2_2 = alpha_bar_coef(1,2,2,theta,alpha)
    bb_1_2_2 = beta_bar_coef(1,2,2,theta,alpha)
    y_parts[0].append( a_1_0_1 )#*st0)
    y_parts[1].append( a_1_1_1 )#*ct0)
    y_parts[2].append( a_1_2_1 )#*st0)
    y_parts[3].append( a_1_0_0 )#*-ct0*ca)
    y_parts[4].append( a_1_0_2 )#*-ct0*ca)
    y_parts[5].append( a_1_1_0 )#*st0*ca)
    y_parts[6].append( a_1_1_2 )#*st0*ca)
    y_parts[7].append( a_1_2_0 )#*-ct0*ca)
    y_parts[8].append( a_1_2_2 )#*-ct0*ca)
    y_parts[9].append( b_1_2_0 )#*sa)
    y_parts[10].append(b_1_2_2 )#*sa)
    y_parts[11].append(ab_1_0_2)#*-ct0*sa)
    y_parts[12].append(ab_1_1_2)#*st0*sa)
    y_parts[13].append(ab_1_2_2)#*-ct0*sa)
    y_parts[14].append(bb_1_2_2)#*-ca)
    
    left_side_ =       a_1_0_1*st0     + a_1_1_1*ct0    + a_1_2_1*st0\
    -a_1_0_0*ct0*ca  - a_1_0_2*ct0*ca  + a_1_1_0*st0*ca + a_1_1_2*st0*ca\
    +b_1_2_0*sa      + b_1_2_2*sa      - a_1_2_0*ct0*ca - a_1_2_2*ct0*ca\
    -ab_1_0_2*ct0*sa + ab_1_1_2*st0*sa - bb_1_2_2*ca    - ab_1_2_2*ct0*sa
    #left2 = 0
    #for i in range(15):
    #  left2 += y_parts[i][-1]
    
    right_side_ = 3 * a_1_1_1
    y_l.append(left_side_)
    #y_l2.append(left2)
    y_r.append(right_side_)
  
  fig = plt.figure()
  axes = fig.add_subplot(1,1,1)
  #axes.plot(theta0,y_r,label=r'$\alpha_{11}^1$')
  axes.scatter(theta0,y_l,s=10,facecolors='none',edgecolors='b',label="Sum")
  #plt.scatter(theta0,y_l2,s=10,facecolors='none',edgecolors='r',label="sum -1")
  formats=[
    ":+","-+","-.+", #m'=1
    ":o","-o", #m=0
    ":x","-x", #m=1
    ":8","-8", #m=2
    ":*","-*", #beta
    ":^","-.^","-^","--^" #bars
  ]
  for i in range(15):
    axes.plot(theta0,y_parts[i],formats[i],label=labels[i])
  axes.legend()
  fig.tight_layout()
  fig.subplots_adjust(left=0.03,bottom=0.03,top=1,right=1)
  mng = plt.get_current_fig_manager()
  mng.window.state('zoomed')
  plt.show()
  
  print("left side = " + str(left_side_))
  print("right side = " + str(right_side_))
  print("ratio = " + str(left_side_/right_side_))
  print("happy? " + str(right_side_ == left_side_))

def Psi_parallel_a(l,m,k,theta0, alpha):
  return alpha_coef(l,1,1,theta0,alpha)

def Psi_parallel_b(l,m,k, theta0, alpha):
  a = alpha_coef(l,m,k,theta0,alpha) 
  if k == 0:
    return a * math.sin(theta0)
  elif k == 1:
    return a * math.cos(theta0)
  elif k == 2:
    return a * math.sin(theta0)

def Psi_parallel_c(l,m,k, theta0, alpha):
  a = alpha_coef(l,m,k,theta0,alpha)
  if k == 0:
    return -a * math.cos(theta0) * math.cos(alpha)
  elif k == 1:
    return a * math.sin(theta0) * math.cos(alpha)
  elif k == 2:
    return beta_coef(l,m,k,theta0,alpha) * math.sin(alpha) \
      - a * math.cos(theta0) * math.cos(alpha)

def Psi_parallel_d(l,m,k,theta0,alpha):
  a = alpha_bar_coef(l,m,k,theta0,alpha)
  if k == 0:
    return -a * math.cos(theta0) * math.sin(alpha)
  elif k == 1:
    return a * math.sin(theta0) * math.sin(alpha)
  elif k == 2:
    return -beta_bar_coef(l,m,k,theta0,alpha) * math.cos(alpha) \
      - a * math.cos(theta0) * math.sin(alpha)

def Psi_orth_a(l,m,k,theta0, alpha):
  return beta_coef(l,1,1,theta0,alpha)

def Psi_orth_b(l,m,k, theta0, alpha):
  b = beta_coef(l,m,k,theta0,alpha)
  if k == 0:
    return b * math.sin(theta0)
  elif k == 1:
    return b * math.cos(theta0)
  elif k == 2:
    return b * math.sin(theta0)

def Psi_orth_c(l,m,k, theta0, alpha):
  a = alpha_coef(l,m,k,theta0,alpha)
  b = beta_coef(l,m,k,theta0,alpha)
  if k == 0:
    return a * math.sin(alpha)
  elif k == 1:
    return b * math.sin(theta0) * math.cos(alpha)
  elif k == 2:
    return -a * math.sin(alpha) - b * math.cos(theta0) * math.cos(alpha)

def Psi_orth_d(l,m,k,theta0,alpha):
  a = alpha_bar_coef(l,m,k,theta0,alpha)
  b = beta_bar_coef(l,m,k,theta0,alpha)
  if k == 0:
    return a * math.cos(alpha)
  elif k == 1:
    return b * math.sin(theta0) * math.sin(alpha)
  elif k == 2:
    return a * math.cos(alpha) - b * math.cos(theta0) * math.sin(alpha)

# n is the complex quantum number where we are in energy spectrum for integration
# z2 is the value of z at the absorption edge of the element we look at
# el_nr is the number of the element in the periodic table
# E_in is the input energy (incoming X-ray)
def calc_S_values(theta0, alpha, l_max, p_max, z, z2, nu_in, el_nr):
  a=0
  b=0
  c=0
  d=0
  for l in range(0,l_max+1):
    for k in range(0,min(l,2)):
      c+= Psi_parallel_c(l,0,k,theta0,alpha) * f_c(el_nr,l,k,z,z2,nu_in,p_max)
      if l == 0:
        continue
      a+=Psi_parallel_a(l,1,k,theta0,alpha) * f_a(el_nr,l,k,z,z2,nu_in,p_max)
      b+=Psi_parallel_b(l,0,k,theta0,alpha) * f_b(el_nr,l,k,z,z2,nu_in,p_max)
      if l == 1:
        continue
      c+=Psi_parallel_c(l,2,k,theta0,alpha) * f_c(el_nr,l,k,z,z2,nu_in,p_max)
      d+=Psi_parallel_d(l,2,k,theta0,alpha) * f_d(el_nr,l,k,z,z2,nu_in,p_max)
  return a,b,c,d

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
  part1 = 512*(z+8.0/3.0)/(3*pow(z,5))
  part2 = sugiura_exps(z,2)
  return part1*part2

def f_p_2_1_EM(z):
  part1 = 512*(z-2)*(z+8.0/3.0)/(15*pow(z,7))
  part2 = sugiura_exps(z,2)
  return part1*part2

def f_p_2_2_EM(z):
  part1 = 1048*(11*z-6)*(z+3)/(45*pow(z,7))
  part2 = sugiura_exps(z,2)
  return part1*part2

def f_p_1(z):
  part1 = 512*(z+8/3)/(3*pow(z,5))
  z_1=z-1
  sqrt_var = numpy.sqrt(z_1)
  part2 = numpy.exp(-8/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-4*math.pi/sqrt_var)
  return part1*part2/part3

def f_l0(z):
  part1 = 512/(9*pow(z,4))
  z_1=z-1
  sqrt_var = numpy.sqrt(z_1)
  part2 = numpy.exp(-8/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-4*math.pi/sqrt_var)
  return part1*part2/part3

def f_l2(z):
  part1 = 2*8192*(z+3)/(45*pow(z,5))
  z_1=z-1
  sqrt_var = numpy.sqrt(z_1)
  part2 = numpy.exp(-8/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-4*math.pi/sqrt_var)
  return part1*part2/part3

def f_l2_c0(z):
  part1 = -pow(2,14)*(z+3)/(45*pow(z,5))
  z_1=z-1
  sqrt_var = numpy.sqrt(z_1)
  part2 = numpy.exp(-8/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-4*math.pi/sqrt_var)
  return part1*part2/part3

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

def x2(nu_in,n_0, el_nr, l_0):
  return pow(n_0* q(nu_in)/b(n_0, l_0, el_nr), 2)

def xn(nu_in,n_0, el_nr, l_0, n):
  return pow(n_0* q(nu_in)/b(n_0, l_0, el_nr), n)

a_result = []
b_result = []
c_result = []
d_result = []
e_result = []
f_result = []
g_result = []
h_result = []
i_result = []
j_result = []
k_result = []
l_result = []
constant_factor = 4*math.pi*math.pi*el_mass/h/h
br = brennan()

def test_for_EM_and_Hoenl(el_nr, z, nu_in):
  #temp_voll = f_a(
  #  Z=el_nr,
  #  l=1,
  #  z=z,
  #  nu_in = nu_in,
  #  n_0 = 1,
  #  p_limit = 3)
##
  #temp1 = f_s_1_hoenl(z) * constant_factor
  #g_result.append(temp1)
  #h_result.append(temp_voll)
  #temp2 = f_s_2_hoenl(z) * constant_factor * x2(nu_in, 1, el_nr, 0)
  #temp_voll_2 = f_a(
  #  Z=el_nr,
  #  l=2,
  #  z=z,
  #  nu_in = nu_in,
  #  n_0 = 1,
  #  p_limit = 5
  #) 
  #e_result.append(temp2)
  #temp_voll_3 = f_a(
  #  Z=el_nr,
  #  l=3,
  #  z=z,
  #  nu_in = nu_in,
  #  n_0 = 1,
  #  p_limit = 5
  #) 
  #f_result.append(temp_voll_2)  

  #temp3 = f_s_1(z) * constant_factor
  #temp_voll_3 = f_a(
  #  Z=el_nr,
  #  l=1,
  #  z=z,
  #  nu_in = nu_in,
  #  n_0 = 2,
  #  p_limit = 3)
  #c_result.append(temp3)
  #d_result.append(temp_voll_3)
#
  #temp4 = f_s_2(z) * constant_factor * x2(nu_in, 2, el_nr, 0)
  temp_voll_4 = f_a(
    Z=el_nr,
    l=1,
    k=0,
    z=z,
    nu_in = nu_in,
    n_0 = 1,
    p_limit = 3) 
  #a_result.append(temp4)
  b_result.append(temp_voll_4)

def test_for_fb_fc_fd(el_nr, z, nu_in):
  a_res = 0
  b_res = 0
  c_res = 0
  d_res = 0
  e_res = 0
  f_res = 0
  g_res = 0
  h_res = 0
  i_res = 0
  p_limit = 3
  #a_res += f_s_1(z) * constant_factor
  #b_res += f_p_1(z) * constant_factor
  
  #c_res += f_l0(z) * constant_factor
  #d_res += f_c_0(el_nr,0,0,z,nu_in,2,1)
  a_res = f_l2_c0(z) * constant_factor
  b_res = -f_l2_c0(z) * constant_factor
  c_res = -0.5*f_l2_c0(z) * constant_factor
  e_res += f_l2(z) * constant_factor
  l = 2
  for k_ in range(0,3):
    f_res += f_c_0(el_nr,l,k_,z,nu_in,2,p_limit)
    g_res += f_b(  el_nr,l,k_,z,nu_in,2,p_limit)
    h_res += f_c_2(el_nr,l,k_,z,nu_in,2,p_limit)
    i_res += f_d(  el_nr,l,k_,z,nu_in,2,p_limit)
  #  if l == 0: continue
  #  c_res += abs(f_b(el_nr,l,1,z,nu_in,2,p_limit))
  #  c_res += abs(f_c_0(el_nr,l,1,z,nu_in,2,p_limit))
  #  c_res += abs(f_c_2(el_nr,l,1,z,nu_in,2,p_limit))
  #  c_res += abs(f_d(el_nr,l,1,z,nu_in,2,p_limit))
  #  if l == 1: continue
  #  c_res += abs(f_b(el_nr,l,2,z,nu_in,2,p_limit))
  #  c_res += abs(f_c_0(el_nr,l,2,z,nu_in,2,p_limit)) 
  #  c_res += abs(f_c_2(el_nr,l,2,z,nu_in,2,p_limit))
  #  c_res += abs(f_d(el_nr,l,2,z,nu_in,2,p_limit))
  #c_res = math.log(c_res)
  d_res = f_res + g_res + h_res + i_res
  a_result.append(a_res)
  b_result.append(b_res)
  
  c_result.append(c_res)
  
  d_result.append(d_res)
  e_result.append(e_res)
  f_result.append(f_res)
  g_result.append(g_res)
  h_result.append(h_res)
  i_result.append(i_res)
def plot_abcd_test():
  fig, axes = plt.subplots(1,1)
  #axes[0].scatter(x,a_result,s=10,facecolors='none',edgecolors='b',label="fs1")
  #axes[0].scatter(x,b_result,s=10,facecolors='none',edgecolors='g',label="fp1")
  #axes[1].scatter(x,c_result,s=10,facecolors='none',edgecolors='r',label="c paper l=0 k=0")
  #axes[1].scatter(x,d_result,s=10,facecolors='none',edgecolors='black',label="c l=0 k=0")
  axes.scatter(x,e_result,s=10,facecolors='none',edgecolors='y',label="l=2 paper")
  axes.scatter(x,a_result,s=30,facecolors='none',edgecolors='y',label="l=2 c_0 paper")
  axes.scatter(x,b_result,s=30,facecolors='none',edgecolors='y',label="l=2 -c_0 paper")
  axes.scatter(x,c_result,s=50,facecolors='none',edgecolors='y',label="l=2 -0.5 c_0 paper")
  axes.scatter(x,f_result,s=10,facecolors='none',edgecolors='b',label="l=2 c_0 code")
  axes.scatter(x,g_result,s=10,facecolors='none',edgecolors='g',label="l=2 b code")
  axes.scatter(x,h_result,s=20,facecolors='none',edgecolors='r',marker='^',label="l=2 c_2 code")
  axes.scatter(x,i_result,s=10,facecolors='none',edgecolors='black',label="l=2 d code")
  axes.plot(x,d_result,color='black',label="l=2 complete code")

  axes.legend()
  plt.show()

def plot_EM_Hoenl_test():
  fig, axes = plt.subplots(1,1)
  #axes.scatter(x,a_result,s=10,facecolors='none',edgecolors='b',marker='^',label="EM l=1")
  axes.scatter(x,b_result,s=30,facecolors='none',edgecolors='g',label="PK l=1")
  #axes[0].scatter(x,c_result,s=10,facecolors='none',edgecolors='r',marker='^',label="EM l=1")
  #axes.scatter(x,d_result,s=30,facecolors='none',edgecolors='y',label="PK l=1")
  
  #axes.scatter(x,e_result,s=10,facecolors='none',edgecolors='b',marker='^',label="Hönl l=2")
  #axes.scatter(x,f_result,s=30,facecolors='none',edgecolors='g',label="PK l=1")
  #axes.scatter(x,g_result,s=10,facecolors='none',edgecolors='r',marker='^',label="Hönl l=1")
  #axes.scatter(x,h_result,s=30,facecolors='none',edgecolors='y',label="PK l=1")
  #axes.plot(x,g_result,'--',label="f_table1")
  axes.legend()
  #axes[1].legend()
  
  plt.show()

def plot_integral_test(edge):
  fig, axes = plt.subplots(1,1)
  
  axes.scatter(x,a_result,s=10,facecolors='none',edgecolors='b',marker='^',label="Hönl l=1")
  #axes.scatter(x,b_result,s=20,facecolors='none',edgecolors='r',label="Hönl good epsilon l=1")
  axes.scatter(x,c_result,s=10,facecolors='none',edgecolors='g',marker='^',label="Hönl complex l=1")
  axes.scatter(x,d_result,color='r',label="Brennan real")
  axes.scatter(x,e_result,color='y',label="Brennan complex")
  axes.vlines(edge,-2,2)
  axes.hlines(0.0,start_nu_in,end_nu_in)
  axes.legend()
  
  plt.show()

def test_integral_hönl(nu_in,el_nr, p_limit):
  z_null = h * nu_in / get_ionization_energy_1s(el_nr)
  z_null_ls = h * nu_in / get_ionization_energy_2s(el_nr)
  z_null_lp1 = h * nu_in / get_ionization_energy_2p1_2(el_nr)
  z_null_lp2 = h * nu_in / get_ionization_energy_2p3_2(el_nr)
  n_0 = 1
  l = 1
  k = 1
  epsilon = 1E-10
  if z_null >= 1:
    integral1 = integrate.quad(
      integrand_matrix_s,
      1,
      z_null-epsilon,
      args=(z_null,el_nr, l, k, nu_in, n_0, p_limit),
      points=z_null,
      limit=200000,
      epsabs=1E-55,
      epsrel=1E-10,
      full_output=2)
    integral2 = integrate.quad(
      integrand_matrix_s,
      z_null+epsilon,
      100,
      args=(z_null,el_nr, l, k, nu_in, n_0, p_limit),
      points=z_null,
      limit=200000,
      epsabs=1E-55,
      epsrel=1E-10,
      full_output=2)
    a_result.append(-2*(integral1[0] + integral2[0])/constant_factor + kpcor[el_nr] - relcor[el_nr])
  else:
    total_integral = integrate.quad(
      integrand_matrix_s,
      1,
      100,
      args=(z_null,el_nr, l, k, nu_in, n_0, p_limit),
      points=z_null,
      limit=200000,
      epsabs=1E-55,
      epsrel=1E-10,
      full_output=2)
    a_result.append(-2*total_integral[0]/constant_factor + kpcor[el_nr] - relcor[el_nr])
    #b_result.append(-2*total_integral[0]/constant_factor)
  if z_null_ls > 1:
    integral1 = integrate.quad(
      integrand_matrix_s,
      1,
      z_null_ls-epsilon,
      args=(z_null_ls,el_nr, l, k, nu_in, 2, p_limit),
      points=z_null_ls,
      limit=200000,
      epsabs=1E-55,
      epsrel=1E-10,
      full_output=2)
    integral2 = integrate.quad(
      integrand_matrix_s,
      z_null_ls+epsilon,
      100,
      args=(z_null_ls,el_nr, l, k, nu_in, 2, p_limit),
      points=z_null_ls,
      limit=200000,
      epsabs=1E-55,
      epsrel=1E-10,
      full_output=2)
    a_result[-1] -= 2*(integral1[0]+integral2[0]) / constant_factor
  else:
    integral = integrate.quad(
      integrand_matrix_s,
      1,
      100,
      args=(z_null_ls,el_nr, 2, k, nu_in, 2, p_limit),
      points=z_null_ls,
      limit=200000,
      epsabs=1E-55,
      epsrel=1E-10,
      full_output=2)
    a_result[-1] -= (integral[0]) / constant_factor
  if z_null_lp1 >1:
    res = 0
    integral1 = integrate.quad(
      integrand_matrix_p,
      1,
      z_null_lp1-epsilon,
      args=(z_null_lp1,el_nr, nu_in, 2, p_limit),
      points=z_null_lp1,
      limit=200000,
      epsabs=1E-55,
      epsrel=1E-10,
      full_output=2)
    integral2 = integrate.quad(
      integrand_matrix_p,
      z_null_lp1+epsilon,
      100,
      args=(z_null_lp1,el_nr, nu_in, 2, p_limit),
      points=z_null_lp1,
      limit=200000,
      epsabs=1E-55,
      epsrel=1E-10,
      full_output=2)
    res += integral1[0]+integral2[0]
    a_result[-1] -= 2 * (res) / constant_factor
  else:
    res = 0
    integral = integrate.quad(
      integrand_matrix_p,
      1,
      100,
      args=(z_null_lp1,el_nr, nu_in, 2, p_limit),
      points=z_null_lp1,
      limit=200000,
      epsabs=1E-55,
      epsrel=1E-10,
      full_output=2)
    res += integral[0]
    a_result[-1] -= 2 * res / constant_factor
  if z_null_lp2 >1:
    res = 0
    integral1 = integrate.quad(
        integrand_matrix_p,
        1,
        z_null_lp2-epsilon,
        args=(z_null_lp2,el_nr, nu_in, 2, p_limit),
        points=z_null_lp2,
        limit=200000,
        epsabs=1E-55,
        epsrel=1E-10,
        full_output=2)
    integral2 = integrate.quad(
        integrand_matrix_p,
        z_null_lp2+epsilon,
        100,
        args=(z_null_lp2,el_nr, nu_in, 2, p_limit),
        points=z_null_lp2,
        limit=200000,
        epsabs=1E-55,
        epsrel=1E-10,
        full_output=2)
    res += integral1[0]+integral2[0]
    a_result[-1] -= 4 * (res) / constant_factor
  else:
    res = 0
    integral = integrate.quad(
      integrand_matrix_p,
      1,
      100,
      args=(z_null_lp2,el_nr, nu_in, 2, p_limit),
      points=z_null_lp2,
      limit=200000,
      epsabs=1E-55,
      epsrel=1E-10,
      full_output=2)
    res += integral[0]
    a_result[-1] -= 4 * res / constant_factor
  fac = 2*math.pi
  c_result.append(f_a(el_nr,l,0,z_null,nu_in,n_0,p_limit).real * fac/constant_factor)
  c_result[-1] += (f_a(el_nr,2,1,z_null_ls,nu_in,n_0,p_limit).real * fac/constant_factor)
  c_result[-1] += 2*((f_c_0(el_nr,0,0,z_null_lp1,nu_in,n_0,p_limit) - 20 * f_c_2(el_nr,2,0,z_null_lp1,nu_in,n_0,p_limit))/3 * fac/constant_factor)
  c_result[-1] += 4*((f_c_0(el_nr,0,0,z_null_lp1,nu_in,n_0,p_limit) - 20 * f_c_2(el_nr,2,0,z_null_lp1,nu_in,n_0,p_limit))/3 * fac/constant_factor)

  lam = speed_of_light / nu_in * 1E10
  fpfdp = br.at_angstrom(lam,'Te')
  d_result.append(fpfdp[0])
  e_result.append(fpfdp[1])

def test_angular(nu,el_nr, theta0, array):
  z = h*nu / get_ionization_energy_2p1_2(el_nr)
  temp1 = f_b(
    Z=el_nr,
    l=2,
    g_k=0,
    z=z,
    nu_in = nu,
    n_0 = 2,
    p_limit = 3) * alpha_coef(2,1,0,theta0,0) * math.sin(theta0)
  temp2 = f_b(
    Z=el_nr,
    l=2,
    g_k=1,
    z=z,
    nu_in = nu,
    n_0 = 2,
    p_limit = 3) * alpha_coef(2,1,1,theta0,0) * math.cos(theta0)
  temp3 = f_b(
    Z=el_nr,
    l=2,
    g_k=2,
    z=z,
    nu_in = nu,
    n_0 = 2,
    p_limit = 3) * alpha_coef(2,1,2,theta0,0) * math.sin(theta0)
  temp4 = f_b(
    Z=el_nr,
    l=2,
    g_k=1,
    z=z,
    nu_in = nu,
    n_0 = 2,
    p_limit = 3) * beta_coef(2,1,1,theta0,0.5*math.pi) * math.cos(theta0)
  temp5 = f_b(
    Z=el_nr,
    l=2,
    g_k=2,
    z=z,
    nu_in = nu,
    n_0 = 2,
    p_limit = 3) * beta_coef(2,1,2,theta0,0.5*math.pi) * math.sin(theta0)
  array.append(math.sqrt(pow(temp1+temp2+temp3,2)+pow(temp4+temp5,2)))

def test_values():
  nu_in = get_ionization_energy_2s(el_nr) * h
  for l in range(0,10):
    b_res = 0
    c0_res = 0
    c2_res = 0
    d_res = 0
    for k in range(0,min(l,2)):
      b_res  += f_b(el_nr,l,k,1.001,nu_in,2,5)
      c0_res += f_c_0(el_nr,l,k,1.001,nu_in,2,5)
      c2_res += f_c_2(el_nr,l,k,1.001,nu_in,2,5)
      d_res  += f_d(el_nr,l,k,1.001,nu_in,2,5)
    print("------------------- l=%d ---------------------"%l)
    print("b: %f"%b_res)
    print("c0: %f"%c0_res)
    print("c2: %f"%c2_res)
    print("d: %f"%d_res)
    null = 0

def test_sum(theta0):
  a = -alpha_coef(0,0,0,theta0,0) * math.cos(theta0)
  b = -(alpha_coef(2,1,0,theta0,0)*math.sin(theta0) \
    - 3*alpha_coef(2,1,1,theta0,0) * math.cos(theta0)\
    - 6*alpha_coef(2,1,2,theta0,0) * math.sin(theta0) 
    )
  c0 = -0.5*( \
    alpha_coef(2,0,0,theta0,0) * math.cos(theta0)\
  + 3*alpha_coef(2,0,1,theta0,0) * math.sin(theta0)\
  + 6*beta_coef(2,0,2,theta0,math.pi/2)\
  - 6*alpha_coef(2,0,2,theta0,0) * math.cos(theta0)
  )
  c2 = 1/2 * 0.5 * ( \
    alpha_coef(2,2,0,theta0,0) * math.cos(theta0)\
  + 3*alpha_coef(2,2,1,theta0,0) * math.sin(theta0)\
  + 6*beta_coef(2,2,2,theta0,math.pi/2)\
  - 6*alpha_coef(2,2,2,theta0,0) * math.cos(theta0)
  )
  d = 1/2 * ( \
    alpha_bar_coef(2,2,0,theta0,math.pi/2) * math.cos(theta0)\
  + 2*alpha_bar_coef(2,2,1,theta0,math.pi/2) * math.sin(theta0)\
  + 6*beta_bar_coef(2,2,2,theta0,0)\
  - 6*alpha_bar_coef(2,2,2,theta0,math.pi/2) * math.cos(theta0)
  )
  return (a + b + c0 +c2 + d)/3


def cos_sin_func(phi, theta0, t, p, al):
  cosa = np.cos(theta0)
  sina = np.sin(theta0)
  return np.sin(t)\
    * pow(np.cos(t) * (cosa - 1) \
    + np.sin(t)\
    * sina * np.cos(phi-al),p)

def integrate_sin_funct(*args):
  #n = theta0
  #m = theta
  return np.vectorize(lambda n,m: integrate.quad(func=cos_sin_func,
      a=0,
      b=2*math.pi,
      args=(n,m)+args,
      epsabs = 1E-8,
      epsrel = 1E-5, 
      limit = 2000)[0])

def A(n,m):
  return (np.cos(n)-1)*np.cos(m)

def B(n,m):
  return np.sin(n)*np.sin(m)

def phi_integrated(m,n,p):
  #m = theta
  #n = theta0
  return (2*math.pi*np.sin(m)
    *(sum(pow(A(n,m),p-x) 
         * scipy.special.binom(p,x)
         * factorial2(x-1)/(factorial2(x)) 
         * pow(B(n,m),x) 
         * delta(x%2,0)
          for x in range(0,p+1)
         )
    )
  )

def integrate_theta(p):
  #n = theta0
  return np.vectorize(lambda n: integrate.quad(func=phi_integrated,
      a=0,
      b=math.pi,
      args=(n,p),
      epsabs = 1E-8,
      epsrel = 1E-5, 
      limit = 2000)[0])

def model(p):
  #n = theta0
  #m = theta
  return np.vectorize(lambda n,m: phi_integrated(m,n,p))

def theta0_function(n,p):
  if p%2 == 0:
    return (pow(
      (pow(np.cos(n)-1,2)
       +pow(np.sin(n),2))
      ,p/2) 
      *4*math.pi / (p+1)
    )
  else: return 0

def model2(p):
  return np.vectorize(lambda n: theta0_function(n,p))

threed_plot = False
test_1s = False
test_beta = False
cmaps = ['blue','orange','green','red','black','blue','orange','green','red','black']
if threed_plot == True:
  x = np.linspace(0,math.pi,30)
  y = np.linspace(0,math.pi,30)
  X,Y = np.meshgrid(x,y)
  fig = plt.figure()
  axes = fig.add_subplot(221,projection='3d')
  axes2 = fig.add_subplot(222,projection='3d')
  axes3 = fig.add_subplot(223,projection='3d')
  for p in range(0,3):
    Z = integrate_sin_funct(p,0)(X,Y)
    Z2 = model(p)(X,Y)
    axes.plot_wireframe(X,Y,Z, color = cmaps[p], label="p = %d"%p)
    axes.scatter3D(X,Y,Z2, cmap = cmaps[p], label="model p=%d"%p)
  axes.legend()
  axes.set_xlabel('theta_0')
  axes.set_ylabel('theta')
  axes.set_zlabel('Integral over theta')
  for p in range(3,6):
    Z = integrate_sin_funct(p,0)(X,Y)
    Z2 = model(p)(X,Y)
    axes2.plot_wireframe(X,Y,Z, color = cmaps[p], label="p = %d"%p)
    axes2.scatter3D(X,Y,Z2, cmap = cmaps[p],label="model p=%d"%p)
  axes2.legend()
  axes2.set_xlabel('theta_0')
  axes2.set_ylabel('theta')
  axes2.set_zlabel('Integral over theta')
  for p in range(6,9):
    Z = integrate_sin_funct(p,0)(X,Y)
    Z2 = model(p)(X,Y)
    axes3.plot_wireframe(X,Y,Z, color = cmaps[p], label="p = %d"%p)
    axes3.scatter3D(X,Y,Z2, cmap = cmaps[p],label="model p=%d"%p)
  axes3.legend()
  axes3.set_xlabel('theta_0')
  axes3.set_ylabel('theta')
  axes3.set_zlabel('Integral over theta')
  mng = plt.get_current_fig_manager()
  mng.window.state('zoomed')
  plt.show()

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
      np.sum(np.fromiter( (pow(-1,x) * pow(pi2_b * n * 1E10, 2*x) * (x+1)  \
        for x in range(p,-1,-1)
      ), dtype=np.float64))
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
    
  x = np.linspace(0,1.5,150)
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
  #axes.plot(thakkar_axis,thakkars[0],"--",label="Thakkar")
  #Z3 = 16/pow(4+pow(2*math.pi*np.sin(x/2)/lam,2),2)
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
  mng = plt.get_current_fig_manager()
  mng.window.state('zoomed')
  plt.show()
  exit()

if test_beta == True:
  d_result = []
  x = np.linspace(0,2*math.pi,200)
  for t in x:
    #a_result.append(-test_sum(t))
    #b_result.append(1/2+pow(math.cos(t),2)/2)
    #c_result.append(-math.cos(t))
    d_result.append(beta_coef(2,1,2,t,math.pi/2))
    e_result.append(-2*pow(math.sin(t),1)*pow(math.cos(t),0))
  
  fig, axes = plt.subplots(1,1)
    
  #axes.scatter(x,a_result,s=10,facecolors='none',edgecolors='b',marker='^',label="sum")
  #axes.plot(x,b_result,color='g',label="(1+cos2) / 2")
  #axes.plot(x,c_result,color='r',label="cos")
  axes.scatter(x,d_result,s=10,facecolors='none',edgecolors='r',marker='^',label="alpha")
  axes.plot(x,e_result,color='black',label="fit")
  axes.legend()
    
  plt.show()

t0 = 0
alp = 0
el_nr = 52
x = None

integration_test = False
EM_Hoenl_test = False
higher_p_test = True

if integration_test == True:
  start_nu_in = 0.8*get_ionization_energy_2s(el_nr)/h
  end_nu_in = 1.2*get_ionization_energy_1s(el_nr)/h
  x = np.linspace(start_nu_in,end_nu_in,150)
  for nu in x:
    test_integral_hönl(nu,el_nr, 1)
  plot_integral_test(get_ionization_energy_1s(el_nr)/h)

if EM_Hoenl_test == True:
  nu_in = 1.05*get_ionization_energy_1s(el_nr)/h
  x = np.linspace(1.0001,5.0001,150)
  for z in x:
    test_for_EM_and_Hoenl(el_nr,z,nu_in)
  plot_EM_Hoenl_test()

if higher_p_test == True:
  nu_in = 1.05*get_ionization_energy_1s(el_nr)/h
  x = np.linspace(1.0001,5.0001,150)
  l_max = 5
  k_ = 0
  x2_ = xn(nu_in,1,el_nr,0,2)
  x2_2 = xn(nu_in,2,el_nr,0, 2)
  t0 = alp = 0
  ct0 = ca = np.cos(t0)
  st0 = sa = np.sin(t0)

  if True:
    a_result = [0.0 for r in range(x.size)]
    b_result = [0.0 for r in range(x.size)]
    c_result = [0.0 for r in range(x.size)]
    b0_result = [0.0 for r in range(x.size)]
    b1_result = [0.0 for r in range(x.size)]
    b2_result = [0.0 for r in range(x.size)]
    d_result = [0.0 for r in range(x.size)]
    e_result = [0.0 for r in range(x.size)]
    e0_result = [0.0 for r in range(x.size)]
    e1_result = [0.0 for r in range(x.size)]
    e2_result = [0.0 for r in range(x.size)]
    f_result = [0.0 for r in range(x.size)]
  
    p0_result = [0.0 for r in range(x.size)]
    p2_0_result = [0.0 for r in range(x.size)]
    p2_1_result = [0.0 for r in range(x.size)]
    p4_0_result = [0.0 for r in range(x.size)]
    p4_1_result = [0.0 for r in range(x.size)]
    p4_2_result = [0.0 for r in range(x.size)]
  
    g_result = [0.0 for r in range(x.size)]
    h_result = [0.0 for r in range(x.size)]
    i_result = [0.0 for r in range(x.size)]
    j_result = [0.0 for r in range(x.size)]
    k_result = [0.0 for r in range(x.size)]
    l_result = [0.0 for r in range(x.size)]
  
    m_result = [0.0 for r in range(x.size)]
    n_result = [0.0 for r in range(x.size)]
    o_result = [0.0 for r in range(x.size)]
  for i,z in enumerate(x):
    for l in range(l_max+1):
      res_0 = f_a_for_p(el_nr,l,k_,z,nu_in,1,0)[0]
      res_2 = f_a_for_p(el_nr,l,k_,z,nu_in,1,2)
      res_4 = f_a_for_p(el_nr,l,k_,z,nu_in,1,4)
      a_result[i] += (res_0)
      b_result[i] += ((res_2[0]+res_2[2]))
      c_result[i] += (res_2[1])
      b0_result[i] += ((res_4[0]+res_4[-1]))
      b1_result[i] += ((res_4[1]+res_4[-2]))
      b2_result[i] += ((res_4[2]))
      res_0 = f_a_for_p(el_nr,l,k_,z,nu_in,2,0)[0]
      res_2 = f_a_for_p(el_nr,l,k_,z,nu_in,2,2)
      res_4 = f_a_for_p(el_nr,l,k_,z,nu_in,2,4)
      d_result[i] += (res_0)
      e_result[i] += ((res_2[0]+res_2[2]))
      f_result[i] += (res_2[1])
      e0_result[i] += ((res_4[0]+res_4[-1]))
      e1_result[i] += ((res_4[1]+res_4[-2]))
      e2_result[i] += ((res_4[2]))

      for _k in range(3):
        b_l   =   f_b_for_p(el_nr,l,_k,z,nu_in,2,0)[0]
        c_0_l = f_c_0_for_p(el_nr,l,_k,z,nu_in,2,0)[0] 
        c_2_l = f_c_2_for_p(el_nr,l,_k,z,nu_in,2,0)[0]
        d_l   =   f_d_for_p(el_nr,l,_k,z,nu_in,2,0)[0]
        if _k == 0:
          b_l   *= -st0 * alpha_coef(l,1,_k,t0,alp)
          c_0_l *= ct0 * alpha_coef(l,0,_k,t0,alp) * ca
          c_2_l *= ct0 * alpha_coef(l,2,_k,t0,alp) * ca
          d_l   *= ct0 * alpha_bar_coef(l,2,_k,t0,alp) * sa
        elif _k == 1:
          b_l   *= ct0 * alpha_coef(l,1,_k,t0,alp)
          c_0_l *= st0 * alpha_coef(l,0,_k,t0,alp) * ca
          c_2_l *= st0 * alpha_coef(l,2,_k,t0,alp) * ca
          d_l   *= st0 * alpha_bar_coef(l,2,_k,t0,alp) * sa
        elif _k == 2:
          b_l   *= -st0 * alpha_coef(l,1,_k,t0,alp)
          c_0_l *= ct0 * alpha_coef(l,0,_k,t0,alp) * ca - sa * beta_coef(l,0,_k,t0,alp)
          c_2_l *= ct0 * alpha_coef(l,2,_k,t0,alp) * ca - sa * beta_coef(l,2,_k,t0,alp)
          d_l   *= ct0 * alpha_bar_coef(l,2,_k,t0,alp) * sa + ca * beta_bar_coef(l,2,_k,t0,alp)
        res_0 = 1.0/3.0 * (b_l + c_0_l + c_2_l + d_l)
        p0_result[i] += res_0
        #res_2 = f_a_for_p(el_nr,l,k,z,nu_in,2,2)
        #res_4 = f_a_for_p(el_nr,l,k,z,nu_in,2,4)

    g_result[i] += (f_s_1_hoenl(z)*constant_factor)
    h_result[i] += (f_s_2_1_hoenl(z)*constant_factor*x2_)
    i_result[i] += (f_s_2_2_hoenl(z)*constant_factor*x2_)

    j_result[i] += (f_s_1_EM(z)*constant_factor)
    k_result[i] += (f_s_2_1_EM(z)*constant_factor*x2_2)
    l_result[i] += (f_s_2_2_EM(z)*constant_factor*x2_2)

    m_result[i] += (f_p_1_EM(z)*constant_factor)
    n_result[i] += (f_p_2_1_EM(z)*constant_factor*x2_2)
    o_result[i] += (f_p_2_2_EM(z)*constant_factor*x2_2)
    p2_0_result[i] = m_result[i] / p0_result[i]
  
  fig, axes = plt.subplots(2,2)
  axes[0,0].scatter(x,g_result,s=20,facecolors='none',edgecolors='b',marker='^',label="Hönl fs_1^(0)")
  axes[0,0].scatter(x,h_result,s=20,facecolors='none',edgecolors='g',marker='^',label="Hönl fs_1^(2)")
  axes[0,0].scatter(x,i_result,s=20,facecolors='none',edgecolors='r',marker='^',label="Hönl fs_2^(2)")
  axes[1,0].scatter(x,j_result,s=20,facecolors='none',edgecolors='b',marker='^',label="EM fs_1^(0)")
  axes[1,0].scatter(x,k_result,s=20,facecolors='none',edgecolors='g',marker='^',label="EM fs_1^(2)")
  axes[1,0].scatter(x,l_result,s=20,facecolors='none',edgecolors='r',marker='^',label="EM fs_2^(2)")
  axes[0,0].plot(x,a_result,color='b',label="p=0, 0/0")
  axes[0,0].plot(x,b_result,color='g',label="p=2, 0/2 + 2/0")
  axes[0,0].plot(x,c_result,color='r',label="p=2, 1/1")
  axes[0,0].plot(x,b0_result,label="p=4, 0/4 + 4/0")
  axes[0,0].plot(x,b1_result,label="p=4, 1/3 + 3/1")
  axes[0,0].plot(x,b2_result,label="p=4, 2/2")
  axes[1,0].plot(x,d_result,color='b',label="p=0, 0/0")
  axes[1,0].plot(x,e_result,color='g',label="p=2, 0/2 + 2/0")
  axes[1,0].plot(x,f_result,color='r',label="p=2, 1/1")
  axes[1,0].plot(x,e0_result,label="p=4, 0/4 + 4/0")
  axes[1,0].plot(x,e1_result,label="p=4, 1/3 + 3/1")
  axes[1,0].plot(x,e2_result,label="p=4, 2/2")
  axes[1,0].legend()
  axes[1,0].axhline(0,1,5,ls='--')
  axes[0,0].legend()
  axes[0,0].axhline(0,1,5,ls='--')

  axes[0,1].scatter(x,m_result,s=20,facecolors='none',edgecolors='b',marker='^',label="EM fp_1^(0)")
  axes[0,1].scatter(x,n_result,s=20,facecolors='none',edgecolors='g',marker='^',label="EM fp_1^(2)")
  axes[0,1].scatter(x,o_result,s=20,facecolors='none',edgecolors='r',marker='^',label="EM fp_2^(2)")
  axes[0,1].plot(x,p0_result,color='b',label="p=0, 0/0")
  axes[0,1].plot(x,p2_0_result,color='g',label="p=2, 0/2 + 2/0")
  #axes[0,1].plot(x,p2_1_result,color='r',label="p=2, 1/1")
  #axes[0,1].plot(x,p4_0_result,label="p=4, 0/4 + 4/0")
  #axes[0,1].plot(x,p4_1_result,label="p=4, 1/3 + 3/1")
  #axes[0,1].plot(x,p4_2_result,label="p=4, 2/2")
  axes[0,1].legend()
  mng = plt.get_current_fig_manager()
  mng.window.state('zoomed')
  plt.show()

