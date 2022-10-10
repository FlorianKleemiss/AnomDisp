from legendre_plynomials import *
from matrix_coefficients_v2 import *
import matplotlib.pyplot as plt

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
  
def f_s_1(z):
  part1 = 512*(z+3)/(3*pow(z,4))
  z_1=z-1
  sqrt_var = numpy.sqrt(z_1)
  part2 = numpy.exp(-8/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-4*math.pi/sqrt_var)
  return part1*part2/part3

def f_s_1_hoenl(z):
  part1 = 64/(3*pow(z,3))
  z_1=z-1
  sqrt_var = numpy.sqrt(z_1)
  part2 = numpy.exp(-4/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-2*math.pi/sqrt_var)
  return part1*part2/part3
  
def f_s_2(z):
  part1 = 2048*pow(z-5,2)*(z+3)/(15*pow(z,7))
  z_1=z-1
  sqrt_var = numpy.sqrt(z_1)
  part2 = numpy.exp(-8/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-4*math.pi/sqrt_var)
  return part1*part2/part3

def f_s_2_hoenl(z):
  part1 = 256*(4*z-3)/(15*pow(z,5))
  z_1=z-1
  sqrt_var = numpy.sqrt(z_1)
  part2 = numpy.exp(-4/sqrt_var*numpy.arctan(sqrt_var))
  part3 = 1-numpy.exp(-2*math.pi/sqrt_var)
  return part1*part2/part3

t0 = math.pi/3
alp = math.pi/1.823
el_nr = 6
nu_in = 3800 * h
nu_2 = get_ionization_energy_2s(el_nr) * h
x = np.linspace(1.001,4.501,200)
a_result = []
b_result = []
c_result = []
d_result = []
#g_result = []
x2 = pow(2* q(nu_in)/b(2,0, el_nr),2)
constant_factor = 4*math.pi*math.pi*el_mass/h/h
for z in x:
  #a_,b_,c_,d_ = calc_S_values(t0, alp, 2,4,z,z2,nu_in,el_nr)
  #temp1 = f_a(el_nr,1,z,nu_in,2,1)
  temp_voll = f_a(
    Z=el_nr,
    l=1,
    z=z,
    nu_in = nu_in,
    n_0 = 1,
    p_limit = 3)

  #temp = temp1 - temp_p
  
  #temp /= x2
  #temp *= 4*math.pi*math.pi*el_mass/h/h
  #temp2 = temp

  temp1 = f_s_1_hoenl(z) * constant_factor
  a_result.append(temp_voll/temp1)
  temp2 = f_s_2_hoenl(z) * constant_factor *x2
  temp_voll_2 = f_a(
    Z=el_nr,
    l=2,
    z=z,
    nu_in = nu_in,
    n_0 = 1,
    p_limit = 5
  )
  b_result.append(temp_voll_2/temp2)
  #temp1 /= x2
  #temp1 *= constant_factor
  #temp *= constant_factor
  
  #b_result.append(temp)
  #c_result.append(temp1)
  null = 0
  #d_result.append(temp2)
  
  #g_result.append(temp)

fig = plt.figure()
axes = fig.add_subplot(1,1,1)
axes.scatter(x,a_result,s=10,facecolors='none',edgecolors='b',label="voll l=1")
axes.scatter(x,b_result,s=10,facecolors='none',edgecolors='g',label="voll l=2")
#axes.scatter(x,c_result,s=10,facecolors='none',edgecolors='r',label="f l=1")
#axes.scatter(x,d_result,s=10,facecolors='none',edgecolors='y',label="f l=2")
#axes.plot(x,g_result,'--',label="f_table1")
axes.legend()

plt.show()