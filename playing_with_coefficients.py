from legendre_plynomials import *
import matplotlib.pyplot as plt
#runny = np.linspace(0,math.pi,100)
#plt.plot(runny, 2*pow(np.cos(0.5*runny),2)-1,'r--', runny, np.cos(runny), 'bs')
#plt.show()
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