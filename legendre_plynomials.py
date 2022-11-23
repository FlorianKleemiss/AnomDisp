import scipy
import scipy.special
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def left_side(l,m,x, theta0,alpha=math.pi/2,Phi=0):
    result = []
    for x_ in x:
        cos_t = np.cos(theta0) * np.cos(x_) - np.sin(theta0) * np.sin(x_) * np.cos(Phi)
        sin_t = np.sqrt(1-cos_t*cos_t)
        eiphi = np.exp(complex(0,1)*alpha) / sin_t * complex(np.sin(theta0)*np.cos(x_)+np.cos(theta0)*np.sin(x_)*np.cos(Phi),np.sin(x_)*np.sin(Phi))
        res = scipy.special.lpmv(m,l,cos_t) * pow(eiphi,m)
        if m%2 == 0:
          result.append(res)
        else:
          result.append(-res)
    return result
  
def right_side11(x,theta0):
    ct0 = np.cos(theta0)
    st0 = np.sin(theta0)
    factors = [(0.5-0.5*ct0),st0,-(0.5+0.5*ct0)]
    res = []
    for i in x:
        arg = np.cos(i)
        res_ = 0
        for m_ in range(-1,2):
            res_ += factors[m_+1] * scipy.special.lpmv(m_,1,arg)
        res.append(res_)
    return res

def right_side21(x,theta0):
    ct0h = np.cos(theta0/2)
    st0h = np.sin(theta0/2)
    ct0h2 = ct0h * ct0h
    st0h2 = st0h * st0h
    ct0h3 = ct0h2 * ct0h
    st0h3 = st0h2 * st0h
    factors = [ct0h*st0h3,
               (3*ct0h2*st0h2-st0h3*st0h),
               6*(ct0h3*st0h-ct0h*st0h3),
               -ct0h3*ct0h-3*ct0h2*st0h2,
               -ct0h3*st0h]
    res = []
    for i in x:
        arg = np.cos(i)
        res_ = 0
        for m_ in range(-2,3):
            res_ += factors[m_+2] * scipy.special.lpmv(m_,2,arg)
        res.append(res_)
    return res

def right_side_abs(l,m,x,theta0,alpha=math.pi/2):
    factors2 = []
    for m_ in range(-l,l+1):
        temp = C_coef(l,m,m_,theta0,math.pi/2)
        #To account for the sign change in even legendre polynomials
        if (m%2 == 0):
            if (m_%2) != 0:
                temp = -temp.real
            else:
                temp = temp.real
        else:
            if (m_%2) != 0:
                temp = -temp.imag
            else:
                temp = temp.imag
        factors2.append(temp)
    res = []
    for i in x:
        arg = np.cos(i)
        res_ = 0
        for m_ in range(-l,l+1):
            factor_ = factors2[m_+l]
            res_ += factor_ * scipy.special.lpmv(m_,l,arg)
        res.append(res_)
    return res
  
def right_side_complex(l,m,Theta,theta0,alpha=math.pi/2,Phi=0):
    factors2 = []
    for m_ in range(-l,l+1):
        temp = C_coef(l,m,m_,theta0,alpha)
        #To account for the sign change in even legendre polynomials
        if (m_%2) != 0:
            temp = -temp
        factors2.append(temp)
    res = []
    for i in Theta:
        arg = np.cos(i)
        res_ = 0
        for m_ in range(-l,l+1):
            factor_ = factors2[m_+l]
            res_ += factor_ * scipy.special.lpmv(m_,l,arg) * np.exp(complex(0,1)*m_*Phi)
        res.append(res_)
    return res

def C_coef(l,m,m_,theta0,alpha=math.pi/2):
    factor1 = np.exp(complex(0,1)*m*alpha) 
    factor1 = round(factor1.real, 10) + round(factor1.imag, 10) * 1j
    factor2 = scipy.special.factorial(l-m_)/scipy.special.factorial(l-m)
    factor3 = 0
    rhos = []
    for rho in range(max(0,m_-m),min(l-m,l+m_)+1):
        temp = pow(-1,rho) * scipy.special.binom(l-m,rho) * scipy.special.binom(l+m,l+m_-rho)
        temp *= pow(np.cos(theta0/2),2*l-m+m_-2*rho)
        temp *= pow(np.sin(theta0/2),m-m_+2*rho)
        factor3 += temp
        rhos.append(rho)
    return factor1 * factor2 * factor3 #/ correction_negative

def alpha_coef(l,m,m_,theta0, alpha):
  st0 = np.sin(theta0/2)
  ct0 = np.cos(theta0/2)
  if m > l or m_ > l:
    return 0
  #special case for m_ = 0
  if m_ == 0:
    part1 = np.cos(m*alpha) * scipy.special.factorial(l) / scipy.special.factorial(l-m) \
      * pow(ct0,2*l)
    sum = 0
    for rho in range(l-m+1):
      sum += pow(-1,rho) * scipy.special.binom(l-m,rho) * scipy.special.binom(l+m,l-rho) \
        * pow(st0/ct0,m+2*rho)
    return part1 * sum
  else:
    part1 = np.cos(m*alpha) * scipy.special.factorial(l-m_) / scipy.special.factorial(l-m)
    if part1 == 0:
      return 0.0
    sum = 0
    #print("lower limit:" +str(max(0,-m-m_,-m+m_)) + " upper limit: " + str(min(l-m, l-m_, l+m_)+1))
    for rho in range(0, l-m+1):
      f0 = scipy.special.binom(l-m,rho)
      f1 = scipy.special.binom(l+m,m+m_+rho)
      f2 = scipy.special.binom(l+m,m-m_+rho)
      if f1 != 0 and f2 != 0:
        q1 = pow(st0,m+m_+2*rho) * pow(ct0,2*l-m-2*rho-m_)
        q2 = pow(st0,m-m_+2*rho) * pow(ct0,2*l-m-2*rho+m_)
        if q1 != 0 and q2 != 0:
          s = (pow(-1,m_) * f1 * q1 + f2 * q2)
        elif q1 != 0:
          s = (pow(-1,m_) * f1 * q1)
        elif q2 != 0:
          s = f2*q2
        else:
          s = 0
      elif f1 == 0:
        q2 = pow(st0,m-m_+2*rho) * pow(ct0,2*l-m-2*rho+m_)
        s = f2 * q2
      elif f2 == 0:
        q1 = pow(st0,m+m_+2*rho) * pow(ct0,2*l-m-2*rho-m_)
        s = pow(-1,m_) * f1 * q1
      else:
        s = 0
      sum += pow(-1,rho) * f0 * s
  return part1 * sum

def alpha_bar_coef(l,m,m_,theta0, alpha):
  if l-m < 0:
    return 0
  st0 = np.sin(theta0/2)
  ct0 = np.cos(theta0/2)
  #special case for m_ = 0
  if m_ == 0:
    part1 = np.sin(m*alpha) * scipy.special.factorial(l) / scipy.special.factorial(l-m) \
      * pow(ct0,2*l)
    sum = 0
    for rho in range(l-m+1):
      sum += pow(-1,rho) * scipy.special.binom(l-m,rho) * scipy.special.binom(l+m,l-rho) \
        * pow(st0/ct0,m+2*rho)
    return part1 * sum
  else:
    part1 = np.sin(m*alpha) * scipy.special.factorial(l-m_) / scipy.special.factorial(l-m)
    if part1 == 0:
      return 0.0
    sum = 0
    #print("lower limit:" +str(max(0,-m-m_,-m+m_)) + " upper limit: " + str(min(l-m, l-m_, l+m_)+1))
    for rho in range(0, l-m+1):
      f0 = scipy.special.binom(l-m,rho)
      f1 = scipy.special.binom(l+m,m+m_+rho)
      f2 = scipy.special.binom(l+m,m-m_+rho)
      if f1 != 0 and f2 != 0:
        q1 = pow(st0,m+m_+2*rho) * pow(ct0,2*l-m-2*rho-m_)
        q2 = pow(st0,m-m_+2*rho) * pow(ct0,2*l-m-2*rho+m_)
        if q1 != 0 and q2 != 0:
          s = (pow(-1,m_) * f1 * q1 + f2 * q2)
        elif q1 != 0:
          s = (pow(-1,m_) * f1 * q1)
        elif q2 != 0:
          s = f2*q2
        else:
          s = 0
      elif f1 == 0:
        q2 = pow(st0,m-m_+2*rho) * pow(ct0,2*l-m-2*rho+m_)
        s = f2 * q2
      elif f2 == 0:
        q1 = pow(st0,m+m_+2*rho) * pow(ct0,2*l-m-2*rho-m_)
        s = pow(-1,m_) * f1 * q1
      else:
        s = 0
      sum += pow(-1,rho) * f0 * s
  return part1 * sum

def beta_coef(l,m,m_,theta0, alpha):
  if l-m < 0:
    return 0
  if m_ == 0:
    return 0.0
  if m == 0:
    return 0.0
  st0 = np.sin(theta0/2)
  ct0 = np.cos(theta0/2)
  part1 = np.sin(m*alpha) * scipy.special.factorial(l-m_) / scipy.special.factorial(l-m)
  if part1 == 0:
    return 0.0
  sum = 0
  #print("lower limit:" +str(max(0,-m-m_,-m+m_)) + " upper limit: " + str(min(l-m, l-m_, l+m_)+1))
  for rho in range(0, l-m+1):
    f0 = scipy.special.binom(l-m,rho)
    f1 = scipy.special.binom(l+m,m+m_+rho)
    f2 = scipy.special.binom(l+m,m-m_+rho)
    if f1 != 0 and f2 != 0:
      q1 = pow(st0,m+m_+2*rho) * pow(ct0,2*l-m-2*rho-m_)
      q2 = pow(st0,m-m_+2*rho) * pow(ct0,2*l-m-2*rho+m_)
      if q1 != 0 and q2 != 0:
        s = (pow(-1,m_) * f1 * q1 - f2 * q2)
      elif q1 != 0:
        s = (pow(-1,m_) * f1 * q1)
      elif q2 != 0:
        s = -f2*q2
      else:
        s = 0
    elif f1 == 0:
      q2 = pow(st0,m-m_+2*rho) * pow(ct0,2*l-m-2*rho+m_)
      s = -f2 * q2
    elif f2 == 0:
      q1 = pow(st0,m+m_+2*rho) * pow(ct0,2*l-m-2*rho-m_)
      s = pow(-1,m_) * f1 * q1
    else:
      s = 0
    sum += pow(-1,rho) * f0 * s
  return part1 * sum

def beta_bar_coef(l,m,m_,theta0, alpha):
  if l-m < 0:
    return 0
  if m_ == 0:
    return 0.0
  st0 = np.sin(theta0/2)
  ct0 = np.cos(theta0/2)
  part1 = np.cos(m*alpha) * scipy.special.factorial(l-m_) / scipy.special.factorial(l-m)
  if part1 == 0:
    return 0.0
  sum = 0
  #print("lower limit:" +str(max(0,-m-m_,-m+m_)) + " upper limit: " + str(min(l-m, l-m_, l+m_)+1))
  for rho in range(0, l-m+1):
    f0 = scipy.special.binom(l-m,rho)
    f1 = scipy.special.binom(l+m,m+m_+rho)
    f2 = scipy.special.binom(l+m,m-m_+rho)
    if f1 != 0 and f2 != 0:
      q1 = pow(st0,m+m_+2*rho) * pow(ct0,2*l-m-2*rho-m_)
      q2 = pow(st0,m-m_+2*rho) * pow(ct0,2*l-m-2*rho+m_)
      if q1 != 0 and q2 != 0:
        s = (pow(-1,m_+1) * f1 * q1 + f2 * q2)
      elif q1 != 0:
        s = (pow(-1,m_+1) * f1 * q1)
      elif q2 != 0:
        s = f2*q2
      else:
        s = 0
    elif f1 == 0:
      q2 = pow(st0,m-m_+2*rho) * pow(ct0,2*l-m-2*rho+m_)
      s = f2 * q2
    elif f2 == 0:
      q1 = pow(st0,m+m_+2*rho) * pow(ct0,2*l-m-2*rho-m_)
      s = pow(-1,m_+1) * f1 * q1
    else:
      s = 0
    sum += pow(-1,rho) * f0 * s
  return part1 * sum
  
def right_side_44_real(l,m,Theta,theta0,alpha=math.pi/2,Phi=0):
    factors2 = []
    for m_ in range(0,l+1):
        temp = float(alpha_coef(l,m,m_,theta0,alpha))
        temp2 = float(beta_coef(l,m,m_,theta0,alpha))
        #To account for the sign change in even legendre polynomials
        if (m_%2) != 0:
            temp = -temp
            temp2 = -temp2
        factors2.append(temp * np.cos(m_*Phi) + temp2 * np.sin(m_*Phi))
    res = []
    for i in Theta:
        arg = np.cos(i)
        res_ = 0
        for m_ in range(0,l+1):
            factor_ = factors2[m_]
            res_ += factor_ * scipy.special.lpmv(m_,l,arg)
        res.append(res_)
    return res
  
def right_side_44_imag(l,m,Theta,theta0,alpha=math.pi/2,Phi=0):
    factors2 = []
    for m_ in range(0,l+1):
        temp = float(alpha_bar_coef(l,m,m_,theta0,alpha))
        temp2 = float(beta_bar_coef(l,m,m_,theta0,alpha))
        #To account for the sign change in even legendre polynomials
        if (m_%2) != 0:
            temp = -temp
            temp2 = -temp2
        factors2.append(temp * np.cos(m_*Phi) + temp2 * np.sin(m_*Phi))
    res = []
    for i in Theta:
        arg = np.cos(i)
        res_ = 0
        for m_ in range(0,l+1):
            factor_ = factors2[m_]
            res_ += factor_ * scipy.special.lpmv(m_,l,arg)
        res.append(res_)
    return res  

def test():
  import random
  random.seed()
  theta0 = random.random()*math.pi
  alpha = random.random()*2*math.pi
  Phi = random.random()*2*math.pi
  #print(theta0)
  #print(alpha)
  #print(Phi)
  
  l=5
  m=3
  Theta = np.linspace(0,math.pi,100)
  
  result_43_right_side = right_side_complex(l,m,Theta,theta0,alpha,Phi)
  result_43_right_sider=[]
  result_43_right_sidei=[]
  for i in result_43_right_side:
    result_43_right_sider.append(i.real)
    result_43_right_sidei.append(i.imag)
  result_44_right_side_real = right_side_44_real(l,m,Theta,theta0,alpha,Phi)
  result_44_right_side_imag = right_side_44_imag(l,m,Theta,theta0,alpha,Phi)
  
  
  plt.plot(Theta,result_43_right_sider)
  plt.plot(Theta,result_43_right_sidei)
  plt.scatter(Theta,result_44_right_side_real)
  plt.scatter(Theta,result_44_right_side_imag)
  plt.show()
  
  
  print("We tested 43 and 44")
  
  x = np.linspace(0,math.pi,100)
  y_1 = left_side(1,1,x,theta0,alpha,Phi)
  y_1r = []
  y_1c = []
  for i in y_1:
    y_1r.append(i.real)
    y_1c.append(i.imag)
  y_2 = right_side_complex(1,1,x,theta0,alpha,Phi)
  y_2r = []
  y_2c = []
  for i in y_2:
    y_2r.append(i.real)
    y_2c.append(i.imag)
    
  y_3 = left_side(2,1,x,theta0)
  y_3 = left_side(2,1,x,theta0,alpha,Phi)
  y_3r = []
  y_3c = []
  for i in y_3:
    y_3r.append(i.real)
    y_3c.append(i.imag)
  #y_4 = right_side21(x,theta0)
  y_4 = right_side_complex(2,1,x,theta0,alpha,Phi)
  y_4r = []
  y_4c = []
  for i in y_4:
    y_4r.append(i.real)
    y_4c.append(i.imag)
  
  y_5 = left_side(2,2,x,theta0,alpha,Phi)
  y_5r = []
  y_5c = []
  for i in y_5:
    y_5r.append(i.real)
    y_5c.append(i.imag)
  y_6 = right_side_complex(2,2,x,theta0,alpha,Phi)
  y_6r = []
  y_6c = []
  for i in y_6:
    y_6r.append(i.real)
    y_6c.append(i.imag)
    
  y_7 = left_side(3,2,x,theta0,alpha,Phi)
  y_7r = []
  y_7c = []
  for i in y_7:
    y_7r.append(i.real)
    y_7c.append(i.imag)
  y_8 = right_side_complex(3,2,x,theta0,alpha,Phi)
  y_8r = []
  y_8c = []
  for i in y_8:
    y_8r.append(i.real)
    y_8c.append(i.imag)
  
  y_9 = left_side(5,3,x,theta0,alpha,Phi)
  y_9r = []
  y_9c = []
  for i in y_9:
    y_9r.append(i.real)
    y_9c.append(i.imag)
  y_10 = right_side_complex(5,3,x,theta0,alpha,Phi)
  y_10r = []
  y_10c = []
  for i in y_10:
    y_10r.append(i.real)
    y_10c.append(i.imag)
    
  y_11 = left_side(5,-3,x,theta0,alpha,Phi)
  y_11r = []
  y_11c = []
  for i in y_11:
    y_11r.append(i.real)
    y_11c.append(i.imag)
  y_12 = right_side_complex(5,-3,x,theta0,alpha,Phi)
  y_12r = []
  y_12c = []
  for i in y_12:
    y_12r.append(i.real)
    y_12c.append(i.imag)
  
  y_m1 = left_side(1,-1,x,theta0,alpha,Phi)
  y_m1r = []
  y_m1c = []
  for i in y_m1:
    y_m1r.append(i.real)
    y_m1c.append(i.imag)
  #y_2 = right_side11(x,theta0)
  y_m2 = right_side_complex(1,-1,x,theta0,alpha,Phi)
  y_m2r = []
  y_m2c = []
  for i in y_m2:
    y_m2r.append(i.real)
    y_m2c.append(i.imag)
  
  fig, axs = plt.subplots(4,2,sharex=True)
  
  axs[2,0].plot(x,y_m1r,'-',label="left_side real")
  axs[2,0].plot(x,y_m1c,'-',label="left_side imag")
  #axs[0,0].plot(x,y_2,'-',label="right_side")
  axs[2,0].scatter(x,y_m2r,s=10,facecolors='none',edgecolors='b',label="right_side real")
  axs[2,0].scatter(x,y_m2c,s=10,facecolors='none',edgecolors='r',label="right_side imag")
  axs[2,0].legend(loc='upper right')
  axs[2,0].set_title('l=1; m=-1')
  axs[2,0].set(xlabel=r'\Theta')
  
  axs[0,0].plot(x,y_1r,'-',label="left_side real")
  axs[0,0].plot(x,y_1c,'-',label="left_side imag")
  #axs[0,0].plot(x,y_2,'-',label="right_side")
  axs[0,0].scatter(x,y_2r,s=10,facecolors='none',edgecolors='b',label="right_side real")
  axs[0,0].scatter(x,y_2c,s=10,facecolors='none',edgecolors='r',label="right_side imag")
  axs[0,0].legend(loc='upper right')
  axs[0,0].set_title('l=m=1')
  axs[0,0].set(xlabel=r'\Theta')
  
  axs[0,1].plot(x,y_3r,'-',label="left_side real")
  axs[0,1].plot(x,y_3c,'-',label="left_side imag")
  #axs[0,1].plot(x,y_4,'-',label="right_side")
  axs[0,1].scatter(x,y_4r,s=10,facecolors='none',edgecolors='b',label="right_side real")
  axs[0,1].scatter(x,y_4c,s=10,facecolors='none',edgecolors='r',label="right_side imag")
  axs[0,1].legend(loc='upper right')
  axs[0,1].set_title('l=2; m=1')
  axs[0,1].set(xlabel=r'\Theta')
  
  axs[1,0].plot(x,y_5r,'-',label="left_side real")
  axs[1,0].plot(x,y_5c,'-',label="left_side imag")
  axs[1,0].scatter(x,y_6r,s=10,facecolors='none',edgecolors='b',label="right_side real")
  axs[1,0].scatter(x,y_6c,s=10,facecolors='none',edgecolors='r',label="right_side imag")
  axs[1,0].legend(loc='upper right')
  axs[1,0].set_title('l=m=2')
  axs[1,0].set(xlabel=r'\Theta')
  
  axs[1,1].plot(x,y_7r,'-',label="left_side real")
  axs[1,1].plot(x,y_7c,'-',label="left_side imag")
  axs[1,1].scatter(x,y_8r,s=10,facecolors='none',edgecolors='b',label="right_side real")
  axs[1,1].scatter(x,y_8c,s=10,facecolors='none',edgecolors='r',label="right_side imag")
  axs[1,1].legend(loc='upper right')
  axs[1,1].set_title('l=3; m=2')
  axs[1,1].set(xlabel=r'\Theta')
  
  axs[2,1].plot(x,y_9r,'-',label="left_side real")
  axs[2,1].plot(x,y_9c,'-',label="left_side imag")
  axs[2,1].scatter(x,y_10r,s=10,facecolors='none',edgecolors='b',label="right_side real")
  axs[2,1].scatter(x,y_10c,s=10,facecolors='none',edgecolors='r',label="right_side imag")
  axs[2,1].legend(loc='upper right')
  axs[2,1].set_title('l=5; m=3')
  axs[2,1].set(xlabel=r'\Theta')
  
  axs[3,1].plot(x,y_11r,'-',label="left_side real")
  axs[3,1].plot(x,y_11c,'-',label="left_side imag")
  axs[3,1].scatter(x,y_12r,s=10,facecolors='none',edgecolors='b',label="right_side real")
  axs[3,1].scatter(x,y_12c,s=10,facecolors='none',edgecolors='r',label="right_side imag")
  axs[3,1].legend(loc='upper right')
  axs[3,1].set_title('l=5; m=-3')
  axs[3,1].set(xlabel=r'\Theta')
  plt.show()
  
  
  print("tada")

if __name__ == "__main__":
  test()