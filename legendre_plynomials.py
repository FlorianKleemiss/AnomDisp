import scipy.special
import numpy as np

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
        temp = m-m_+2*rho
        if temp < 0 and st0 == 0:
          q2 = 0
        else:
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
