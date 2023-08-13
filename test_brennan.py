from brennan import brennan
from constants_and_atomic_properties import *

#energy = float(input("Please enter Wavelength in Angstrom: "))
#z = input("Please Enter Element Name: ")
#energy = angstrom2eV / wavelength / 1000

from scipy import constants as conts

def ret_energy(wavelength):
  return (conts.h * conts.c) / (wavelength*1.60218e-19) * 1E10

def ret_wl(energy):
  return (conts.h * conts.c) / (energy*1.60218e-19) * 1E10


import matplotlib.pyplot as plt
Z = elements[83]
bren = brennan()
minimal = 200
maximal = 2000
stepsize = 1
for i in range(50,60):
  Z = elements[i]
  x = []
  y = []
  z = []
  r = []
  c = []
  for step in range(int(minimal),int(maximal),int(stepsize)):
    x.append(step*0.001)
    temp = bren.at_angstrom(wavelength=step*0.001,element=Z)
    #ray = bren.get_ray_at_angstrom_in_electrons(step*0.001,element=Z)
    #comp = bren.get_comp_at_angstrom_in_electrons(step*0.001,element=Z)
    y.append(temp[0])
    #z.append(temp[1])
    #r.append(ray)
    #c.append(comp)
  
  plt.plot(x,y,':',label="f' "+Z)
  #plt.plot(x,z,':',label="f\"")
  #plt.plot(x,c,'--',label="compton")
  #plt.plot(x,r,'--',label="rayleigh")
plt.xlim(maximal*0.001,minimal*0.001)
plt.legend()
plt.xlabel("Wavelength /A")
plt.ylabel("contribution")
plt.subplots_adjust(left=0.04, bottom=0.05, right=1.0, top=1.0, wspace=0.15, hspace=0.05)
mng = plt.get_current_fig_manager()
mng.window.state('zoomed')
plt.show()

x = []
m = []
mp = []
for step in range(int(minimal),int(maximal),int(stepsize)):
  x.append(step*0.001)
  temp = bren.get_mu_at_angstrom(wavelength=step*0.001,element=Z)
  m.append(temp)
  #temp = bren.get_mu_pure_at_angstrom(wavelength=step*0.001,element=Z)
  #mp.append(temp)
plt.xlim(maximal*0.001,minimal*0.001)
plt.plot(x,m,':',label=Z)
#plt.plot(x,mp,'--',label="fdp only")
plt.xlabel("Wavelength /A")
plt.ylabel("µ /(barns/atom)")
plt.legend()
plt.show()