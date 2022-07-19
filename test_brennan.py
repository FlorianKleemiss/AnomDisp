from brennan import brennan
from constants_and_atomic_properties import *

#energy = float(input("Please enter Wavelength in Angstrom: "))
#z = input("Please Enter Element Name: ")
#energy = angstrom2eV / wavelength / 1000

Z = elements[60]
bren = brennan()
steps = []
minimal = 100
maximal = 3000
stepsize = 1
x = []
y = []
z = []
for step in range(int(minimal),int(maximal),int(stepsize)):
  x.append(step*0.001)
  temp = bren.at_angstrom(wavelength=step*0.001,element=Z)
  y.append(temp[0])
  z.append(temp[1])
import matplotlib.pyplot as plt
plt.plot(x,y,'+:',label="f'")
plt.plot(x,z,'+:',label="f\"")
plt.xlim(maximal*0.001,minimal*0.001)
plt.show()
