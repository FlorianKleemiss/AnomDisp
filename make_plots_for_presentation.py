from playing_with_coefficients import *
from constants_and_atomic_properties import *
el_nr = 42
t0 = 0.48 * math.pi
alp = 0.25 * math.pi

fig, axes = plt.subplots(1,1)
z_temp = np.linspace(0.800, 5.000, 250)
s1t = -2.0*apply_angle_part_s_parallel(f_s_1_hoenl(z_temp)  *constant_factor     ,t0, alp, 1)
for i,e in enumerate(z_temp):
  if e < 1.0:
    s1t[i] = 0
z_temp = z_temp*get_ionization_energy_1s(el_nr)/1000
axes.plot(z_temp,-s1t,color='black',label="f''")
axes.legend()
axes.set_xlabel('Energy (in keV)')
axes.axhline(y=0,linestyle='dashed',color='gray')
fig.savefig("K-edge.png",dpi=300,transparent=False,bbox_inches='tight')

x = np.linspace(1.0001, 5.0001, 100)
z = x*get_ionization_energy_1s(el_nr)/1000
x12 = xn(x*get_ionization_energy_1s(el_nr)/h,1,el_nr,0,2)
s1 = -2*apply_angle_part_s_parallel(f_s_1_hoenl(x)  *constant_factor     ,t0, alp, 1)
s2 = -2*apply_angle_part_s_parallel(f_s_2_1_hoenl(x)*constant_factor*x12,t0, alp, 1) + apply_angle_part_s_parallel(f_s_2_2_hoenl(x)*constant_factor*x12,t0, alp, 2)
s3 = numpy.zeros_like(s2)

for i,zl in enumerate(x):
  for l in range(15):
    res_4 = f_a_for_p(el_nr,l,0,zl,zl*get_ionization_energy_1s(el_nr)/h,1,4)
    for runny in range(len(res_4)):
      res_4[runny] = apply_angle_part_s_parallel(res_4[runny].real,t0,alp, l)
    s3[i] += 2*sum(res_4)


fig, axes = plt.subplots(1,1)
axes.scatter(z,-s1,s=10,facecolors='none',edgecolors='b',marker='^',label="K-shell p=0")# type:ignore
axes.scatter(z,-s2,s=10,facecolors='none',edgecolors='g',marker='^',label="K-shell p=2")# type:ignore
axes.scatter(z,-s3,s=10,facecolors='none',edgecolors='r',marker='^',label="K-shell p=4")# type:ignore
axes.plot(z,-s1-s2-s3,color='black',label="Sum")# type:ignore
axes.legend()# type:ignore
axes.set_xlabel('Energy (in keV)')# type:ignore
axes.axhline(y=0,linestyle='dashed',color='gray')# type:ignore
fig.savefig("s-electrons.png",dpi=300,transparent=False,bbox_inches='tight')
#plt.show()

x = np.linspace(1.0001, 10.000, 100)
z2 = x*get_ionization_energy_2s(el_nr)/1000
z3 = x*get_ionization_energy_2p1_2(el_nr)/1000
x221 = xn(x*get_ionization_energy_2s(el_nr)/h,2,el_nr,0,2)
x222 = xn(x*get_ionization_energy_2p1_2(el_nr)/h,2,el_nr,0,2)

s1 = 2*apply_angle_part_s_parallel(f_s_1_EM(x)  *constant_factor     ,t0, alp, 1)
s2 = 2*(apply_angle_part_s_parallel(f_s_2_1_EM(x)*constant_factor*x221,t0, alp, 1) + apply_angle_part_s_parallel(f_s_2_2_EM(x)*constant_factor*x221,t0, alp, 2))
s3 = 6*apply_angle_part_s_parallel(f_p_1_EM(x)  *constant_factor     ,t0, alp, 1)
s4 = 6*(apply_angle_part_s_parallel(f_p_2_1_EM(x)*constant_factor*x222,t0, alp, 1) + apply_angle_part_s_parallel(f_s_2_2_EM(x)*constant_factor*x222,t0, alp, 2))


fig, axes = plt.subplots(1,1)
axes.scatter(z2,s1,s=10,facecolors='none',edgecolors='b',marker='^',label="L-shell s-electron p=0")# type:ignore
axes.scatter(z2,s2,s=10,facecolors='none',edgecolors='g',marker='^',label="L-shell s-electron p=2")# type:ignore
axes.scatter(z3,s3,s=10,facecolors='none',edgecolors='magenta',marker='o',label="L-shell p-electron p=0")# type:ignore
axes.scatter(z3,s4,s=10,facecolors='none',edgecolors='cyan',marker='o',label="L-shell p--electron p=2")# type:ignore
axes.plot(z2,(s1+s2),color='black',label="Sum s-electrons")# type:ignore
axes.plot(z3,(s3+s4),color='orange',label="Sum p-electrons")# type:ignore
axes.legend()# type:ignore
axes.set_xlabel('Energy (in keV)')# type:ignore
axes.axhline(y=0,linestyle='dashed',color='gray')# type:ignore
fig.savefig("p-electrons.png",dpi=300,transparent=False,bbox_inches='tight')
#plt.show()