
import numpy as np
import matplotlib.pyplot as plt
import sys
import linecache
r_wire = 0.125 # in mm


dat = sys.argv[1]

nelectrons = linecache.getline(dat, 14)[25]
print('Processing {} electrons'.format(nelectrons))

if dat == None:
    dat = "av_out.dat"

ELdat = np.loadtxt(dat, skiprows=17)

X_av = ELdat[:,0]*10
Y_av = ELdat[:,1]*10
Z_av = ELdat[:,2]*10
T_av = ELdat[:,3]
type_av = ELdat[:,4]
lvl_av = ELdat[:,5]

fig = plt.figure(1)
fig.set_figheight(5.0)
fig.set_figwidth(8.0)

ax = fig.add_subplot(111)
ax.plot(-Z_av[type_av == 4],Y_av[type_av == 4],'*',c='blue',markersize=1.5)
ax.plot(-Z_av[type_av == -1],Y_av[type_av == -1],'--',c='gray',alpha=0.2)
ax.set_xlim(-0.5,1.5)
ax.set_ylim(-0.5, 1.0)
ax.set_xlabel('z (mm)')
ax.set_ylabel('y (mm)')
ax.axvspan(-r_wire, r_wire, color='gray')
plt.text(-0.35, 0.42, "Buffer region", size=14, rotation=90.)
plt.text(0.8, 0.2, "Active region", size=14, rotation=0.)
plt.text(0.9, -0.4, "{} excitations".format(len(Z_av[type_av == 4])), size=14, rotation=0.,color='blue')
plt.savefig("EL_wire.pdf", bbox_inches='tight')
