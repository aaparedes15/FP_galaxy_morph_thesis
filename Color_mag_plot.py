# This program will make a color (g-r) magnitude (M_r) diagram for the selected galaxies in our sample.

# import needed modules
import numpy as np
import matplotlib.pyplot as plt

# import the data
u,g,r,i,z = np.loadtxt('/home/ant/Research/fp_proj/data/Quiescents_mag.csv',delimiter=',',unpack=True) # selected galaxies
u0,g0,r0,i0,z0 =  np.loadtxt('/home/ant/Research/fp_proj/data/z0.025_0.100_mag.csv',delimiter=',',unpack=True) # all the data
# plot the thing
plt.scatter(r0,g0-r0,c='b',alpha=0.4,label='All Galaxies')
plt.scatter(r,g-r,c='r',label='Chosen Data')
plt.grid(True)
plt.legend(loc=4)
plt.ylabel('$g-r$')
plt.xlabel('$M_r$')
plt.ylim(0.0,1.0)
plt.xlim(-17.0,-23.0)
plt.savefig('/home/ant/Documents/Research/fp_proj/hwgn/figs/color_mag.png')
plt.savefig('/home/ant/Documents/Research/fp_proj/hwgn/report/present/color_mag.png')
plt.cla()
