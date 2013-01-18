import numpy as np
import matplotlib.pyplot as plt

# Get the data
u,g,r,i,z = np.loadtxt('/home/ant/Research/fp_proj/data/Quiescents_mag.csv',
                       delimiter=',',unpack=True) 
u0,g0,r0,i0,z0 =  np.loadtxt('/home/ant/Research/fp_proj/data/z0.025_0.100_mag.csv',
                             delimiter=',',unpack=True) 

all_x = r0
all_y = g0-r0

density = np.histogram2d(all_x, all_y, bins=40, range=[[-23,-17],[0,1]])

# define bin centers
X = 0.5 * (density[1][0:-1] + density[1][1:])
Y = 0.5 * (density[2][0:-1] + density[2][1:])
# parameters for plt.contour
Z = density[0].T
N = 10
levels = np.linspace(20,Z.max(),N)

# make the figure
plt.scatter(r,g-r,alpha=.4 ,c='r')
plt.contour(X, Y, Z, N, levels=levels)
plt.grid()
plt.ylabel('$g-r$',fontsize = 16)
plt.xlabel('$r$',fontsize = 16)
plt.ylim(0.2,1.0)
plt.xlim(-17.0,-23.0)

# save it
plt.savefig('color_mag.png')
