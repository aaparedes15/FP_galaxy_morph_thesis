import numpy as np
import matplotlib.pyplot as plt
from pylab import *

# import the needed data
D_spot, D_Age_Hb = np.loadtxt('/home/ant/Research/fp_proj/hwgn/data/Disc_spec.dat',unpack=True,usecols=(0,1),delimiter=',')
B_spot, B_Age_Hb = np.loadtxt('/home/ant/Research/fp_proj/hwgn/data/Bulge_spec.dat',unpack=True,usecols=(0,1),delimiter=',')
D_spot_err, D_age_plus_err, D_age_minus_err = np.loadtxt('/home/ant/Research/fp_proj/hwgn/data/Disc_spec.dat_err',unpack=True,usecols=(0,1,2),delimiter=',')
B_spot_err, B_age_plus_err, B_age_minus_err = np.loadtxt('/home/ant/Research/fp_proj/hwgn/data/Bulge_spec.dat_err',unpack=True,usecols=(0,1,2),delimiter=',')

# import the fp data for bulge type galaxies
bulge = csv2rec('/home/ant/Research/fp_proj/pros/data4GG/Bulge_data4GG.csv')

# split the data into high, mid and low Surface Brightnesses
low = (bulge['delta_logie'] < -0.06)
mid = (bulge['delta_logie'] > -0.06) & (bulge['delta_logie'] < 0.06)
high = (bulge['delta_logie'] > 0.06)

low_bulge = bulge[np.where(low)]
mid_bulge = bulge[np.where(mid)]
high_bulge = bulge[np.where(high)]

# import the stellar population data
#B_Age, B_Fe_H, B_Mg_Fe, B_C_Fe, B_N_Fe, B_Ca_Fe, B_O_Fe, B_Na_Fe, B_SB_Fe, B_Cr_Fe, B_TB_Fe, B_iso, B_x_imf, B_Age_Hb, B_Age_Hg, B_Age_Hd, B_Spot = np.loadtxt('/home/ant/Documents/Research/fp_proj/data/stack_spec_dat/Bulge_spec_dat.csv',delimiter='\t',unpack=True)


# repeat everything from the top, but for Disc galaxies.
disc = csv2rec('/home/ant/Research/fp_proj/pros/data4GG/Disc_data4GG.csv')

# split the data into high, mid and low Surface Brightnesses
low = (disc['delta_logie'] < -0.06)
mid = (disc['delta_logie'] > -0.06) & (disc['delta_logie'] < 0.06)
high = (disc['delta_logie'] > 0.06)

low_disc = disc[np.where(low)]
mid_disc = disc[np.where(mid)]
high_disc = disc[np.where(high)]

# import the stellar population data
#D_Age, D_Fe_H, D_Mg_Fe, D_C_Fe, D_N_Fe, D_Ca_Fe, D_O_Fe, D_Na_Fe, D_SB_Fe, D_Cr_Fe, D_TB_Fe, D_iso, D_x_imf, D_Age_Hb, D_Age_Hg, D_Age_Hd, D_Spot = np.loadtxt('/home/ant/Documents/Research/fp_proj/data/stack_spec_dat/Disc_spec_dat.csv',delimiter='\t',unpack=True)

Delta_Age_Hb = np.array([])
Spot = np.array([])
Age_Disc = np.array([])
Age_Bulge = np.array([])

for i in np.arange(len(B_spot)):
	for j in np.arange(len(D_spot)):
		 if B_spot[i] == D_spot[j]:
			Delta_Age_Hb = np.append(Delta_Age_Hb, B_Age_Hb[i] - D_Age_Hb[j])
			Age_Disc = np.append(Age_Disc,D_Age_Hb[j])
			Age_Bulge = np.append(Age_Bulge, B_Age_Hb[i])
			Spot = np.append(Spot,B_spot[i])

Low = np.concatenate((low_bulge,low_disc))
Mid = np.concatenate((mid_bulge,mid_disc))
High = np.concatenate((high_bulge,high_disc))

low_dAge_Hb = np.array([])
mid_dAge_Hb = np.array([])
high_dAge_Hb = np.array([])

low_re = np.array([])
mid_re = np.array([])
high_re = np.array([])

low_vdisp = np.array([])
mid_vdisp = np.array([])
high_vdisp = np.array([])

low_Disc_age_Hb = np.array([])
mid_Disc_age_Hb = np.array([])
high_Disc_age_Hb = np.array([])

low_Bulge_age_Hb = np.array([])
mid_Bulge_age_Hb = np.array([])
high_Bulge_age_Hb = np.array([])

for i in np.arange(len(Spot)):
	for j in np.arange(len(Low['position'])):
		if Spot[i] == Low['position'][j]:
			low_dAge_Hb = np.append(low_dAge_Hb, Delta_Age_Hb[i])
			low_re = np.append(low_re, Low['logre'][j])
			low_vdisp = np.append(low_vdisp, Low['logvdisp'][j])
			low_Disc_age_Hb = np.append(low_Disc_age_Hb, Age_Disc[i])
			low_Bulge_age_Hb = np.append(low_Bulge_age_Hb, Age_Bulge[i])
	for j in np.arange(len(Mid['position'])):
		if Spot[i] == Mid['position'][j]:
			mid_dAge_Hb = np.append(mid_dAge_Hb, Delta_Age_Hb[i])
			mid_re = np.append(mid_re, Mid['logre'][j])
			mid_vdisp = np.append(mid_vdisp, Mid['logvdisp'][j])
			mid_Disc_age_Hb = np.append(mid_Disc_age_Hb, Age_Disc[i])
			mid_Bulge_age_Hb = np.append(mid_Bulge_age_Hb, Age_Bulge[i])
	for j in np.arange(len(High['position'])):
		if Spot[i] == High['position'][j]:
			high_dAge_Hb = np.append(high_dAge_Hb, Delta_Age_Hb[i])
			high_re = np.append(high_re, High['logre'][j])
			high_vdisp = np.append(high_vdisp, High['logvdisp'][j])
			high_Disc_age_Hb = np.append(high_Disc_age_Hb, Age_Disc[i])
			high_Bulge_age_Hb = np.append(high_Bulge_age_Hb, Age_Bulge[i])


low_re = np.delete(low_re, np.where(np.isnan(low_dAge_Hb)))
mid_re = np.delete(mid_re, np.where(np.isnan(mid_dAge_Hb)))
high_re = np.delete(high_re, np.where(np.isnan(high_dAge_Hb)))

low_vdisp = np.delete(low_vdisp, np.where(np.isnan(low_dAge_Hb)))
mid_vdisp = np.delete(mid_vdisp, np.where(np.isnan(mid_dAge_Hb)))
high_vdisp = np.delete(high_vdisp, np.where(np.isnan(high_dAge_Hb)))

Age_Disc = np.delete(Age_Disc, np.where(np.isnan(Delta_Age_Hb)))
Age_Bulge = np.delete(Age_Bulge, np.where(np.isnan(Delta_Age_Hb))) 

Delta_Age_Hb = np.delete(Delta_Age_Hb, np.where(np.isnan(Delta_Age_Hb)))

low_Disc_Age_Hb = np.delete(low_Disc_age_Hb, np.where(np.isnan(low_dAge_Hb)))
mid_Disc_Age_Hb = np.delete(mid_Disc_age_Hb, np.where(np.isnan(mid_dAge_Hb)))
high_Disc_Age_Hb = np.delete(high_Disc_age_Hb, np.where(np.isnan(high_dAge_Hb)))

low_Bulge_Age_Hb = np.delete(low_Bulge_age_Hb, np.where(np.isnan(low_dAge_Hb)))
mid_Bulge_Age_Hb = np.delete(mid_Bulge_age_Hb, np.where(np.isnan(mid_dAge_Hb)))
high_Bulge_Age_Hb = np.delete(high_Bulge_age_Hb, np.where(np.isnan(high_dAge_Hb)))

low_dAge_Hb = np.delete(low_dAge_Hb, np.where(np.isnan(low_dAge_Hb)))
mid_dAge_Hb = np.delete(mid_dAge_Hb, np.where(np.isnan(mid_dAge_Hb)))
high_dAge_Hb = np.delete(high_dAge_Hb, np.where(np.isnan(high_dAge_Hb)))

subplot(3,3,1)
scatter(low_vdisp,low_re,vmin=3,vmax=12,c=low_Disc_Age_Hb,cmap='bwr')
xlim(1.8,2.5); ylim(-.2,1.0)
ylabel('Disc Age\n\nlog $R_e$ [kpc]',fontsize=16)
title('low-SB',fontsize=16)
grid(True)

subplot(3,3,2)
scatter(mid_vdisp,mid_re,vmin=3,vmax=12,c=mid_Disc_Age_Hb,cmap='bwr')
xlim(1.8,2.5); ylim(-.2,1.0)
title('mid-SB',fontsize=16)
grid(True)

subplot(3,3,3)
scatter(high_vdisp,high_re,vmin=3,vmax=12,c=high_Disc_Age_Hb,cmap='bwr')
colorbar()
xlim(1.8,2.5); ylim(-.2,1.0)
title('high-SB',fontsize=16)
grid(True)

subplot(3,3,4)
scatter(low_vdisp,low_re,vmin=Age_Bulge.min(),vmax=12,c=low_Bulge_Age_Hb,cmap='bwr')
xlim(1.8,2.5); ylim(-.2,1.0)
ylabel('Bulge Age\n\nlog $R_e$ [kpc]',fontsize=16)
grid(True)

subplot(3,3,5)
scatter(mid_vdisp,mid_re,vmin=3,vmax=12,c=mid_Bulge_Age_Hb,cmap='bwr')
xlim(1.8,2.5); ylim(-.2,1.0)
grid(True)

subplot(3,3,6)
scatter(high_vdisp,high_re,vmin=3,vmax=12,c=high_Bulge_Age_Hb,cmap='bwr')
colorbar()
xlim(1.8,2.5); ylim(-.2,1.0)
grid(True)

subplot(3,3,7)
scatter(low_vdisp,low_re,vmin=-3,vmax=3,c=low_dAge_Hb,cmap='bwr')
xlim(1.8,2.5); ylim(-.2,1.0)
ylabel('Age$_{Bulge}$ - Age$_{Disc}$ \n\n log $R_e$ [kpc]',fontsize=16);xlabel('log $\sigma$ [km/s]',fontsize=16)
grid(True)

subplot(3,3,8)
scatter(mid_vdisp,mid_re,vmin=-3,vmax=3,c=mid_dAge_Hb,cmap='bwr')
xlim(1.8,2.5); ylim(-.2,1.0)
xlabel('log $\sigma$ [km/s]',fontsize=16)
grid(True)

subplot(3,3,9)
scatter(high_vdisp,high_re,vmin=-3,vmax=3,c=high_dAge_Hb,cmap='bwr')
colorbar()
xlim(1.8,2.5); ylim(-.2,1.0)
xlabel('log $\sigma$ [km/s]',fontsize=16)
grid(True)
'''
plt.hist(Delta_Age_Hb)
plt.title('Age$_{Disc}$ - Age$_{Bulge}$')
plt.annotate('$\mu$ =%10.3f\n$\sigma$ =%10.3f\n$\sigma_{mean}$ =%10.3f' % (np.mean(Delta_Age_Hb),np.std(Delta_Age_Hb),np.std(Delta_Age_Hb)/(len(Delta_Age_Hb))),(-4.5,9.5))
plt.show()

'''

show()
