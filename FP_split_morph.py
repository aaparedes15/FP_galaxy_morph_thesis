import matplotlib.pyplot as plt
import scipy as scipy
import numpy as np
import csv

# import the data from csv files

FP = csv.reader(open("/home/ant/Research/fp_proj/data/Quiescents_FP.csv"))
EW = csv.reader(open("/home/ant/Research/fp_proj/data/Quiescents_EW.csv"))
ID = csv.reader(open("/home/ant/Research/fp_proj/data/Quiescents_ID.csv"))
AB = csv.reader(open("/home/ant/Research/fp_proj/data/Quiescents_AB.csv"))
Luc = csv.reader(open("/home/ant/Research/fp_proj/data/Quiescents_Luc.csv"))
petro = csv.reader(open("/home/ant/Research/fp_proj/data/Quiescents_petro.csv"))
petro_phy = csv.reader(open("/home/ant/Research/fp_proj/data/Quiescents_petro_phy.csv"))
morph = csv.reader(open("/home/ant/Research/fp_proj/data/Quiescents_morph.csv"))

# put the data into seperate tuples

Re, Re_err, Ie, Ie_err, Vdisp, Vdisp_err = zip(*FP)
EW_Halpha, EW_Halpha_err, EW_OII, EW_OII_err = zip(*EW)
RA, Dec, plate, MJD, fiber = zip(*ID)
AB_deV, AB_deV_err = zip(*AB)
BT, s2 = zip(*Luc)
R50, R90, cons = zip(*petro)
r50_phy, r90_phy = zip(*petro_phy)
Sersic_n, fracDeV = zip(*morph)

# put the data into seperate arrays

Re = np.array(Re, dtype='float')
Re_err = np.array(Re_err, dtype='float')
Ie = np.array(Ie, dtype='float')
Ie_err = np.array(Ie_err, dtype='float')
Vdisp = np.array(Vdisp, dtype='float')
Vdisp_err = np.array(Vdisp_err, dtype='float')

EW_Halpha = np.array(EW_Halpha, dtype='float')
EW_Halpha_err = np.array(EW_Halpha_err, dtype='float')
EW_OII = np.array(EW_OII, dtype='float')
EW_OII_err = np.array(EW_OII_err, dtype='float')

RA = np.array(RA, dtype='float')
Dec = np.array(Dec, dtype='float')
plate = np.array(plate, dtype='float')
MJD = np.array(MJD, dtype='float')
fiber = np.array(fiber, dtype='float')

AB_deV = np.array(AB_deV, dtype='float')
AB_deV_err = np.array(AB_deV_err, dtype='float')

BT = np.array(BT, dtype='float')
s2 = np.array(s2, dtype='float')

R50 = np.array(R50, dtype='float')
R90 = np.array(R90, dtype='float')
cons = np.array(cons, dtype='float')

r50_phy = np.array(r50_phy, dtype='float')
r90_phy = np.array(r90_phy, dtype='float')

Sersic_n = np.array(Sersic_n, dtype='float')
fracDeV = np.array(fracDeV, dtype='float')

Index = np.arange(len(RA))

# Now get a best fit plane of logIe vs. logRe & logVdisp

from scipy import optimize

logVdisp = np.log10(Vdisp)
logRe = np.log10(Re)
logIe = np.log10(Ie)
logIe_err = Ie_err / Ie

fitfunc = lambda p, x1, x2: p[0] + p[1]*x1 + p[2]*x2
errfunc = lambda p, x1, x2, y, err: (y - fitfunc(p,x1,x2))/err

guess = [1.0, -1.0, -1.0]

out = optimize.leastsq(errfunc, guess, args=(logRe, logVdisp, logIe, logIe_err), full_output = 1)

coeff = out[0]
covar = out[1]

print "Plane of best fit:   log(Ie) = {0} log(Re) + {1} log(Vdisp) + {2}".format(round(coeff[1],2),round(coeff[2],2),round(coeff[0],2))

# get a line of best fit

best_fit = coeff[0] + coeff[1]*logRe + coeff[2]*logVdisp

xR_range = max(logRe) - min(logRe)
xV_range = max(logVdisp) - min(logVdisp)

xR_step = 0.1
xV_step = xR_step/(xR_range/xV_range)

xR = np.arange(min(logRe),max(logRe),xR_step)
xV = np.arange(min(logVdisp),max(logVdisp),xV_step)

best_fit_line = coeff[0] + coeff[1]*xR + coeff[2]*xV

# from now on logIe ---> Ie,  logRe ---> Re,   logVdisp ---> Vdisp


A = BT

not_B = np.where(A <= 0.5)[0]
B_Index = np.delete(Index, not_B)
B = np.delete(s2, not_B)
B_AB_deV = np.delete(AB_deV, not_B)
B_cons = np.delete(cons, not_B)

not_K = np.where(B > 0.08)[0]
K_Index = np.delete(B_Index, not_K)
K = np.delete(B_AB_deV, not_K)
KB_cons = np.delete(B_cons, not_K)

not_D = np.where(K <= 0.65)[0]
D_Index = np.delete(K_Index, not_D)
D = np.delete(KB_cons, not_D)

not_Bulge = np.where(D <= 2.9)[0]
Bulge = np.delete(D_Index, not_Bulge)  # Indecies of Bulge Galaxies

dummy_K = (K > 0.45) & (K <= 0.65)
Int_K = K_Index[np.where(dummy_K)]

Int_D = D_Index[np.where(D <= 2.29)]

Intermediate = np.concatenate([Int_K, Int_D]) # Indecies of Intermediate Galaxies

Disc = np.delete(Index, np.concatenate([Bulge,Intermediate]))  # Indecies of Disc Galaxies

not_Disc = np.delete(Index, Disc)
s2_Disc = np.delete(s2, not_Disc)

Unsmooth_Disc = Disc[np.where(s2_Disc > 0.08)[0]] # Indecies of Unsmooth_Disc Galaxies

Smooth_Disc = Disc[np.where(s2_Disc <= 0.08)[0]] # Indecies of Smooth_Disc Galaxies
'''
loc_disc = np.array([RA[Disc],DEC[Disc]]).T
loc_int = np.array([RA[Intermediate],DEC[Intermediate]]).T
loc_bulge = np.array([RA[Bulge],DEC[Bulge]]).T

np.savetxt('loc_Jdisc.dat',loc_disc,fmt='%10.6f')
np.savetxt('loc_Jintermediat.dat',loc_int,fmt='%10.6f')
np.savetxt('loc_Jbulge.dat',loc_bulge,fmt='%10.6f')
'''
# Get FP parameters for each morphological type

Re_disc = np.delete(logRe, np.concatenate([Bulge,Intermediate]))
Re_int = np.delete(logRe, np.concatenate([Bulge,Disc]))
Re_bulge = np.delete(logRe, np.concatenate([Intermediate,Disc]))

Ie_disc = np.delete(logIe, np.concatenate([Bulge,Intermediate]))
Ie_int = np.delete(logIe, np.concatenate([Bulge,Disc]))
Ie_bulge = np.delete(logIe, np.concatenate([Intermediate,Disc]))

Vdisp_disc = np.delete(logVdisp, np.concatenate([Bulge,Intermediate]))
Vdisp_int = np.delete(logVdisp, np.concatenate([Bulge,Disc]))
Vdisp_bulge = np.delete(logVdisp, np.concatenate([Intermediate,Disc]))

best_fit_disc = np.delete(best_fit, np.concatenate([Bulge,Intermediate]))
best_fit_int = np.delete(best_fit, np.concatenate([Bulge,Disc]))
best_fit_bulge = np.delete(best_fit, np.concatenate([Intermediate,Disc]))

# Split each morph type into low, mid and high SB. call them a, b, and c respectively

# disc

Delta_Ie_disc = Ie_disc - best_fit_disc

a_cuts = (Delta_Ie_disc >= -0.06) | (Delta_Ie_disc <= -0.18)
not_a = np.where(a_cuts)[0]
a_Ie_disc = np.delete(Ie_disc, not_a)
a_Re_disc = np.delete(Re_disc, not_a)
a_Vdisp_disc = np.delete(Vdisp_disc, not_a)
a_plate = np.delete(plate, not_a)
a_MJD = np.delete(MJD, not_a)
a_fiber = np.delete(fiber, not_a)
a_Index = np.delete(Index, not_a)

b_cuts = (Delta_Ie_disc < -0.06) | (Delta_Ie_disc > 0.06)
not_b = np.where(b_cuts)[0]
b_Ie_disc = np.delete(Ie_disc, not_b)
b_Re_disc = np.delete(Re_disc, not_b)
b_Vdisp_disc = np.delete(Vdisp_disc, not_b)
b_plate = np.delete(plate, not_b)
b_MJD = np.delete(MJD, not_b)
b_fiber = np.delete(fiber, not_b)
b_Index = np.delete(Index, not_b)

c_cuts = (Delta_Ie_disc <= 0.06) | (Delta_Ie_disc >= 0.18)
not_c = np.where(c_cuts)[0]
c_Ie_disc = np.delete(Ie_disc, not_c)
c_Re_disc = np.delete(Re_disc, not_c)
c_Vdisp_disc = np.delete(Vdisp_disc, not_c)
c_plate = np.delete(plate, not_c)
c_MJD = np.delete(MJD, not_c)
c_fiber = np.delete(fiber, not_c)
c_Index = np.delete(Index, not_c)

# intermediate

Delta_Ie_int = Ie_int - best_fit_int

a_cuts = (Delta_Ie_int >= -0.06) | (Delta_Ie_int <= -0.18)
not_a = np.where(a_cuts)[0]
a_Ie_int = np.delete(Ie_int, not_a)
a_Re_int = np.delete(Re_int, not_a)
a_Vdisp_int = np.delete(Vdisp_int, not_a)
a_plate = np.delete(plate, not_a)
a_MJD = np.delete(MJD, not_a)
a_fiber = np.delete(fiber, not_a)
a_Index = np.delete(Index, not_a)

b_cuts = (Delta_Ie_int < -0.06) | (Delta_Ie_int > 0.06)
not_b = np.where(b_cuts)[0]
b_Ie_int = np.delete(Ie_int, not_b)
b_Re_int = np.delete(Re_int, not_b)
b_Vdisp_int = np.delete(Vdisp_int, not_b)
b_plate = np.delete(plate, not_b)
b_MJD = np.delete(MJD, not_b)
b_fiber = np.delete(fiber, not_b)
b_Index = np.delete(Index, not_b)

c_cuts = (Delta_Ie_int <= 0.06) | (Delta_Ie_int >= 0.18)
not_c = np.where(c_cuts)[0]
c_Ie_int = np.delete(Ie_int, not_c)
c_Re_int = np.delete(Re_int, not_c)
c_Vdisp_int = np.delete(Vdisp_int, not_c)
c_plate = np.delete(plate, not_c)
c_MJD = np.delete(MJD, not_c)
c_fiber = np.delete(fiber, not_c)
c_Index = np.delete(Index, not_c)

# bulge

Delta_Ie_bulge = Ie_bulge - best_fit_bulge

a_cuts = (Delta_Ie_bulge >= -0.06) | (Delta_Ie_bulge <= -0.18)
not_a = np.where(a_cuts)[0]
a_Ie_bulge = np.delete(Ie_bulge, not_a)
a_Re_bulge = np.delete(Re_bulge, not_a)
a_Vdisp_bulge = np.delete(Vdisp_bulge, not_a)
a_plate = np.delete(plate, not_a)
a_MJD = np.delete(MJD, not_a)
a_fiber = np.delete(fiber, not_a)
a_Index = np.delete(Index, not_a)

b_cuts = (Delta_Ie_bulge < -0.06) | (Delta_Ie_bulge > 0.06)
not_b = np.where(b_cuts)[0]
b_Ie_bulge = np.delete(Ie_bulge, not_b)
b_Re_bulge = np.delete(Re_bulge, not_b)
b_Vdisp_bulge = np.delete(Vdisp_bulge, not_b)
b_plate = np.delete(plate, not_b)
b_MJD = np.delete(MJD, not_b)
b_fiber = np.delete(fiber, not_b)
b_Index = np.delete(Index, not_b)

c_cuts = (Delta_Ie_bulge <= 0.06) | (Delta_Ie_bulge >= 0.18)
not_c = np.where(c_cuts)[0]
c_Ie_bulge = np.delete(Ie_bulge, not_c)
c_Re_bulge = np.delete(Re_bulge, not_c)
c_Vdisp_bulge = np.delete(Vdisp_bulge, not_c)
c_plate = np.delete(plate, not_c)
c_MJD = np.delete(MJD, not_c)
c_fiber = np.delete(fiber, not_c)
c_Index = np.delete(Index, not_c)

# graph the FP with colors representing the diffrent morph types
# red:bulge, green:intermediate, blue:disc

plt.plot(best_fit_disc, Ie_disc, 'b.', markersize=2)
plt.plot(best_fit_int, Ie_int, 'g.', markersize=2)
plt.plot(best_fit_bulge, Ie_bulge, 'r.', markersize=2)
plt.plot(best_fit_line, -0.18 + best_fit_line, 'k-',
         best_fit_line, -0.06 + best_fit_line, 'k-',
        best_fit_line, 0.06 + best_fit_line, 'k-',
         best_fit_line, 0.18 + best_fit_line, 'k-')
plt.axis([1.5, 4.0, 1.0, 4.0])
plt.xlabel('{0} log $R_e$ + {1} log $\sigma$ + {2}'.format(round(coeff[1], 2), round(coeff[2], 2), round(coeff[0], 2)))
plt.ylabel('log $I_e$')
plt.title('Red:Bulge, Green:Intermediate, Blue:Disc')
#plt.savefig('/home/ant/Documents/Research/fp_proj/figs/FP_morph.png')
plt.show()
plt.cla()

# now graph the FP all split up with the colors

x = np.arange(1.85,2.45,0.15)
y = np.arange(-0.3,1.2,0.001)

plt.subplot(3,3,1)
plt.plot(a_Vdisp_disc,a_Re_disc,'b.', markersize=2)
plt.plot(x, x*0 - 0.3, 'k-', x, x*0,
         'k-', x, x*0 + 0.3, 'k-', x, x*0 + 0.6, 'k-', x, x*0 + 0.9,
         'k-', x, x*0 +1.2 , 'k-', y*0 + 1.85, y, 'k-', y*0 + 2.00, y,
         'k-', y*0 + 2.15, y, 'k-', y*0 + 2.30, y, 'k-', y*0 + 2.45, y,'k-')
plt.ylim(-0.5,1.25)
plt.ylabel(' $R_e$ (kpc)')
plt.title('low-SB')

plt.subplot(3,3,2)
plt.plot(b_Vdisp_disc,b_Re_disc,'b.', markersize=2)
plt.plot(x, x*0 - 0.3, 'k-', x, x*0,
         'k-', x, x*0 + 0.3, 'k-', x, x*0 + 0.6, 'k-', x, x*0 + 0.9,
         'k-', x, x*0 +1.2 , 'k-', y*0 + 1.85, y, 'k-', y*0 + 2.00, y,
         'k-', y*0 + 2.15, y, 'k-', y*0 + 2.30, y, 'k-', y*0 + 2.45, y,'k-')
plt.ylim(-0.5,1.25)
plt.title('midplane')

plt.subplot(3,3,3)
plt.plot(c_Vdisp_disc,c_Re_disc,'b.', markersize=2)
plt.plot(x, x*0 - 0.3, 'k-', x, x*0,
         'k-', x, x*0 + 0.3, 'k-', x, x*0 + 0.6, 'k-', x, x*0 + 0.9,
         'k-', x, x*0 +1.2 , 'k-', y*0 + 1.85, y, 'k-', y*0 + 2.00, y,
         'k-', y*0 + 2.15, y, 'k-', y*0 + 2.30, y, 'k-', y*0 + 2.45, y,'k-')
plt.ylim(-0.5,1.25)
plt.title('high-SB')

plt.subplot(3,3,4)
plt.plot(a_Vdisp_int,a_Re_int,'g.', markersize=2)
plt.plot(x, x*0 - 0.3, 'k-', x, x*0,
         'k-', x, x*0 + 0.3, 'k-', x, x*0 + 0.6, 'k-', x, x*0 + 0.9,
         'k-', x, x*0 +1.2 , 'k-', y*0 + 1.85, y, 'k-', y*0 + 2.00, y,
         'k-', y*0 + 2.15, y, 'k-', y*0 + 2.30, y, 'k-', y*0 + 2.45, y,'k-')
plt.ylim(-0.5,1.25)
plt.ylabel(' $R_e$ (kpc)')

plt.subplot(3,3,5)
plt.plot(b_Vdisp_int,b_Re_int,'g.', markersize=2)
plt.plot(x, x*0 - 0.3, 'k-', x, x*0,
         'k-', x, x*0 + 0.3, 'k-', x, x*0 + 0.6, 'k-', x, x*0 + 0.9,
         'k-', x, x*0 +1.2 , 'k-', y*0 + 1.85, y, 'k-', y*0 + 2.00, y,
         'k-', y*0 + 2.15, y, 'k-', y*0 + 2.30, y, 'k-', y*0 + 2.45, y,'k-')
plt.ylim(-0.5,1.25)

plt.subplot(3,3,6)
plt.plot(c_Vdisp_int,c_Re_int,'g.', markersize=2)
plt.plot(x, x*0 - 0.3, 'k-', x, x*0,
         'k-', x, x*0 + 0.3, 'k-', x, x*0 + 0.6, 'k-', x, x*0 + 0.9,
         'k-', x, x*0 +1.2 , 'k-', y*0 + 1.85, y, 'k-', y*0 + 2.00, y,
         'k-', y*0 + 2.15, y, 'k-', y*0 + 2.30, y, 'k-', y*0 + 2.45, y,'k-')
plt.ylim(-0.5,1.25)

plt.subplot(3,3,7)
plt.plot(a_Vdisp_bulge,a_Re_bulge,'r.', markersize=2)
plt.plot(x, x*0 - 0.3, 'k-', x, x*0,
         'k-', x, x*0 + 0.3, 'k-', x, x*0 + 0.6, 'k-', x, x*0 + 0.9,
         'k-', x, x*0 +1.2 , 'k-', y*0 + 1.85, y, 'k-', y*0 + 2.00, y,
         'k-', y*0 + 2.15, y, 'k-', y*0 + 2.30, y, 'k-', y*0 + 2.45, y,'k-')
plt.ylim(-0.5,1.25)
plt.xlabel(' $\sigma$ (km/s)')
plt.ylabel(' $R_e$ (kpc)')

plt.subplot(3,3,8)
plt.plot(b_Vdisp_bulge,b_Re_bulge,'r.', markersize=2)
plt.plot(x, x*0 - 0.3, 'k-', x, x*0,
         'k-', x, x*0 + 0.3, 'k-', x, x*0 + 0.6, 'k-', x, x*0 + 0.9,
         'k-', x, x*0 +1.2 , 'k-', y*0 + 1.85, y, 'k-', y*0 + 2.00, y,
         'k-', y*0 + 2.15, y, 'k-', y*0 + 2.30, y, 'k-', y*0 + 2.45, y,'k-')
plt.ylim(-0.5,1.25)
plt.xlabel(' $\sigma$ (km/s)')

plt.subplot(3,3,9)
plt.plot(c_Vdisp_bulge,c_Re_bulge,'r.', markersize=2)
plt.plot(x, x*0 - 0.3, 'k-', x, x*0,
         'k-', x, x*0 + 0.3, 'k-', x, x*0 + 0.6, 'k-', x, x*0 + 0.9,
         'k-', x, x*0 +1.2 , 'k-', y*0 + 1.85, y, 'k-', y*0 + 2.00, y,
         'k-', y*0 + 2.15, y, 'k-', y*0 + 2.30, y, 'k-', y*0 + 2.45, y,'k-')
plt.ylim(-0.5,1.25)
plt.xlabel(' $\log{\sigma}$ (km/s)')
plt.show()
#plt.savefig('/home/ant/Documents/Research/fp_proj/figs/FP_split_morph.png')








