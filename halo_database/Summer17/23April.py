import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pynbody
import numpy as np
#import tangos as db
import math
from sys import exit, argv
from astropy.io import ascii

f_b = 0.17
h = 0.678

data = ascii.read('./AllRomulus.txt')
data['fgas'] *= f_b
Sun09 = ascii.read('Sun09scalingdata.txt')
Planck15 = ascii.read('Planck15scalingdata.txt')


#test to see ratio
sph = data["YSZ_sph"]
d01 = data["YSZ_cyl01"]
d1 = data["YSZ_cyl1"]
d10 = data["YSZ_cyl10"]
L = len(sph)
ratio01 = [sph[i]/d01[i] for i in xrange(L)]
ratio1 = [sph[i]/d1[i] for i in xrange(L)]
ratio10 = [sph[i]/d10[i] for i in xrange(L)]
print ratio01
print ratio1
print ratio10
#exit()


def scaling_YSZ_Planck(mass):
	M = 3e14 #MSun
	a_M = 5.0/3.0 #fitting exponent
	d = 68.123075601149 #distance to z = 0.0157 in Mpc
	Y_M = 0.73e-3 * math.pow((math.pi/(60*180)), 2.0) * d * d #0.73e-3 in arcmin^2, need to convert to Mpc^2
	Y_M *= 70
	data = [Y_M*math.pow(i/M, a_M) for i in mass]

	return data



MIN = 1e13/h
mlow = np.linspace(5e11/h, MIN, 100)
mhigh = np.linspace(MIN, 1e14/h, 100)
mhigher = np.linspace(1e14/h, 1e15/h, 100)
mall = np.linspace(5e11/h, 1e15/h, 100)


fit = scaling_YSZ_Planck(mall)
sph = data["YSZ_sph"]
d01 = data["YSZ_cyl01"]
d1 = data["YSZ_cyl1"]
d10 = data["YSZ_cyl10"]
plt.plot(mall, fit, 'k-', label='Planck+13')
plt.plot(data["M500"], sph, 'bo', label='sph, '+r'R$_{500}$')
#plt.plot(data["M500"], d01, 'o', markerfacecolor='none', label='cyl, 0.1 Mpc')
plt.plot(data["M500"], d1, 'yo', markerfacecolor='none', label='cyl, 1 Mpc')
plt.plot(data["M500"], d10, 'go', markerfacecolor='none', label='cyl, 10 Mpc')

plt.plot(Planck15["M500"], Planck15["YSZ"], 'go', label='Planck+15')

plt.xlabel(r"$M_{500}$(M$_{Sun}$)", fontsize=14)
plt.ylabel(r"YSZ (Mpc$^2$)", fontsize=14)
plt.xscale('log')
plt.yscale('log')
plt.legend()


plt.savefig('23April_YSZwithplanck.png')


plt.clf()
plt.close()

