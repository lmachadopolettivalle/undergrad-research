import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sys import argv
from astropy.io import ascii
import numpy as np
import math

if len(argv) != 2:
	print "Usage: ./plotTemperatures.py PLOT(nocut, midcut, highcut, halfcut)"
	exit()

program, PLOT = argv

h = 0.678

data = ascii.read('MikeLuisVikTemperatures.txt')

DATA = {'Luis_no': data['Luis_nocut'], 'Luis_mid': data['Luis_midcut'], 'Luis_high': data['Luis_highcut'], 'Mike_no': data['Mike_nocut'], 'Mike_mid': data['Mike_midcut'], 'Mike_high': data['Mike_highcut'], 'Luis_half': data['Luis_halfcut'], 'Mike_half': data['Mike_halfcut'], 'Vik_high': data['Vik_highcut'], 'Vik_half': data['Vik_halfcut']}


mass = ascii.read('500and2500info.txt')

MASS = {'Mgas500': mass['Mgas500'], 'Mtot500': mass['Mtot500']}

whichMass = 'Mtot500'


#r = [Mike[i]/Luis[i] for i in xrange(6)]
#plt.plot(Luis, r, 'o')

def scaling_Sun09(mass, a, V):
	M = 1.26e14 / h # MSun

	data = [V*math.pow(i/M, 1.0/a) for i in mass]
	
	return data

MIN = 1e13/h
mlow = np.linspace(5e11/h, MIN, 100)
mhigh = np.linspace(MIN, 1e14/h, 100)
mhigher = np.linspace(1e14/h, 2e14/h, 100)

a = 1.65
V = 3.0
fitlow = scaling_Sun09(mlow, a, V)
fithigh = scaling_Sun09(mhigh, a, V)
fithigher = scaling_Sun09(mhigher, a, V)

plt.plot(mhigher, fithigher, 'k--')
plt.plot(mhigh, fithigh, 'k-', label='Sun+09')
plt.plot(mlow, fitlow, 'k--', label='Sun+09 (extrapolated)')

plt.ylabel(r'$T (keV)$', fontsize=16)
plt.xlabel(r'$M_{500}$', fontsize=16)
plt.xscale('log')
plt.ylim([0.1, 2])
plt.yscale('log')

if PLOT == 'nocut':
	plt.plot(MASS[whichMass], DATA['Luis_no'], 'bo', label='Mazzotta+04')
	plt.plot(MASS[whichMass], DATA['Mike_no'], 'rs', label='Emission Weighted')
	plt.title(r'T vs. $M_{500}$ - no density cut', fontsize=16)
	plt.legend()
	plt.savefig('nocut.png')
elif PLOT == 'midcut':
	plt.plot(MASS[whichMass], DATA['Luis_mid'], 'bo', label='Mazzotta+04')
	plt.plot(MASS[whichMass], DATA['Mike_mid'], 'rs', label='Emission Weighted')
	plt.title(r'T vs. $M_{500}$: $\rho < 500\cdot \rho_c$', fontsize=16)
	plt.legend()
	plt.savefig('midcut.png')
elif PLOT == 'highcut':
	plt.plot(MASS[whichMass], DATA['Luis_high'], 'bo', label='Mazzotta+04')
	plt.plot(MASS[whichMass], DATA['Mike_high'], 'rs', label='Emission Weighted')
	plt.plot(MASS[whichMass], DATA['Vik_high'], 'gv', label='Vikhlinin+05')
	plt.title(r'T vs. $M_{500}$: $\rho < 1000\cdot \rho_c$', fontsize=16)
	plt.legend()
	plt.savefig('highcut.png')
elif PLOT == 'halfcut': #compare halfcuts to highcuts, as decided on October 16
	plt.plot(MASS[whichMass], DATA['Luis_half'], 'bo', label='Mazzotta+04')
	plt.plot(MASS[whichMass], DATA['Mike_half'], 'rs', label='Emission Weighted')
	plt.plot(MASS[whichMass], DATA['Vik_half'], 'gv', label='Vikhlinin+05')
	plt.title(r'T vs. $M_{500}$: $\rho < 1000\cdot \rho_c$, R $< 0.5*R_{500}$', fontsize=16)
	plt.legend()
	plt.savefig('halfcut.png')

