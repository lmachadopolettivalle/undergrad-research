from sys import argv
from astropy.io import ascii

if len(argv) < 2:
	print "Usage: python scaling.py HALOID\n"
	print "Also, remember to vim scaling.py and choose what function will be used - by commenting out the other functions at the bottom\n"
	exit()

program, haloid = argv

# Inputs
DIR = '/home/fas/nagai/lm643/project/'
SimDir = DIR+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'
if haloid != 'C':
	SimDir = DIR+'cosmo25/cosmo25p.768sg1bwK1BHe75.008192'

Centers = ascii.read('RomulusCenters.txt')
Center = 0
if haloid == 'C':
	Center = (Centers['Centerx'][0], Centers['Centery'][0], Centers['Centerz'][0])
elif haloid in ['1', '2', '3', '4', '5']:
	Center = (Centers['Centerx'][int(haloid)], Centers['Centery'][int(haloid)], Centers['Centerz'][int(haloid)])
#center = {'C': (2546.87015086, -97.94119966, -2971.89876526), '1': (7446.7740934127805, 8478.3427413597001, 7965.3717622399463), '2': (5888.7597183651769, -3920.9626619405344, 9269.1704642761524), '3': (5338.6136956516802, -9639.365229590805, 7312.5757339895217), '4': (11731.934544501413, 10293.593262519084, 10325.876620841929), '5': (1976.5326767473034, -6942.4867710401186, 4635.9382550076662)}
#Center = center[haloid]

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pynbody
import numpy as np
#import tangos as db
import math
from sys import exit

#sigma_T = 0.665 * 1e-24 #cm**2
sigma_T = 6.98426558 * 1e-74 #Mpc**2
m_ec2 = 511.0 #keV
h = 0.678
k = 8.6173303 * 1e-8 #Boltzmann in keV K**-1
Zsun = 0.0134
mu = 8.0 / 7.0
m_p = 8.4089382 * 1e-58 #proton mass in solar masses
kpc3 = 2.93799895 * 1e64 #kpc**3 in cm**3
f = 500 	#for R500
f_b = 0.17	#cosmological bayon fraction


#from Nagai et al. 2007
#keV
def T500(M500, z=0):
	omegaL = 0.6914
	omegaM = 0.3086
	rhocrit0 = 127.131573017 # from Wiki, in Msun kpc**-3
	Ea = omegaM*math.pow((1+z), 3.0) + omegaL

	m = (M500 * h * 1e-15)
	
	return 11.05 * math.pow(m, 2.0/3.0) * math.pow(Ea, 2.0/3.0)
	
#keV cm**2
def K500(M500, z=0):
	omegaL = 0.6914
	omegaM = 0.3086
	rhocrit0 = 127.131573017 # from Wiki, in Msun kpc**-3
	Ea = omegaM*math.pow((1+z), 3.0) + omegaL

	m = (M500 * h * 1e-15)

	return 1963 * math.pow(m, 2.0/3.0) / math.pow(Ea, 2.0/3.0)

def crit_dens(z=0):
	h = 0.678
	omegaL = 0.6914
	omegaM = 0.3086
	rhocrit0 = 127.131573017 # from Wiki, in Msun kpc**-3

	Ea = omegaM*math.pow((1+z), 3.0) + omegaL

	rhocrit = Ea * rhocrit0

	return rhocrit


#find M500 for given haloid
#haloid=1 is main halo
def M500(haloid=1, F=500, step=1):
	data = ascii.read('./mass_C.txt')
	gas_mass = data["gas"]
	dm_mass = data["dm"]
	star_mass = data["star"]
	rbins = data["rbins"]
	N = len(rbins)
	z = 0
	tot_mass = [gas_mass[i]+star_mass[i]+dm_mass[i] for i in xrange(N)]
	#tot_mass = H.calculate('tot_mass_profile')
	#gas_mass = H.calculate('gas_mass_profile')
	#star_mass = H.calculate('star_mass_profile')
	#rbins =  H.calculate_abcissa_values('tot_mass_profile')
	
	j = 0
	N = len(rbins)
	crit = crit_dens(z)
	while j < N and (3 * tot_mass[j] / (4*math.pi*math.pow(rbins[j], 3.0))) > F * crit:
		#print 3 * tot_mass[j] / (4*math.pi*math.pow(rbins[j], 3.0)) / crit
		j += 1
	#exit()	
	#returned:	R500	Mgas		M500		fgas
	return [rbins[j], gas_mass[j], tot_mass[j], gas_mass[j]/(f_b * tot_mass[j])]

#print M500(F=2500)
#exit()

def YSZ(haloid=1, R500=672.5, option='sph', L=0):
	s = pynbody.load(SimDir)
	s.physical_units()
	
	#Center = center[haloid]
	if option == 'sph':
		halo = s.gas[pynbody.filt.Sphere(str(R500)+' kpc', cen=Center)]
	elif option == 'cyl':
		halo = s.gas[pynbody.filt.Disc(str(R500)+' kpc', str(L/2)+' Mpc', cen=Center)]
		print "L = ", L
	mass = halo['mass'] * (1.0/mu) * (1.0/m_p)  
	#dens = halo['rho']
	temp = halo['temp']
	N = len(halo)
	Mass = 0
	T = 0
	MT = 0
	for i in xrange(N):
		Mass += mass[i]
		T += k*temp[i]
		MT += mass[i]*k*temp[i]
	
	print "C ", (sigma_T/m_ec2)*MT
	exit()
	return 0







#my estimate, from Mazzotta+04
def Tspec(haloid, R500, f1=0.15, f2=1):
	CUT = 0.1

	s = pynbody.load(SimDir)
	s.physical_units()
	
	Center = center[haloid]
	
	halo1 = s.gas[pynbody.filt.Sphere(str(f1*R500)+' kpc', cen=Center)]
	halo2 = s.gas[pynbody.filt.Sphere(str(f2*R500)+' kpc', cen=Center)]

	N1 = len(halo1)
	N2 = len(halo2)

	dens1 = halo1['rho']
	dens2 = halo2['rho']
	temp1 = halo1['temp']
	temp2 = halo2['temp']
	mass1 = halo1['mass']
	mass2 = halo2['mass']

	num1 = 0
	num2 = 0
	den1 = 0
	den2 = 0
	m1 = 0
	m2 = 0
	for i in xrange(N1):
		if k * temp1[i] > CUT:
			num1 += dens1[i]*math.pow(temp1[i], 0.25)
			den1 += dens1[i]/math.pow(temp1[i], 0.75)
			m1 += mass1[i]
	for i in xrange(N2):
		if k * temp2[i] > CUT:
			num2 += dens2[i]*math.pow(temp2[i], 0.25)
			den2 += dens2[i]/math.pow(temp2[i], 0.75)
			m2 += mass2[i]
	
	print str(haloid)+" "+str(k * (num2-num1)/(den2-den1))+" "+str(m2-m1)
	#print "Hot Mass: ", m2-m1
	#return k * (num2-num1)/(den2-den1)
	return 0

def tcool(rho, T):
	C1 = 3.88*1e11 * 2.94e64
	C2 = 5*1e7

	return C1 * 0.59 * m_p * math.pow(T, 0.5) / (rho * (1 + C2*1/T))

#Michael's estimate
#This can be done either:
# i) from 0.15*R500 to R500 or R2500
# ii) dr=1kpc around R500 or R2500
def T_Mike(haloid, R500, R2500, f1=0.15, f2=1, entropy_ring='0'):
	#CUT = 5e5 #in K
	CUT = 0.1 / k

	s = pynbody.load(SimDir)
	s.physical_units()
	
	#Center = center[haloid]

	r1 = f1*R500
	r2 = f2*R500
	if entropy_ring == '500':
		r1 = R500 - 1
		r2 = R500 + 1
	elif entropy_ring == '2500':
		r1 = R2500 - 1
		r2 = R2500 + 1


	halo1 = s.gas[pynbody.filt.Sphere(str(r1)+' kpc', cen=Center)]
	halo2 = s.gas[pynbody.filt.Sphere(str(r2)+' kpc', cen=Center)]

	N1 = len(halo1)
	N2 = len(halo2)

	dens1 = halo1['rho']
	dens2 = halo2['rho']
	temp1 = halo1['temp']
	temp2 = halo2['temp']
	#mass1 = halo1['mass']
	#mass2 = halo2['mass']

	em1 = 0
	em2 = 0
	emT1 = 0
	emT2 = 0
	crit = crit_dens()
	for i in xrange(N1):
		T = temp1[i]
		rho = dens1[i]
		Tcool = tcool(rho, T)
		if T > CUT and rho < 500 * crit:
		#if T > CUT:
			em1 += rho * T / Tcool
			emT1 += rho * (T**2) / Tcool
	for i in xrange(N2):
		T = temp2[i]
		rho = dens2[i]
		Tcool = tcool(rho, T)
		if T > CUT and rho < 500 * crit:
		#if T > CUT:
			em2 += rho * T / Tcool
			emT2 += rho * (T**2) / Tcool


	if entropy_ring != '0':
		#K = (k/1000) * (emT2-emT1)/(em2-em1) / (math.pow(n, 2.0/3.0))
		#print K
		print "entropy to be done soon"
	else:
		print str(haloid)+" in keV: "+str((k) * (emT2-emT1)/(em2-em1))
	
	return 0
"""
data = ascii.read('500and2500info.txt')
R500 = 0
R2500 = 0
if haloid == 'C':
	R500 = data['R500'][0]
	R2500 = data['R2500'][0]
else:
	R500 = data['R500'][int(haloid)]
	R2500 = data['R2500'][int(haloid)]

T_Mike(haloid, R500, R2500)
if haloid == 'C' or haloid == '31':
	exit(0)
exit(1)
"""



#WITHOUT TEMPERATURE CUT
TSPECS = [0.00167099906569, 0.00557755099381, 0.24142546818594204, 0.001948260745492353, 0.00233285475624127, 0.001672193555942395, 0.007300485045825053, 0.17639107195594234, 0.001620174219408887, 0.002773139918885975, 0.0016191916433756009, 0.0013115779807522715, 0.00176948866126, 0.00130391376848, 0.00204635170793, 0.0020186744662469909, 0.157791815546, 0.13684163547889383, 0.00134304008588, 0.0176793156385, 0.00126131660901, 0.0111824193352, 0.0471436927976, 0.0013149560207, 0.0013327655132, 0.0455617935842, 0.00272283280422, 0.00779220292752, 0.00124489864133, 0.00124808565022, 0.0011171941755]



def YX(haloid):
	data = ascii.read('./AllRomulus.txt')
	Mgasring = float(data["Mgasring"][haloid])
	Mhot = float(data["Mhot"][haloid])
	Tspec = float(data["Tspec"][haloid])

	#return Mhot * Tspec
	return Mgasring * Tspec


def scaling_Sun09(mass, a, V):
	if a > 1: #Tspec
		M = 1.26e14 / h # MSun
	else: #YX
		M = 1.14e14 / h # MSun

	data = [V*math.pow(i/M, 1.0/a) for i in mass]
	
	return data

def scaling_YSZ_Planck(mass):
	M = 3e14 #MSun
	a_M = 5.0/3.0 #fitting exponent
	d = 68.123075601149 #distance to z = 0.0157 in Mpc
	Y_M = 0.73e-3 * math.pow((math.pi/(60*180)), 2.0) * d * d #0.73e-3 in arcmin^2, need to convert to Mpc^2
	Y_M *= 70
	data = [Y_M*math.pow(i/M, a_M) for i in mass]

	return data


def plot_vsmass():
	program, option = argv
	if option is None:
		exit()
	if option in ['Mgas', 'fgas', 'Tspec', 'YX', 'YSZ']:
		data = ascii.read('./AllRomulus.txt')
		data['fgas'] *= f_b
		Sun09 = ascii.read('Sun09scalingdata.txt')
		Planck15 = ascii.read('Planck15scalingdata.txt')
	else:
		print "Call plot function with option in [Mgas, fgas, Tspec, YX, YSZ]"
		exit()

	MIN = 1e13/h
	mlow = np.linspace(5e11/h, MIN, 100)
	mhigh = np.linspace(MIN, 1e14/h, 100)
	mhigher = np.linspace(1e14/h, 1e15/h, 100)
	mall = np.linspace(5e11/h, 1e15/h, 100)

	if option == 'YSZ':
		fit = scaling_YSZ_Planck(mall)
		sph = data["YSZ_sph"]
		d01 = data["YSZ_cyl01"]
		d1 = data["YSZ_cyl1"]
		d10 = data["YSZ_cyl10"]
		plt.plot(mall, fit, 'r-', label='Planck+13')
		plt.plot(data["M500"], sph, 'bo', label='sph')
		#plt.plot(data["M500"], d01, 'o', markerfacecolor='none', label='cyl, 0.1 Mpc')
		plt.plot(data["M500"], d1, 'yo', markerfacecolor='none', label='cyl, 1 Mpc')
		plt.plot(data["M500"], d10, 'go', markerfacecolor='none', label='cyl, 10 Mpc')

		plt.plot(Planck15["M500"], Planck15["YSZ"], 'go', label='Planck+15')

		plt.xlabel(r"$M_{500}$(M$_{Sun}$)", fontsize=14)
		plt.ylabel(r"YSZ (Mpc$^2$)", fontsize=14)
		plt.xscale('log')
		plt.yscale('log')
		plt.legend()
		#plt.savefig('./allYSZ.png')
		#plt.savefig('./someYSZ.png')
		plt.savefig('YSZwithplanck.png')
		exit()
	

	MIN = 1e13/h
	mlow = np.linspace(5e11/h, MIN, 100)
	mhigh = np.linspace(MIN, 1e14/h, 100)
	mhigher = np.linspace(1e14/h, 1e15/h, 100)

	plt.plot(data["M500"], data[option], 'o', label='Romulus')
	plt.xlabel(r"$M_{500}$(M$_{Sun}$)", fontsize=14)
	plt.xscale('log')
	if option == 'Mgas':
		plt.plot(Sun09["M500"], Sun09["Mgas"], 'ro', label='Sun+09')
		plt.plot(Planck15["M500"], Planck15["Mgas"], 'go', label='Planck+15')
		plt.yscale('log')
		plt.legend()
		plt.ylabel(r'Gas Mass (M$_{Sun}$)', fontsize=14)
	elif option == 'fgas':
		plt.plot(Sun09["M500"], Sun09["fgas"], 'ro', label='Sun+09')
		plt.plot(Planck15["M500"], Planck15["fgas"], 'go', label='Planck+15')
		plt.yscale('linear')
		plt.ylabel(r'f$_{gas}$', fontsize=14)
		plt.legend()
		#plt.ylim([0, 1.0])
	elif option == 'Tspec':
		#Sun+09: 1.65, 3.0
		a = 1.65
		V = 3.0
		fitlow = scaling_Sun09(mlow, a, V)
		fithigh = scaling_Sun09(mhigh, a, V)
		fithigher = scaling_Sun09(mhigher, a, V)
		
		plt.plot(mhigher, fithigher, 'r--')
		plt.plot(mhigh, fithigh, 'r-', label='Sun+09')
		plt.plot(mlow, fitlow, 'r--', label='Sun+09 (extrapolated)')
		plt.plot(Planck15["M500"], Planck15["Tspec"], 'go', label='Planck+15')
		plt.legend()
		plt.yscale('log')
		plt.ylabel(r'T$_X$ (keV)', fontsize=14)
	elif option == 'YX':
		#Sun+09: 0.571, 4.0e13/h
		a = 0.571
		V = 4.0e13/h
		fitlow = scaling_Sun09(mlow, a, V)
		fithigh = scaling_Sun09(mhigh, a, V)
		fithigher = scaling_Sun09(mhigher, a, V)
		
		plt.plot(mhigher, fithigher, 'r--')
		plt.plot(mhigh, fithigh, 'r-', label='Sun+09')
		plt.plot(mlow, fitlow, 'r--', label='Sun+09 (extrapolated)')
		plt.plot(Planck15["M500"], Planck15["YX"], 'go', label='Planck+15')
		plt.legend()
		
		plt.yscale('log')
		plt.ylabel(r'YX (M$_{Sun}$ $\cdot$ keV)', fontsize=14)

	#plt.savefig('./Scaling/'+option+'.png')
	plt.savefig(option+'.png')
	plt.clf()
	plt.close()

plot_vsmass()

#writeall(np.arange(1, 100))


