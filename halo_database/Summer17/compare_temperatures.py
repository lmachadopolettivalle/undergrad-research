from sys import argv, exit
from astropy.io import ascii

if len(argv) != 4:
	print len(argv)
	print "Usage: python compare_temperatures.py HALOID FUNCTION HALF\n"
	print "HALOID is 0-31. FUNCTION is Mike, Luis or Vik. HALF is True or False.\n"
	exit(0)

program, haloid, function, half = argv
F = 1000 #decided by Nagai on 10/16/2017
# F = 1000 based on Zhuraleva et al. who compared F = 0, 500, 1000.

# Inputs
DIR = '/home/fas/nagai/lm643/project/'
SimDir = DIR+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'
if haloid != '0':
	SimDir = DIR+'cosmo25/cosmo25p.768sg1bwK1BHe75.008192'

Centers = ascii.read('RomulusCenters.txt')
Center = (Centers['Centerx'][int(haloid)], Centers['Centery'][int(haloid)], Centers['Centerz'][int(haloid)])

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pynbody
import numpy as np
import math
from yt.utilities.physical_constants import mh

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

def crit_dens(z=0):
	h = 0.678
	omegaL = 0.6914
	omegaM = 0.3086
	rhocrit0 = 127.131573017 # from Wiki, in Msun kpc**-3

	Ea = omegaM*math.pow((1+z), 3.0) + omegaL

	rhocrit = Ea * rhocrit0

	return rhocrit



#my estimate, from Mazzotta+04
def T_Luis(haloid, R500, F=500, f1=0.15, f2=1):
	T_CUT = 0.1
	crit = crit_dens()
	#rho_CUT = 500 * crit
	rho_CUT = F * crit

	s = pynbody.load(SimDir)
	s.physical_units()
	
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
		if F == 0:
			if k * temp1[i] > T_CUT:
				num1 += dens1[i]*math.pow(temp1[i], 0.25)
				den1 += dens1[i]/math.pow(temp1[i], 0.75)
				m1 += mass1[i]
		else:
			if k * temp1[i] > T_CUT and dens1[i] < rho_CUT:
				num1 += dens1[i]*math.pow(temp1[i], 0.25)
				den1 += dens1[i]/math.pow(temp1[i], 0.75)
				m1 += mass1[i]
	for i in xrange(N2):
		if F == 0:
			if k * temp2[i] > T_CUT:
				num2 += dens2[i]*math.pow(temp2[i], 0.25)
				den2 += dens2[i]/math.pow(temp2[i], 0.75)
				m2 += mass2[i]
		else:
			if k * temp2[i] > T_CUT and dens2[i] < rho_CUT:
				num2 += dens2[i]*math.pow(temp2[i], 0.25)
				den2 += dens2[i]/math.pow(temp2[i], 0.75)
				m2 += mass2[i]
	
	#print str(haloid)+" "+str(k * (num2-num1)/(den2-den1))+" "+str(m2-m1)
	#print "Hot Mass: ", m2-m1
	#return k * (num2-num1)/(den2-den1)
	return k * (num2-num1)/(den2-den1)





#helper function for T_Mike
def tcool(rho, T):
	C1 = 3.88*1e11 * 2.94e64
	C2 = 5*1e7

	return C1 * 0.59 * m_p * math.pow(T, 0.5) / (rho * (1 + C2*1/T))



#Michael's estimate
#This can be done either:
# i) from 0.15*R500 to R500
# ii) dr=1kpc around R500 or R2500
def T_Mike(haloid, R500, R2500, F=500, f1=0.15, f2=1, entropy_ring='0'):
	#T_CUT = 5e5 #in K
	T_CUT = 0.1 / k
	crit = crit_dens()
	#rho_CUT = 500 * crit
	rho_CUT = F * crit

	s = pynbody.load(SimDir)
	s.physical_units()
	
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
	for i in xrange(N1):
		T = temp1[i]
		rho = dens1[i]
		Tcool = tcool(rho, T)
		if F == 0:
			if T > T_CUT:
				em1 += rho * T / Tcool
				emT1 += rho * (T**2) / Tcool
		else:
			if T > T_CUT and rho < rho_CUT:
				em1 += rho * T / Tcool
				emT1 += rho * (T**2) / Tcool
	for i in xrange(N2):
		T = temp2[i]
		rho = dens2[i]
		Tcool = tcool(rho, T)
		if F == 0:
			if T > T_CUT:
				em2 += rho * T / Tcool
				emT2 += rho * (T**2) / Tcool
		else:
			if T > T_CUT and rho < rho_CUT:
				em2 += rho * T / Tcool
				emT2 += rho * (T**2) / Tcool


	#if entropy_ring != '0':
		#K = (k/1000) * (emT2-emT1)/(em2-em1) / (math.pow(n, 2.0/3.0))
		#print K
		#print "entropy to be done soon"
	#else:
		#print str(haloid)+" in keV: "+str((k) * (emT2-emT1)/(em2-em1))
	
	return (k) * (emT2-emT1)/(em2-em1)


def hunt (xx, n, x):
	INCREASING = 'True'
	if xx[n-1] < xx[0]:
		INCREASING = 'False'
		xx = xx[::-1]
	
	if x < xx[0]:
		if INCREASING == 'True':
			return 0
		else:
			return n-1

	i = 0
	while i < n:
		if xx[i] > x:
			if INCREASING == 'True':
				return i-1
			else:
				return n-i
	return n-1


def lin_interpolate (x,y,n,x0):
	y0 = 0
	ibin = 0

	if x0 <= x[0]:
		ibin = 0
	elif x0 >= x[n-1]:
		ibin = n-1
	else:
		if ibin < 1 or ibin > n-1:
			ibin=n/2
		ibin = hunt (x,n-1,x0)
		if ibin < 1 or ibin > n-1:
			exit(1)

	if x[ibin+1] == x[ibin]:
		y0 = y[ibin]
	else:
		y0 = y[ibin] + (y[ibin+1] - y[ibin]) * (x0 - x[ibin]) / (x[ibin+1] - x[ibin])
	
	return y0


def calc_Txspec(T, Z, E, N):
	data = ascii.read('/home/fas/nagai/lm643/mixT/tcal.dat')

	tcal = data['col1'] #T
	fcont_cal = data['col2'] #fluxcont
	fline_cal = data['col3'] #fluxline
	emean_cal = data['col4'] #emean
	fluxlinetot = data['col5'] #fluxlinetot
	logemean = data['col6'] #logemean
	ncal = len(tcal)
	
	acont = 0.875
	beta = 1
	delta1 = 0.19
	delta2 = 0.25

	wline = 0
	wcont = 0
	Tline = 0
	Tcont = 0
	fluxline = 0
	fluxcont = 0

	for i in xrange(N):
		#Find continuum and line flux, and average energy of the line emission
		print "hehehe"
		fcont = lin_interpolate (tcal,fcont_cal,ncal,T[i])
		fline = lin_interpolate (tcal,fline_cal,ncal,T[i])
		emean = lin_interpolate (tcal,emean_cal,ncal,T[i])
		#Multiply fluxes by emission measure and metallicity
		fcont = fcont * E[i]
		fline = fline * E[i] * Z[i]
		print "second"
		#eq.[6] in the paper
		Tcont = Tcont + T[i] * fcont/math.pow(T[i], acont)
		wcont = wcont + fcont/math.pow(T[i], acont)
		fluxcont = fluxcont + fcont
		print "almost all the way"
		#eq.[2] in the paper
		Tline = Tline + emean * fline
		wline = wline + fline
	
	print "end of for!"
	#eq.[4] in the paper
	Tcont = Tcont / wcont

	#eq.[3] in the paper
	emean = Tline / wline
	Tline = lin_interpolate (emean_cal,tcal,ncal,emean)
	fluxline = wline

	#eq.[7,8,12] in the paper
	x = fluxcont/(fluxcont+fluxline)
	wcont = math.exp(-(((x-1)**2)/delta1**2)**beta) * math.exp(-(((x-1)**2)/delta2**2)**4.0)
	Tmean = wcont*Tcont + (1-wcont)*Tline

	print Tmean
	exit(0)


#calc_Txspec([1],[1],[1],1)

#Using Vikhlinin's code estimate, from Mazzotta+04
# basically Mazzotta's method, but with different "a" (0.875) and with prefactor
# prefactor = c(T) = 3.49521971 for z=0.0157, NHI=4e20, emin=0.1, emax=10
def T_Vik(haloid, R500, F=500, f1=0.15, f2=1):
	print "started"
	T_CUT = 0.1
	crit = crit_dens()
	rho_CUT = F * crit

	s = pynbody.load(SimDir)
	s.physical_units()
	
	ring = s.gas[pynbody.filt.Annulus(f1*R500, f2*R500, cen=Center)]

	N = len(ring)

	dens = ring['rho'].in_units('g cm**-3')
	temp = ring['temp'].in_units('K')
	mass = ring['mass'].in_units('g')
	met = ring['metals']
	print "got data"
	D = []
	T = []
	E = []
	Z = []
	X_H = 0.76
	f = open('/home/fas/nagai/lm643/mixT/dataforweightT_'+str(haloid)+'.txt', 'w')
	for i in xrange(N):
		if k * temp[i] > T_CUT and dens[i] < rho_CUT:
			D.append(dens[i])
			T.append(k*temp[i])
			Z.append(met[i]/Zsun)

			#nenh = D[-1]/mh.value
			nenh = D[-1]/mh
			nenh *= nenh
			cell_v = mass[i]/D[-1]
			nenh *= 0.5*(1.+X_H)*X_H*cell_v
			E.append(nenh/2.93799894e64)
			
			f.write("%.5f %.5f %.10f\n" % (T[-1],Z[-1],E[-1]))
	f.close()


if function == 'Mike':
	data = ascii.read('500and2500info.txt')
	R500 = 0
	R2500 = 0
	
	R500 = data['R500'][int(haloid)]
	R2500 = data['R2500'][int(haloid)]

	if half == 'True':
		R500 *= 0.5
		R2500 *= 0.5

	print haloid+' '+str(T_Mike(haloid, R500, R2500, F))
	#if haloid == 'C' or haloid == '31':
	if haloid == '31':
		exit(0)
	exit(1)

if function == 'Luis':
	data = ascii.read('500and2500info.txt')
	R500 = 0

	R500 = data['R500'][int(haloid)]

	if half == 'True':
		R500 *= 0.5

	print haloid+' '+str(T_Luis(haloid, R500, F))
	#if haloid == 'C' or haloid == '31':
	if haloid == '31':
		exit(0)
	exit(1)

if function == 'Vik':
	data = ascii.read('500and2500info.txt')

	R500 = data['R500'][int(haloid)]

	if half == 'True':
		R500 *= 0.5
	T_Vik(haloid, R500, F)

	if haloid == '0' or haloid == '31':
		exit(0)
	exit(1)


