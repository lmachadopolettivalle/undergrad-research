from sys import argv
from sys import exit
import matplotlib
matplotlib.use('Agg')
import pynbody
#import tangos as db
import matplotlib.pyplot as plt
import math
from astropy.io import ascii

if len(argv) != 2:
	print "Usage: ./maps.py OPTION"
	exit()
program, option = argv
haloid = 0
DIR = '/home/fas/nagai/lm643/halo_database/Summer17'

#Constants
k = 8.6173303 * 1e-8 #Boltzmann in keV K**-1
Zsun = 0.0134
mu = 8.0 / 7.0
m_p = 8.4089382 * 1e-58 #electron mass in solar masses
kpc3 = 2.93799895 * 1e64 #kpc**3 in cm**3


#R500 = [672.5, 375.15, 359.05, 332.35, 325.35, 304.45]
R500 = 672.5

#1: zoom in (5*R500)
#width = 2*5*R500, rad = 7*R500
#r500 = float(R500[haloid-1])
#rad = str(7*r500)+' kpc'
#width = str(10*r500)+' kpc'

#2: zoom out (8*R500)
#width = 2*8*R500, rad = 12*R500
rad = 12*R500
width = str(16*R500)+' kpc'


s = pynbody.load('/home/fas/nagai/lm643/project/h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096')

s.physical_units()

Center = (2546.87015086,   -97.94119966, -2971.89876526) #from database
#Center = (18666.1120692, 16872.9956994, 14925.66691555) #from s.halos() properties


halo = s.gas[pynbody.filt.Sphere(str(rad)+' kpc', cen=Center)]
halo['pos'] -= Center


#Density Map
if option == 'dens':
	#halo['rho'] *= (1.0/mu) * (1.0/m_p) * (1.0/kpc3)
	pynbody.plot.sph.image(halo.g, qty='rho', width=width, vmin=1000, vmax=1e8, units='Msol kpc^-2', filename=DIR+'/densitymap_'+str(haloid)+'.png', cmap='jet')
	exit(1)

#Temperature (mass weighted) Map
elif option == 'temp':
	halo['temp'] *= k
	pynbody.plot.sph.image(halo.g, qty='temp', width=width, vmin=0.0001, vmax=20, qtytitle=r'$log_{10} keV$', filename=DIR+'/tempmap_'+str(haloid)+'.png', cmap='jet')
	exit(1)

elif option ==	'tempmass':
	halo['temp'] *= k
	pynbody.plot.sph.image(halo.g, qty='temp', width=width, av_z='mass', vmin=0.0001, vmax=20, qtytitle=r'$log_{10} keV$', filename=DIR+'/tempmap_massav_'+str(haloid)+'.png', cmap='jet')
	exit(1)

#Metallicity (mass weighted) Map
elif option == 'met':
	halo['metals'] *= (1.0 / Zsun)
	pynbody.plot.sph.image(halo.g, qty='metals', width=width, vmin=0.001, vmax=10, qtytitle=r'$log_{10} Z_{sun}$', filename=DIR+'/metmap_'+str(haloid)+'.png', cmap='jet')
	exit(1)

elif option == 'metmass':
	halo['metals'] *= (1.0 / Zsun)
	pynbody.plot.sph.image(halo.g, qty='metals', width=width, av_z='mass', vmin=0.001, vmax=10, qtytitle=r'$log_{10} Z_{sun}$', filename=DIR+'/metmap_massav_'+str(haloid)+'.png', cmap='jet')
	exit(1)

elif option == 'entr':
	dens = halo['rho'] * (1.0/mu) * (1.0/m_p) * (1.0/kpc3)
	temp = halo['temp'] * k
	N = len(dens)
	halo['metals'] = [temp[i] / math.pow(dens[i], 2.0/3.0) for i in range(N)]
	pynbody.plot.sph.image(halo.g, qty='metals', width=width, vmin=0.01, vmax=1e6, qtytitle=r'$log_{10} keV cm^2$', filename=DIR+'/entrmap_'+str(haloid)+'.png', cmap='jet')
	exit(1)

elif option == 'entrmass':
	dens = halo['rho'] * (1.0/mu) * (1.0/m_p) * (1.0/kpc3)
	temp = halo['temp'] * k
	N = len(dens)
	halo['metals'] = [temp[i] / math.pow(dens[i], 2.0/3.0) for i in range(N)]
	pynbody.plot.sph.image(halo.g, qty='metals', width=width, av_z='mass', vmin=10, vmax=1e6, qtytitle=r'$log_{10} keV cm^2$', filename=DIR+'/entrmap_massav_'+str(haloid)+'.png', cmap='jet')
	exit(0)

elif option == 'XRay':
	dens = halo['rho'] * (1.0/mu) * (1.0/m_p) * (1.0/kpc3)
	temp = halo['temp'] * k
	dens2 = [i*i for i in dens]
	temp05 = [math.sqrt(i) for i in temp]
	N = len(dens)
	halo['metals'] = [dens2[i]*temp05[i] for i in xrange(N)]
	pynbody.plot.sph.image(halo.g, qty='metals', width=width, vmin=1e-24, vmax=1e-2, qtytitle=r'$log_{10} keV^{1/2} \cdot cm^{-6}$', filename=DIR+'/XRAY_'+str(haloid)+'.png', cmap='jet')
	exit(0)

elif option == 'TSZ':
	#mass density -> number density
	dens = halo['rho'] * (1.0/mu) * (1.0/m_p) * (1.0/kpc3)
	temp = halo['temp'] * k
	N = len(dens)
	halo['metals'] = [dens[i]*temp[i] for i in xrange(N)]
	pynbody.plot.sph.image(halo.g, qty='metals', width=width, qtytitle=r'$log_{10} keV \cdot cm^{-3}$', filename=DIR+'/TSZ_'+str(haloid)+'.png', cmap='jet', title='Thermal SZ - RomulusC')
	exit(0)

elif option == 'kSZ':
	#mass density -> number density
	dens = halo['rho'] * (1.0/mu) * (1.0/m_p) * (1.0/kpc3)
	N = len(dens)
	vel = halo['vel']
	v = [math.sqrt(math.pow(vel[i][0],2) + math.pow(vel[i][1],2) + math.pow(vel[i][2],2)) for i in xrange(N)]
	
	halo['metals'] = [dens[i]*v[i] for i in xrange(N)]
	pynbody.plot.sph.image(halo.g, qty='metals', width=width, qtytitle=r'$log_{10} km/s \cdot cm^{-3}$', filename=DIR+'/kSZ_'+str(haloid)+'.png', cmap='jet', title='Kinetic SZ - RomulusC')
	exit(0)

