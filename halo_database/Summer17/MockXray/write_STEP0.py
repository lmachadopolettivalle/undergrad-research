#################################
# 								#
# write_particlesh5.py			#
# 								#
# Makes STEP0					#
# Create.hdf5 file with			#
# particle data					#
#								#
#################################
from sys import argv, exit

if len(argv) != 3:
	print "Usage: ./write_STEP0.py HALOID CORE(YES or NO)"
	exit(0)

program, haloid, core = argv

print "*** Starting haloid =", haloid, "***"

R500 = {'C': 672.5, '1': 375.15, '2': 359.05, '3': 332.35, '4': 325.35, '5': 304.45}
R200 = {'C': 1200.0, '1': 551.55, '2': 549.05, '3': 488.35, '4': 481.85, '5': 451.75}

rad = 2*R200[haloid]
L = rad / 2
width = str(2*L)+' kpc'
#N = 128
#rad = 50
#L = rad / 2
#width = str(2*L)+' kpc'

FILE = core+'CORE_'+str(rad)+'kpc_Rom'+haloid+'Particles.hdf5'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pynbody
import h5py
#import tangos as db
import numpy as np
from astropy.io import ascii
#import yt

print 'Importing... done'

sigma_T = 6.98426558 * 1e-74 #Mpc**2
m_ec2 = 511.0 #keV
h = 0.678
k = 8.6173303 * 1e-8 #Boltzmann in keV K**-1
Zsun = 0.0134
mu = 8.0 / 7.0
m_p = 8.4089382 * 1e-58 #proton mass in solar masses
kpc3 = 2.93799895 * 1e64 #kpc**3 in cm**3
f = 500         #for R500
f_b = 0.17      #cosmological bayon fraction


DIR = '/home/fas/nagai/lm643/project/'
SimDir = DIR+'cosmo25/cosmo25p.768sg1bwK1BHe75.008192'
if haloid == 'C':
	SimDir = DIR+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'

s = pynbody.load(SimDir)
s.physical_units()

print 'pynbody.load()... done'

CenterFile = '/home/fas/nagai/lm643/halo_database/Summer17/RomulusCenters.txt'
i = 0
if haloid != 'C':
	i = int(haloid)
Centers = ascii.read(CenterFile)
Center = (Centers['Centerx'][i], Centers['Centery'][i], Centers['Centerz'][i])


halo = 0

if core == 'YES':
	halo = s.gas[pynbody.filt.Sphere(str(rad)+' kpc', cen=Center)]
	print "Number of particles within ", rad, "kpc is", len(halo)
elif core == 'NO':
	r1 = 20 #in kpc - I tested this on Oct 22; Annulus assumes kpc
	halo = s.gas[pynbody.filt.Annulus(r1, rad, cen=Center)]
	print "Number of particles between", r1, 'kpc and', rad, "kpc is", len(halo)
else:
	print core, " is neither YES nor NO!"
	exit(0)

print 'pynbody.filt... done'

halo['pos'] -= Center

print 'halo[pos] -= Center... done'

#pynbody.plot.sph.image(halo.g, qty='rho', width=width, units='g cm^-2', filename='/home/fas/nagai/lm643/halo_database/Summer17/MockXray/'+FILE+'.png', cmap='jet')
#exit()
f = h5py.File("/home/fas/nagai/lm643/halo_database/Summer17/MockXray/STEP0/"+FILE, "a")

pos = halo['pos'] #kpc
f.create_dataset("position",data=pos)
print 'create POSITION... done'
dens = halo['rho'] #Msun/kpc**3
f.create_dataset("density",data=dens)
print 'create DENSITY... done'
temp = halo['temp'] #Kelvin
f.create_dataset("temperature",data=temp)
print 'create TEMPERATURE... done'
met = halo['metals'] * (1.0/Zsun)
f.create_dataset("metallicity",data=met)
print 'create METALLICITY... done'
vel = halo['vel'] #km/s
velx = [i[0] for i in vel]
f.create_dataset("x-velocity",data=velx)
print 'create VELZ... done'
vely = [i[1] for i in vel]
f.create_dataset("y-velocity",data=vely)
print 'create VELY... done'
velz = [i[2] for i in vel]
f.create_dataset("z-velocity",data=velz)
print 'create VELZ... done'
mass = halo['mass'] #Msun
f.create_dataset("mass",data=mass)
print 'create MASS... done'

if haloid == '5':
	exit(0)
exit(1)
