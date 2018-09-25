#################################
# 								#
# write_h5.py					#
# 								#
# Create.hdf5 file with			#
# regular grid data				#
#								#
#################################
import h5py
N = 256
rad = 1000
L = rad / 2
width = str(2*L)+' kpc'

old = h5py.File("/home/fas/nagai/lm643/halo_database/Summer17/MockXray/RomCParticles.hdf5", "r")
pos = old['position'].value

FILE = str(rad)+'kpcData.hdf5'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pynbody
#import h5py
#import tangos as db
import numpy as np
from scipy.interpolate import griddata
from sys import argv

program, option = argv

grid_x, grid_y, grid_z = np.mgrid[-1*L:L:N*1j, -1*L:L:N*1j, -1*L:L:N*1j]

print "about to enter ifs and elifs"
if option == 'density':
	f = h5py.File("/home/fas/nagai/lm643/halo_database/Summer17/MockXray/"+FILE, "w")
	dens = old['density'].value
	griddens = griddata(pos, dens, (grid_x, grid_y, grid_z), method='linear')
	f.create_dataset("density",data=griddens)
	exit(0)

f = h5py.File("/home/fas/nagai/lm643/halo_database/Summer17/MockXray/"+FILE, "a")

if option == 'temp':
	temp = old['temperature'].value
	gridtemp = griddata(pos, temp, (grid_x, grid_y, grid_z), method='linear')
	f.create_dataset("temperature",data=gridtemp)
	exit(1)
elif option == 'met':
	met = old['metallicity'].value
	gridmet = griddata(pos, met, (grid_x, grid_y, grid_z), method='linear')
	f.create_dataset("metallicity",data=gridmet)
	exit(1)
elif option == 'velx':
	velx = old['x-velocity'].value
	gridvelx = griddata(pos, velx, (grid_x, grid_y, grid_z), method='linear')
	f.create_dataset("x-velocity",data=gridvelx)
	exit(1)
elif option == 'vely':
	vely = old['y-velocity'].value
	gridvely = griddata(pos, vely, (grid_x, grid_y, grid_z), method='linear')
	f.create_dataset("y-velocity",data=gridvely)
	exit(0)
elif option == 'velz':
	velz = old['z-velocity'].value
	gridvelz = griddata(pos, velz, (grid_x, grid_y, grid_z), method='linear')
	f.create_dataset("z-velocity",data=gridvelz)
	exit(1)
elif option == 'mass':
	mass = old['mass'].value
	gridmass = griddata(pos, mass, (grid_x, grid_y, grid_z), method='linear')
	f.create_dataset("mass",data=gridmass)
	exit(0)


