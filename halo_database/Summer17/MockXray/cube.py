import numpy as np
import h5py
from sys import argv
from sys import exit

if len(argv) != 4:
	print "Usage: cube.py HALOID OPTION CORE(YES or NO)"
	exit(0)
program, haloid, OPTION, core = argv

R500 = {'C': 672.5, '1': 375.15, '2': 359.05, '3': 332.35, '4': 325.35, '5': 304.45}
R200 = {'C': 1200.0, '1': 551.55, '2': 549.05, '3': 488.35, '4': 481.85, '5': 451.75}

rad = 2*R200[haloid]
L = rad / 2
width = str(2*L)+' kpc'

BOX_SIZE = L #box goes from -1*BOX_SIZE to +1*BOX_SIZE; in kpc
NUM_PIXELS = 512 #number of pixels from -L to +L

PARTICLE_FILE = str(rad)+'kpc_Rom'+haloid+'Particles.hdf5'
DATA_FILE = str(rad)+'kpc_Rom'+haloid+'Data.hdf5'
if core == 'NO':
	PARTICLE_FILE = core + 'CORE_' + PARTICLE_FILE
	DATA_FILE = core + 'CORE_' + DATA_FILE

#h5py.File("/home/fas/nagai/lm643/halo_database/Summer17/MockXray/"+DATA_FILE, "w")
#exit()

def find_index(x, L=BOX_SIZE, N=NUM_PIXELS):
	return int(np.floor((x+L)*N/(2.*L)))

#props = ['density', 'temperature', 'metallicity', 'x-velocity', 'y-velocity', 'z-velocity', 'mass']
def main(PROPERTY):
	f = h5py.File("/home/fas/nagai/lm643/halo_database/Summer17/MockXray/STEP0/"+PARTICLE_FILE, "r")
	
	data = np.zeros((NUM_PIXELS, NUM_PIXELS, NUM_PIXELS))
	weight = np.ones((NUM_PIXELS, NUM_PIXELS, NUM_PIXELS))
	part_property = f[PROPERTY].value
	part_pos = f['position'].value
	
	NUM_PARTICLES = len(part_property)
	part_weight = f['mass'].value
	if PROPERTY == 'density' or PROPERTY == 'mass':
		part_weight = [1]*NUM_PARTICLES

	for p in xrange(NUM_PARTICLES):
		x = part_pos[p][0]
		y = part_pos[p][1]
		z = part_pos[p][2]
		w = part_weight[p]
		if x > BOX_SIZE or x < -1*BOX_SIZE:
			continue
		if y > BOX_SIZE or y < -1*BOX_SIZE:
			continue
		if z > BOX_SIZE or z < -1*BOX_SIZE:
			continue
		print "#particle number", p

		i = find_index(x)
		j = find_index(y)
		k = find_index(z)
		i1 = i+1
		j1 = j+1
		k1 = k+1

		fx = (x+BOX_SIZE)*NUM_PIXELS/(2*BOX_SIZE) - i
		fy = (y+BOX_SIZE)*NUM_PIXELS/(2*BOX_SIZE) - j
		fz = (z+BOX_SIZE)*NUM_PIXELS/(2*BOX_SIZE) - k
		tx = 1. - fx
		ty = 1. - fy
		tz = 1. - fz
		
		assert fx >= 0 and tx >= 0, "f or t not between 0 and 1"
		assert fy >= 0 and ty >= 0, "f or t not between 0 and 1"
		assert fz >= 0 and tz >= 0, "f or t not between 0 and 1"

		part_prop = part_property[p]
		if i >= 0 and i < NUM_PIXELS and j >= 0 and j < NUM_PIXELS and k >= 0 and k < NUM_PIXELS:
			data[i][j][k] += part_prop * fx * fy * fz * w
			weight[i][j][k] += w
		if i >= 0 and i < NUM_PIXELS and j >= 0 and j < NUM_PIXELS and k1 >= 0 and k1 < NUM_PIXELS:
			data[i][j][k1] += part_prop * fx * fy * tz * w
			weight[i][j][k1] += w
		if i >= 0 and i < NUM_PIXELS and j1 >= 0 and j1 < NUM_PIXELS and k >= 0 and k < NUM_PIXELS:
			data[i][j1][k] += part_prop * fx * ty * fz * w
			weight[i][j1][k] += w
		if i >= 0 and i < NUM_PIXELS and j1 >= 0 and j1 < NUM_PIXELS and k1 >= 0 and k1 < NUM_PIXELS:
			data[i][j1][k1] += part_prop * fx * ty * tz * w
			weight[i][j1][k1] += w
		if i1 >= 0 and i1 < NUM_PIXELS and j >= 0 and j < NUM_PIXELS and k >= 0 and k < NUM_PIXELS:
			data[i1][j][k] += part_prop * tx * fy * fz * w
			weight[i1][j][k] += w
		if i1 >= 0 and i1 < NUM_PIXELS and j >= 0 and j < NUM_PIXELS and k1 >= 0 and k1 < NUM_PIXELS:
			data[i1][j][k1] += part_prop * tx * fy * tz * w
			weight[i1][j][k1] += w
		if i1 >= 0 and i1 < NUM_PIXELS and j1 >= 0 and j1 < NUM_PIXELS and k >= 0 and k < NUM_PIXELS:
			data[i1][j1][k] += part_prop * tx * ty * fz * w
			weight[i1][j1][k] += w
		if i1 >= 0 and i1 < NUM_PIXELS and j1 >= 0 and j1 < NUM_PIXELS and k1 >= 0 and k1 < NUM_PIXELS:
			data[i1][j1][k1] += part_prop * tx * ty * tz * w
			weight[i1][j1][k1] += w
	if PROPERTY != 'density' and PROPERTY != 'mass':
		data /= weight
	print data.shape
	print data[128][128]
	g = h5py.File("/home/fas/nagai/lm643/halo_database/Summer17/MockXray/STEP1/"+DATA_FILE, "a")
	g.create_dataset(PROPERTY, data=data)
	print g.items()

#properties = ['density', 'temperature', 'metallicity', 'x-velocity', 'y-velocity', 'z-velocity', 'mass']
main(OPTION)
exit(1)

