#################################
# 								#
# nocore_photonlist.py			#
# 								#
# Read .hdf5 file and load 		#
# into yt, then make photon		#
# list and write .hdf5 file		#
#								#
#################################

import yt
#yt.enable_parallelism()
import pyxsim
import h5py
import numpy as np
from yt.utilities.cosmology import Cosmology
from yt.utilities.physical_constants import mh
from yt.units import Msun, parsec, g, s, km
from yt.units.yt_array import YTArray
from sys import argv, exit

#program, haloid, count, originalcount = argv
#count = int(count)
#originalcount = int(originalcount)
#if count - originalcount >= 4:
	#exit(0)
program, haloid, core = argv


R500 = {'C': 672.5, '1': 375.15, '2': 359.05, '3': 332.35, '4': 325.35, '5': 304.45}
R200 = {'C': 1200.0, '1': 551.55, '2': 549.05, '3': 488.35, '4': 481.85, '5': 451.75}

rad = 2*R200[haloid]
L = rad / 2
width = str(2*L)+' kpc'

BOX_SIZE = L #box goes from -1*BOX_SIZE to +1*BOX_SIZE; in kpc
NUM_PIXELS = 512 #number of pixels from -L to +L

PARTICLE_FILE = str(rad)+'kpc_Rom'+haloid+'Particles.hdf5'
DATA_FILE = str(rad)+'kpc_Rom'+haloid+'Data.hdf5'
PHOTON_FILE = 'Rom'+haloid+'_Photon.hdf5'
if core == 'NO':
	PARTICLE_FILE = 'NOCORE_' + PARTICLE_FILE
	DATA_FILE = 'NOCORE_' + DATA_FILE
	PHOTON_FILE = 'NOCORE_' + PHOTON_FILE

X_H = 0.76


f = h5py.File("/home/fas/nagai/lm643/halo_database/Summer17/MockXray/STEP1/"+DATA_FILE, "a")

def emission_measure():
	dens = YTArray(f['density'].value, "Msun/kpc**3").in_units('g/cm**3')
	mass = YTArray(f['mass'].value, "Msun").in_units('g')

	nenh = dens/mh
	nenh *= nenh
	cell_v = mass/dens
	nenh *= 0.5*(1.+X_H)*X_H*cell_v

	f.create_dataset('emission_measure',data=nenh)
	exit()

emission_measure()



units = ["Msun/kpc**3", "cm**(-3)", "Msun", "1", "K", "km/s", "km/s", "km/s"]
data = {k:YTArray(v.value,u) for (k,v), u in zip(f.items(),units)}
data['metallicity'] = f['metallicity'].value

bbox = np.array([[-L, L], [-L, L], [-L, L]])
ds = yt.load_uniform_grid(data, data["density"].shape, length_unit="kpc", bbox=bbox)

ad = ds.all_data()
sp = ds.sphere([0, 0, 0],(L, "kpc"))
#sp = ds.sphere("c",(L, "kpc"))
center = sp.center

#yt.ProjectionPlot(ds, "z", "density", width=(L, "kpc")).save('./projectionOct29.png')
#yt.SlicePlot(ds, 'x', 'density').save('./sliceplotOct29.png')
#q = yt.ProfilePlot(sp, "radius", "density")
#q.set_unit('radius', 'kpc')
#q.save('./profileOct29.png')
#exit()


cosmo = Cosmology(hubble_constant=0.67769, omega_matter=0.3086, omega_lambda=0.6914)

redshift = 0.0157
#ULTIMATE LARGE VALUES: 5000 cm2, 1000 ks
exp_time = (300., "ks")
area = (1000., "cm**2")


source_model = pyxsim.ThermalSourceModel("mekal", 0.1, 10.0, 300, Zmet="metallicity", temperature_field=('stream', 'temperature'),emission_measure_field=('stream', 'emission_measure'))
#source_model = pyxsim.ThermalSourceModel("apec", 0.1, 10.0, 300, Zmet="metallicity", temperature_field=('stream', 'temperature'),emission_measure_field=('stream', 'emission_measure'))

#photons = pyxsim.PhotonList.from_data_source(sp, redshift, area, exp_time, source_model, center=center, cosmology=cosmo,velocity_fields=["x-velocity", "y-velocity", "z-velocity"])
photons = pyxsim.PhotonList.from_data_source(ad, redshift, area, exp_time, source_model, cosmology=cosmo,velocity_fields=["x-velocity", "y-velocity", "z-velocity"])

photons.write_h5_file('/home/fas/nagai/lm643/scratch60/STEP2/'+str(exp_time[0])+'ks_'+PHOTON_FILE)

#exit(1)
exit(0)

