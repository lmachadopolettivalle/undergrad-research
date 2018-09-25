#################################
# 								#
# events_list.py					#
# 								#
# Read .hdf5 file and load 		#
# into yt, then make photon		#
# list and write .hdf5 file		#
#								#
#################################

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
import numpy as np
import h5py
#from yt.utilities.cosmology import Cosmology
#from yt.utilities.physical_constants import mh
#from yt.units import Msun, parsec, g, s, km
#from yt.units.yt_array import YTArray
import pyxsim
import soxs

rad = 100

z = 0.0157
area = (300., "cm**2")
time = (1., "ks")

fov = (1.2, "arcmin") # the field of view / width of the image
nx = 1024 # The resolution of the image on a side

def events():
	photons = pyxsim.PhotonList.from_file("Rom1smallphotonlist.hdf5")
	events = photons.project_photons("z", (30.0, 45.0), area_new=area, exp_time_new=time, redshift_new=z)
	events.write_fits_file(str(rad)+"kpc_SMALL_events.fits", fov, nx, overwrite=True)
	events.write_simput_file(str(rad)+"kpc_SMALL", overwrite=True, emin=0.1, emax=9.0)
	

def showfits():
	events = pyxsim.event_list.EventList.from_fits_file("100kpc_TINY_events.fits")
	#events = pyxsim.event_list.EventList.from_fits_file(str(rad)+"kpc_events.fits")

	events.write_fits_image("myimage.fits", fov, nx, overwrite=True, emin=0.5, emax=7.0)
	hdu_list = fits.open('./myimage.fits')
	image_data = hdu_list[0].data
	#image_data = fits.getdata(image_file)
	plt.imshow(image_data, cmap='jet',norm=LogNorm(vmin=0.6, vmax=2))
	plt.colorbar()
	plt.savefig('100kpc_TINY.png')
	plt.clf()
	fig = soxs.plot_spectrum(str(rad)+"kpc_TINY.pha")
	fig.axes[0].set_xscale("log")
	fig.axes[0].set_yscale("log")
	fig.axes[0].set_xlim(0.4, 9.0)
	fig.savefig('tryingspectrumtuesday.png')

def usesoxs():
	from soxs import instrument_simulator

	events = pyxsim.event_list.EventList.from_fits_file(str(rad)+"kpc_SMALL_events.fits")
	simput_file = str(rad)+"kpc_SMALL_simput.fits" # SIMPUT file to be read
	out_file = "OBSERVED.fits" # event file to be written
	exp_time = 45000. # The exposure time in seconds
	instrument = "acisi_cy18" # short name for instrument to be used
	sky_center = [30., 45.] # RA, Dec of pointing in degrees
	instrument_simulator(simput_file, out_file, exp_time, instrument, sky_center, overwrite=True) 
	soxs.write_spectrum(out_file, str(rad)+"kpc_SMALL.pha", overwrite=True)
	exit()

events()
usesoxs()
showfits()

