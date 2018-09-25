#!/usr/bin/env python

import sys
import numpy as np
import scipy
from scipy import loadtxt
from xspec import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec 
from pylab import *

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times'],'size':12})
rc('text', usetex=True)
rc('axes',labelsize='larger')
rc('xtick',labelsize='large')
rc('ytick',labelsize='large')
rc('xtick.major',size=10)
rc('xtick.minor',size=6)
rc('ytick.major',size=10)
rc('ytick.minor',size=6)
rc('axes.formatter',limits="-2, 3" )

import os.path

clight = 3.0e5
Xset.abund = "angr"
Xset.cosmo = "67.769 0.0 0.6914"
Xset.xsect = "bcmc"

Xset.addModelString("APECROOT","/home/fas/nagai/etl28/programs/Xrays/atomdb/atomdb_v3.0.8/apec")
Xset.addModelString("APECTHERMAL","yes")
Xset.addModelString("APEC_TRACE_ABUND","0.3")
#Xset.addModelString("APEC_TRACE_ABUND","Fe")
#Xset.addModelString("APEC_TRACE_ABUND","1.0")

#Fit.statMethod = "chi"
Fit.statMethod = "cstat"
Fit.query = "no"
Fit.method = "leven  100 .05"
Fit.nIterations = 10000

fig = figure( figsize=(7,7) )
grid = gridspec.GridSpec(3, 3 , hspace = 0.13)
grid.update(left=0.15, right=0.95, bottom=0.15, top=0.95)

DIR = "/home/fas/nagai/lm643/yt-conda/lib/python2.7/site-packages/pyxsim/response_files/"
rmffile = DIR+"aciss_aimpt_cy18.rmf"
arffile = DIR+"aciss_aimpt_cy18.arf"



def weighted_avg_and_std (values, weights) :

    average = np.average(values, weights = weights)
    variance = np.average( (values-average)**2, weights = weights )

    return (average, math.sqrt(variance))

def sigma_to_energy_width ( sigma, E ):

    #DeltaE/E = v/c
    c = 3.0e5 #speed of light in km/s
    return sigma*E/c

#TODO -- DONE
paramfile = open("./tinyfit.param", "w")
print >>paramfile, "first line"
phafile = '100kpc_TINY.pha'

if ( os.path.isfile(phafile)): 
	#to make this next line work, i needed to do the "ln -s" trick to
	#create a symbolic link to the response files
	s = Spectrum(phafile)
	s.response = rmffile
	s.response.arf = arffile 
	#TODO -- DONE
	#s.ignore("**-5. 10.-**")
	s.ignore("**-.1 5.-**")
	#s.ignore("**-.1 10.-**")

	#check pyxspec documentation to see which model!
	#TODO
	m = Model("wabs(apec)")
	#(wabs factor;
	#plasma temperature keV;
	#metal abundances;
	#redshift;
	#normalization)
	m.setPars(0.138,11.0,0.30,0.0157,0.01)
	#m.setPars(0.138,2.42,0.30,0.0157,0.01)
	#ORIGINAL MODEL FROM ERWIN:
	#m = Model("wabs(bapec)")
	#m.setPars(0.138,11.0,0.30,0.017284,80.0,0.01)

	for i in range(m.nParameters):
		m(i+1).frozen = False
	#m(1).frozen = True
	#m(2).frozen = True
	#m(3).frozen = True
	#m(4).frozen = True
	pfrozenlist = [ m(i+1).name for i in range(m.nParameters) if m(i+1).frozen == True ]

	Fit.perform()

	dof = float(len(s.noticed))-float(m.nParameters)+float(len(pfrozenlist))
	redchisq = Fit.statistic/dof

	#TODO -- DONE
	Tx = m(2).values[0]
	Abund = m(3).values[0]
	Redshift = m(4).values[0]
	#Velocity = m(5).values[0]
	Norm = m(5).values[0]	

	#TODO -- DONE
	TxErr = m(2).sigma	
	AbundErr = m(3).sigma	
	RedshiftErr = m(4).sigma	
	#VelocityErr = m(5).sigma	
	NormErr = m(5).sigma	

	xm = [ 0.5*(s.energies[i][0] + s.energies[i][1]) for i in range(len(s.energies)) ]
	delx = [ s.energies[i][1] - s.energies[i][0] for i in range(len(s.energies)) ]
	xerr = delx 
	delx = 1.0

	xd = xm
	yd = array(s.values)/delx
	ym = array(m.folded(1))/delx
	res = yd - ym
	ratio = yd/ym
	yerr = sqrt(array(s.variance))/delx
	ratioerr = (yerr)/ym 

	ax1 = plt.subplot(grid[:-1,:])
	ax2 = plt.subplot(grid[-1,:])
	ax1.errorbar(xd, yd, xerr=xerr, yerr=yerr, fmt='none', capsize=0, color='b',label="data")
	ax1.plot(xm, ym, 'r',label="fit")
	ax2.errorbar(xd, ratio, xerr=xerr, yerr=ratioerr, fmt='none', capsize=0, color='b')
	ax2.axhline(1.0,color='k',ls='--')

	plt.setp(ax1.get_xticklabels(), visible=False)	
	ax2.set_xlabel(r'Energy [keV]')
	ax1.set_ylabel(r'Normalized counts sec$^{-1}$ keV$^{-1}$')
	ax2.set_ylabel(r'Ratio = Data/Model')

	ax1.yaxis.labelpad = 20
	ax2.yaxis.labelpad = 10
	ax2.xaxis.labelpad = 20

	#velstr = "{0:.2f}".format(Velocity)
	#velerrstr = "{0:.2f}".format(VelocityErr)
	shiftstr = "{0:.2f}".format((Redshift-0.017284)*3e5)
	shifterrstr = "{0:.2f}".format(RedshiftErr*3e5)
	chisqstr = "{0:.3f}".format(redchisq)

	#TODO
	#change these texts to be appropriate to my spectra

	#ax1.text(0.62,0.80,r"CL "+str(id)+" "+proj,transform=ax1.transAxes, horizontalalignment='left') 
	#ax1.text(0.62,0.74,r"$\sigma ="+velstr+r"\pm"+velerrstr+r"$ km/s",transform=ax1.transAxes, horizontalalignment='left')
	#ax1.text(0.62,0.68,r"$v ="+shiftstr+r"\pm"+shifterrstr+r"$ km/s",transform=ax1.transAxes, horizontalalignment='left')

	tempstr = "{0:.2f}".format(Tx)
	temperrstr = "{0:.2f}".format(TxErr)
	ax1.text(0.62,0.62,r"$T_x ="+tempstr+r"\pm"+temperrstr+r"$ keV",transform=ax1.transAxes, horizontalalignment='left')

	#ax1.set_xlim(6.4,8)
	#ax2.set_xlim(6.4,8)
#ax1.set_xlim(6.78,6.88)
	#ax2.set_xlim(6.78,6.88)

	fig.add_subplot(ax1)
	fig.add_subplot(ax2,sharex=ax1)

	fig.savefig("secondattempt.png",dpi=300)
	#fig.savefig("firstattempt.png",dpi=300)
	#fig.savefig("CL"+str(id)+"_"+proj+"_bapec.png",dpi=300)
	fig.clf()

	bulk = (Redshift-0.017284)*3e5
	bulkErr = (RedshiftErr)*3e5
	#TODO -- DONE
	print >>paramfile,Tx,TxErr,Abund,AbundErr,bulk,bulkErr,Norm,NormErr,redchisq

	AllModels.clear()
	AllData.clear()


else :
	print("File "+phafile+" doesn't exist!")

paramfile.close()
	
exit()

