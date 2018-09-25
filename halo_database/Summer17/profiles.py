import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sys import argv
import math
import pynbody
#import tangos as db
from astropy.io import ascii


center = {'C': (2546.87015086, -97.94119966, -2971.89876526), 1: (7446.7740934127805, 8478.3427413597001, 7965.3717622399463), 2: (5888.7597183651769, -3920.9626619405344, 9269.1704642761524), 3: (5338.6136956516802, -9639.365229590805, 7312.5757339895217), 4: (11731.934544501413, 10293.593262519084, 10325.876620841929), 5: (1976.5326767473034, -6942.4867710401186, 4635.9382550076662)}

#Constants
k = 8.6173303 * 1e-8 #Boltzmann in keV K**-1
Zsun = 0.0134
mu = 8.0 / 7.0
m_p = 8.4089382 * 1e-58 #electron mass in solar masses
kpc3 = 2.93799895 * 1e64 #kpc**3 in cm**3

def write_profiles(haloid):
	# Inputs
	rad = 1 #in Mpc
	CUT = 0.3 #only consider particles with T > CUT (keV)

	DIR = '/home/fas/nagai/lm643/project/cosmo25/'
	SimDir = DIR+'/cosmo25p.768sg1bwK1BHe75.008192'
	if haloid == 'C':
		DIR = '/home/fas/nagai/lm643/project/h1.cosmo50/'
		SimDir = DIR+'/h1.cosmo50PLK.1536gst1bwK1BH.004096'
	
	s = pynbody.load(SimDir)
	s.physical_units()
	#h = s.halos()
	#h1 = h[1]
	Center = center[haloid]


	part = 'gas'
	#part = 'dm'
	#part = 'star'
	halo = []
	if part == 'gas':
		#halo = s.gas[pynbody.filt.Sphere(str(rad)+' Mpc', cen=Center)]
		halo = s.gas[(pynbody.filt.Sphere(str(rad)+' Mpc', cen=Center)) & (pynbody.filt.HighPass('temp', str(CUT/k)+' K'))]
	elif part == 'dm':
		halo = s.dm[pynbody.filt.Sphere(str(rad)+' Mpc', cen=Center)]
	elif part == 'star':
		halo = s.s[pynbody.filt.Sphere(str(rad)+' Mpc', cen=Center)]
	halo['pos'] -= Center


	#plog = pynbody.analysis.profile.Profile(halo.gas,min=1,max=1000*rad,ndim=3,type='log')
	if part == 'gas':
		p = pynbody.analysis.profile.Profile(halo.gas,min=0,max=1000*rad,ndim=3,nbins=200)
	elif part == 'dm':
		p = pynbody.analysis.profile.Profile(halo.dm,min=0,max=1000*rad,ndim=3,nbins=200)
	elif part == 'star':
		p = pynbody.analysis.profile.Profile(halo.s,min=0,max=1000*rad,ndim=3,nbins=200)


	#masses = p['mass']
	#rlogbins = plog['rbins']
	rbins = p['rbins']
	dens = p['density'] * (1.0/mu) * (1.0/m_p) * (1.0/kpc3)
	temp = p['temp'] * k
	#met = p['metals'] * (1.0/Zsun)
	#clump = plog['clump']
	N = len(rbins)
	entr = [temp[i]/ (math.pow(dens[i], 2.0/3.0)) for i in range(N)]


	#b = open('./'+part+'mass_'+str(haloid)+'.txt', 'w+')
	#b.write('rbins\tmass\n')
	#c = open('./clump_'+str(haloid)+'.txt', 'w+')
	#c.write('rbins\tclump\n')
	f = open('./hotdensity_'+str(haloid)+'.txt', 'w+')
	f.write('rbins\tdensity\n')
	g = open('./hottemp_'+str(haloid)+'.txt', 'w+')
	g.write('rbins\ttemp\n')
	#h = open('./met_'+str(haloid)+'.txt', 'w+')
	#h.write('rbins\tmet\n')
	i = open('./hotentr_'+str(haloid)+'.txt', 'w+')
	i.write('rbins\tentr\n')


	for j in range(N):
		#b.write(str(rbins[j])+'\t'+str(masses[j])+'\n')
		#c.write(str(rlogbins[j])+'\t'+str(clump[j])+'\n')
		f.write(str(rbins[j])+'\t'+str(dens[j])+'\n')
		g.write(str(rbins[j])+'\t'+str(temp[j])+'\n')
		#h.write(str(rbins[j])+'\t'+str(met[j])+'\n')
		i.write(str(rbins[j])+'\t'+str(entr[j])+'\n')

	#b.close()
	#c.close()
	f.close()
	g.close()
	#h.close()
	i.close()
	return 0


h = 0.678
# H(z)**2 / H_0**2
def Ea(z):
	omegaL = 0.6914
	omegaM = 0.3086

	return omegaM*math.pow((1+z), 3.0) + omegaL

#from Nagai et al. 2007
#keV
def T500(M500, z=0):
	m = (M500 * h * 1e-15)
	
	return 11.05 * math.pow(m, 2.0/3.0) * math.pow(Ea(z), 2.0/3.0)

#keV cm**2
def K500(M500, z=0):
    m = (M500 * h * 1e-15)

    return 1963 * math.pow(m, 2.0/3.0) / math.pow(Ea(z), 2.0/3.0)



def plot_profiles():
	program, option = argv
	#option = str(raw_input('density, temp, met, entr, c, hottemp, hotentr, hotdensity, pressure? '))
	if option == 'c':
		option = 'clump'
	#norm = str(raw_input('Normalize by R500, T500, K500? (y / n) '))
	norm = 'y'
	
	
	NGC_R500 = 456.0
	NGC_R500_plus = 456.0 + 28.0
	NGC_R500_minus = 456.0 - 28.0
	#Following values of M500 calculated with four methods:
	#
	#1 and 2 from Sun et al. 2009, Table 6
	#M500 = M * (T500 / 3 keV)^a, with:
	#T500 ~ T_X = 1.2 keV
	#1: M = 1.27+-0.12*1e14*h**-1, a = 1.67+-0.15, "Tier 1 + 2"
	#2: M = 1.26+-0.07*1e14*h**-1, a = 1.65+-0.04, "Tier 1 + 2 + clusters"
	#
	#3 and 4 from Morandi et al. 2017, Figure C1, with:
	#3: Top Line
	#4: Bottom Line
	#
	#NGC_M500 = [4.055, 4.098, 3.196, 2.804]
	#NGC_M500 = [i * 1e13 for i in NGC_M500]
	#NGC_T500 = T500(NGC_M500)
	#NGC_K500 = K500(NGC_M500)

	NGC_M500 = (4.055 + 4.098) * 1e13 / 2
	NGC_T500 = 1.01012170948
	NGC_K500 = 179.445150742
	
	#data about Romulus25 groups
	#computed with functions above
	#data = ascii.read('./fivehundred.txt')
	data = ascii.read('./500and2500info.txt')
	R500 = data['R500']
	M500 = data['M500']
	T500 = data['T500']
	K500 = data['K500']

	#more data about Romulus25 groups
	#computed with calc_profiles.py
	data = {}
	for i in range(1,6):
		data[i] = ascii.read('./Profilestxt/'+option+'_'+str(i)+'.txt')
	
	data[0] = ascii.read('./Profilestxt/'+option+'_C.txt')

	#more data about NGC
	try:
		NGC = ascii.read('./Profilestxt/NGC'+option+'.txt')
	except IOError:
		NGC = ascii.read('./Profilestxt/NGC'+option[3:]+'.txt')
	NGCr = NGC["rbins"]
	L = len(NGCr)
	try:
		NGCy = NGC[option]
		NGCupper = [NGC["eplus"][i] - NGC[option][i] for i in range(L)]
		NGClower = [NGC[option][i] - NGC["eminus"][i] for i in range(L)]
	except:
		NGCy = NGC[option[3:]]
		NGCupper = [NGC["eplus"][i] - NGC[option[3:]][i] for i in range(L)]
		NGClower = [NGC[option[3:]][i] - NGC["eminus"][i] for i in range(L)]

	#errorbars for R/R500, left (l) and right (r)
	NGCrl = [(1.0/NGC_R500 - 1.0/NGC_R500_plus)*i for i in NGCr]
	NGCrr = [(1.0/NGC_R500_minus - 1.0/NGC_R500)*i for i in NGCr]
	
	if norm == 'y':
		NGCr = NGCr / NGC_R500
		if option == 'temp' or option == 'hottemp' or option == 'pressure':
			NGCy = [j / NGC_T500 for j in NGCy]
			NGCupper = [j / NGC_T500 for j in NGCupper]
			NGClower = [j / NGC_T500 for j in NGClower]
		elif option == 'entr' or option == 'hotentr':
			NGCy = [j / NGC_K500 for j in NGCy]
			NGCupper = [j / NGC_K500 for j in NGCupper]
			NGClower = [j / NGC_K500 for j in NGClower]


	rbins = {}
	y = {}
	for i in range(0, 6):
		rbins[i] = data[i]["rbins"]
		y[i] = data[i][option]
		if norm == 'y':
			rbins[i] = [j / R500[i] for j in rbins[i]]
			if option == 'temp' or option == 'hottemp' or option == 'pressure':
				y[i] = [j / T500[i] for j in y[i]]
			elif option == 'entr' or option == 'hotentr':
				y[i] = [j / K500[i] for j in y[i]]
	 
	 
	plt.errorbar(NGCr, NGCy, xerr=[NGCrl, NGCrr], yerr=[NGClower, NGCupper], fmt='o', label='NGC2563')
	#plt.plot(NGCr, NGCy, 'o', label='NGC2563')

	for i in range(0, 6):
		if option == 'clump':
			for j in range(len(rbins[i])):
				if rbins[i][j] < 1.5 and y[i][j] > 1.5:
						y[i][j] = 1.2
				if rbins[i][j] > 1.5 and y[i][j] > 4:
						y[i][j] = y[i][j-1]
		if i == 0:
			plt.plot(rbins[i], y[i], label='RomC')
		else:
			#plt.plot(rbins[i], y[i], 'o', label=str(i))
			plt.plot(rbins[i], y[i], label='Rom25 #'+str(i))

	if option in ['density', 'hotdensity', 'entr', 'hotentr', 'pressure']:
		plt.yscale('log')
	
	if norm == 'y':
		plt.xlabel(r'$R / R_{500}$', fontsize=20)
		if option == 'temp' or option == 'hottemp':
			if option == 'hottemp':
				plt.title(r'Hot Gas ($>0.3 keV$) Temperature', fontsize=20)
			else:
				plt.title(r'Gas Temperature', fontsize=20)
			plt.ylabel(r'$T / T_{500}$', fontsize=20)
			plt.ylim([0, 2.5])
		elif option == 'pressure':
			plt.title('Electron Pressure', fontsize=20)
			plt.ylabel(r'$n_e \cdot T / T_{500} (cm^{-3})$', fontsize=20)
		elif option == 'entr' or option == 'hotentr':
			if option == 'hotentr':
				plt.title(r'Hot Gas ($>0.3 keV$) Entropy', fontsize=20)
			else:
				plt.title(r'Gas Entropy', fontsize=20)
			plt.ylabel(r'$K / K_{500}$', fontsize=20)
			plt.ylim([0.06, 10])
	else:
		plt.xlabel(r'$R (kpc)$', fontsize=20)
		if option == 'temp' or option == 'hottemp':
			plt.ylabel(r'$T (keV)$', fontsize=20)
		elif option == 'entr' or option == 'hotentr':
			plt.ylabel(r'$K (keV cm^2)$', fontsize=20)
	
	if option == 'density':
		plt.ylim([1e-6, 1e-1])
		plt.ylabel(r'$n_e (cm^{-3})$', fontsize=20)
		plt.title('Electron Number Density', fontsize=20)
	elif option == 'hotdensity':
		plt.ylim([1e-6, 1e-1])
		plt.ylabel(r'$n_e (cm^{-3})$', fontsize=20)
		plt.title('Hot Gas ($>0.3 keV$) Electron Number Density', fontsize=20)
	elif option == 'met':
		plt.title('Gas Metallicity', fontsize=20)
		plt.ylabel(r'$Z (Z_{Sun})$', fontsize=20)
		plt.ylim([0,1])
	elif option == 'clump':
		plt.title(r'Clumping Factor', fontsize=20)
		plt.ylabel('C', fontsize=20)
		if norm == 'y':
			plt.xlim([0.2,3.0])
		else:
			plt.xscale('log')
			plt.xlim([10,1000])
		plt.ylim([0.9,5])

	plt.legend()
	plt.savefig('./FINAL_'+option+'_profile.png')

	return 0

#plot_profiles()
write_profiles(1)
