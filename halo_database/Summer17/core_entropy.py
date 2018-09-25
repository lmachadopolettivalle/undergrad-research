# Inputs
DIR = '/home/fas/nagai/etl28/group_scratch/Tremmel/Romulus'
#DIR = DIR+'C/'
#SimDir = DIR+'h1.cosmo50/h1.cosmo50PLK.1536gst1bwK1BH.004096'
DIR = DIR+'25/'
SimDir = DIR+'cosmo25/cosmo25p.768sg1bwK1BHe75.00'

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import ascii
import os
import math
from scipy.interpolate import interp1d
import pynbody
import tangos as db
from sys import argv
from sys import exit


core_entropy = {}

def halonames():
	for filename in os.listdir(DIR+'cosmo25/'):
		if len(filename) == 31:
			HALOS.append(filename[27:])
	print HALOS
	print "done with filenames"

#Constants
k = 8.6173303 * 1e-8 #Boltzmann in keV K**-1
Zsun = 0.0134
mu = 8.0 / 7.0
m_p = 8.4089382 * 1e-58 #electron mass in solar masses
kpc3 = 2.93799895 * 1e64 #kpc**3 in cm**3
CUT = int(0.1 / k)


def crit_dens(z=0):
    h = 0.678
    omegaL = 0.6914
    omegaM = 0.3086
    rhocrit0 = 127.131573017 # from Wiki, in Msun kpc**-3

    Ea = omegaM*math.pow((1+z), 3.0) + omegaL

    rhocrit = Ea * rhocrit0

    return rhocrit

#find R500 for given haloid
def R500(H, z):
	F = 500
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
	#rbins = H.calculate_abcissa_values('tot_mass_profile')

	j = 0
	N = len(rbins)
	crit = crit_dens(z)
	while j < N and (3 * tot_mass[j] / (4*math.pi*math.pow(rbins[j], 3.0))) > F * crit:
		j += 1

	#returned: R500
	return rbins[j]
    #returned:  R500    Mgas        M500              fgas
    #return [rbins[j], gas_mass[j], tot_mass[j], gas_mass[j]/(f_b * tot_mass[j])]



#haloid=1 is main halo
def calculate_coreK(f=0):
	haloid = int(argv[1])
	HALO = argv[2]
	count = int(argv[3])

	#data = ascii.read('./CoreEntropy'+str(haloid)+'.txt')
	#if HALO in data["id"]:
	#	exit(1)
	

	Sim = db.get_simulation('cosmo25')
	steps = Sim.timesteps
	steps = steps[::-1]
	steps = steps[:115]


	dir = SimDir+HALO
	step = steps[count]
	try:
		s = pynbody.load(dir)
	except IOError:
		if HALO == '1270':
			exit(0)
		exit(2)

	z = step.redshift
	if z < 0:
		z = 0.0

	s.physical_units()
	H = step.halos[haloid - 1]
	
	#find rad for entropy core
	if f <= 0:
		rad = 10
	else:
		rad = f * R500(H, z)	
	
	Center = H['SSC']
	Center = (Center[0], Center[1], Center[2])
	

	halo = s.gas[(pynbody.filt.Sphere(str(2*rad)+' kpc', cen=Center)) & (pynbody.filt.HighPass('temp', str(CUT)+' K'))]
	
	halo['pos'] -= Center

	#plog = pynbody.analysis.profile.Profile(halo.gas,min=1,max=1000*rad,ndim=3,type='log')
	p = pynbody.analysis.profile.Profile(halo.gas,min=0,max=2*rad,ndim=3)


	#masses = p['mass']
	rbins = p['rbins']
	dens = p['density'] * (1.0/mu) * (1.0/m_p) * (1.0/kpc3)
	temp = p['temp'] * k
	#met = p['metals'] * (1.0/Zsun)
	#clump = plog['clump']
	N = len(rbins)
	entr = [temp[i]/ (math.pow(dens[i], 2.0/3.0)) for i in range(N)]

	for j in range(N):
		if rbins[j] > rad:
			break

	x = [rbins[j-1], rbins[j]]
	y = [entr[j-1], entr[j]]
	f = interp1d(x, y)
	print HALO+" "+str(z)+" "+str(f(rad))

	if HALO == '1270':
		exit(0)

	exit(1)
	return 0

z = {8192: 0, 7936: 0.030527243817950023, 3072: 1.1853049846852821, 3478: 0.9998100146757658, 4096: 0.774398375005892, 1543: 2.509348147588138, 3336: 1.0606344086283346, 2304: 1.6712331299362209, 5271: 0.4630143029385907, 2690: 1.3992051629582094, 4352: 0.6958070914669776, 4111: 0.769587321565466, 1536: 2.5201195050422855, 4608: 0.6239681449700463, 5529: 0.40830521105884254, 2840: 1.309758252352847, 6937: 0.16435041362805203, 4864: 0.5579284854772784, 5022: 0.5197198785548629, 1826: 2.1315048812575337, 5795: 0.3557019492001796, 1270: 2.999730532662382, 5120: 0.4969122099805612, 3328: 1.0641819980905391, 7212: 0.1249144540624576, 4781: 0.5787550553108292, 5376: 0.44028086756423757, 7680: 0.06246651955101412, 1945: 1.999866242206278, 2560: 1.4834355967903856, 1458: 2.6458287033774055, 2547: 1.492236174108259, 6069: 0.30508457947303547, 5632: 0.3875040573920594, 2998: 1.2233366039282516, 1792: 2.171713042305934, 7869: 0.03874529195980081, 5888: 0.338137317038685, 3905: 0.8381287532902209, 1302: 2.9336081107837964, 1931: 2.0146662275081915, 4549: 0.6399790373658016, 6144: 0.29180526368008475, 7241: 0.12088008658345895, 2159: 1.7932867461336923, 4173: 0.7499841583835849, 6400: 0.24818857919248094, 1378: 2.78673944283867, 3584: 0.9568846984293302, 6350: 0.25650758630744686, 6656: 0.20701384803795775, 2816: 1.3235575219650384, 1550: 2.498656405479553, 3163: 1.140479124986304, 6912: 0.1680455397893872, 1632: 2.379029013172149, 2048: 1.8960673164146993, 7394: 0.09996738303039998, 7779: 0.04994019580335607, 4326: 0.7034631603233255, 7168: 0.13107962307367194, 1280: 2.978776227441135, 2536: 1.4997393169017066, 2281: 1.689754445252182, 2411: 1.588839351771342, 7424: 0.09593843407808444, 3517: 0.9837816305553304, 6640: 0.20952035371375488, 7552: 0.07900276297443654, 5107: 0.4999001078658487, 1726: 2.253427491427582, 6390: 0.24984487724165594, 2042: 1.9018800717153437, 3707: 0.9094935825547465, 3840: 0.8609266989851763}

def plot_coreK(option):
	data = {}
	for i in range(1, 6):
		data[i] = ascii.read('./Core'+option+'_'+str(i)+'.txt')
		if option == '10':
			data[i]["redshift"] = [0.0]*len(data[i])
			for j in range(len(data[i])):
				data[i]["redshift"][j] = z[int(data[i]["id"][j])]
	
	
	for i in range(1, 6):
		plt.plot(data[i]["redshift"], data[i]["CoreK"], label=str(i))
	plt.xlabel("Redshift")
	if option == '10':
		plt.ylabel(r"Core Entropy, $K_{10kpc}$")
	else:
		plt.ylabel(r"Core Entropy, $K_{0."+option+"*R500}$")
	plt.title("Core Entropy vs z - Romulus25 Groups")
	plt.yscale('log')
	plt.legend(title="Rom25 Group")
	plt.savefig('./entropy'+option+'.png')
	return 0

calculate_coreK(f=float(argv[4]))
#plot_coreK("03")
