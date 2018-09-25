from astropy.io import ascii
from sys import argv
import math

dens = ascii.read('./NGCdensity.txt')
temp = ascii.read('./NGCtemp.txt')
Dens = [float(k) for k in dens["density"]]
EPD = [float(k) for k in dens["eplus"]]
EMD = [float(k) for k in dens["eminus"]]
Temp = [float(k) for k in temp["temp"]]
EPT = [float(k) for k in temp["eplus"]]
EMT = [float(k) for k in temp["eminus"]]

f = open('./NGCpressure.txt', "w+")
f.write('rbins pressure eplus eminus\n')
for j in range(len(temp)):
	p = Dens[j]*Temp[j]
	ep = p*math.pow((math.pow((EPD[j]/Dens[j]), 2) + math.pow((EPT[j]/Temp[j]), 2)), 0.5)
	em = p*math.pow((math.pow((EMD[j]/Dens[j]), 2) + math.pow((EMT[j]/Temp[j]), 2)), 0.5)
	f.write(str(dens["rbins"][j])+' '+str(p)+' '+str(ep)+' '+str(em)+'\n')

exit()
array = ['C', '1', '2', '3', '4', '5']
for i in array:
	dens = ascii.read('./density_'+i+'.txt')
	temp = ascii.read('./temp_'+i+'.txt')
	Dens = [float(k) for k in dens["density"]]
	Temp = [float(k) for k in temp["temp"]]
	f = open('./pressure_'+i+'.txt', "w+")
	f.write('rbins pressure\n')
	for j in range(len(dens)):
		f.write(str(dens["rbins"][j])+' '+str(Dens[j]*Temp[j])+'\n')

	f.close()

