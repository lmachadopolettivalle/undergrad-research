# afile has all data about temperature estimates
# bfile comes from some recent run of Mike.sh
# newfile will be the new file containing afile + bfile data

from astropy.io import ascii
from sys import argv, exit

if len(argv) != 4:
	print "\nUsage: ./writeTemperatures.py ALLDATA SOMEDATA NEWNAME"
	print "Write file names WITH .txt ending\n"
	exit()

program, a, b, newfile = argv

A = open(a, 'r')
B = open(b, 'r')
aline0 = A.readline().strip()
bline0 = B.readline().strip()

field = bline0.split()[-1]

new = open(newfile, 'w')
new.write(aline0+' '+field+'\n')

alines = A.readlines()
blines = B.readlines()

N = len(alines)
for i in xrange(N):
	new.write(alines[i].strip()+' '+blines[i].strip().split()[-1]+'\n')
new.close()
A.close()
B.close()


