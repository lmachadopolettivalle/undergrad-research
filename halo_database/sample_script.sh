#!/bin/bash
## job name
#PBS -N SampleJob
## queue type
#PBS -q astro_prod
## total cpu resource allocation
#PBS -l nodes=1:ppn=8,mem=35gb
## run time allocation
#PBS -l walltime=3:00:00
## combine stderr &120 stdout into one file
#PBS -j oe
## where to put run time information
#PBS -o ./outputs/$PBS_JOBNAME
## email results: a = aborted, b = begin, e = end
#PBS -m abe
#PBS -M luisfernando.machadopolettivalle@yale.edu
#PBS -V
cd $PBS_O_WORKDIR

source /home/fas/nagai/lm643/yt-x86_64/bin/activate
source ~/halo_database/environment.sh

mpirun -n 4 /home/fas/nagai/lm643/halo_database/Summer17/calc_profiles.py
#mpirun /home/fas/nagai/lm643/halo_database/Summer17/calc_profiles.py
python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 dens
python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 temp
python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 tempmass
python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 entr
python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 entrmass
python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 met
