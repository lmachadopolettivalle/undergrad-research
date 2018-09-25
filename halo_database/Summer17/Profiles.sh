#!/bin/bash
#SBATCH --job-name=RomCProfile
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem=100000
#SBATCH --output="./RomCProfile"
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luisfernando.machadopolettivalle@yale.edu

source /home/fas/nagai/lm643/yt-conda/bin/activate
source ~/halo_database/environment.sh

python /home/fas/nagai/lm643/halo_database/Summer17/profiles.py
#mpirun -n 4 /home/fas/nagai/lm643/halo_database/Summer17/calc_profiles.py
#mpirun /home/fas/nagai/lm643/halo_database/Summer17/calc_profiles.py
#python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 dens
#python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 temp
#python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 tempmass
#python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 entr
#python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 entrmass
#python /home/fas/nagai/lm643/halo_database/Summer17/maps.py 1 met
