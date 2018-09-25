#!/bin/bash
#SBATCH --job-name=CUBEC
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem=100000
#SBATCH --output="./CUBEC"
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luisfernando.machadopolettivalle@yale.edu

source /home/fas/nagai/lm643/yt-conda/bin/activate
source ~/halo_database/environment.sh

ARRAY=("density" "temperature" "metallicity" "x-velocity" "y-velocity" "z-velocity" "mass")
haloid='C'
count=6
core='NO'
#python /home/fas/nagai/lm643/halo_database/Summer17/MockXray/cube.py $haloid ${ARRAY[$count]} $core
until python /home/fas/nagai/lm643/halo_database/Summer17/MockXray/cube.py $haloid ${ARRAY[$count]} $core; do
	#ARRAY+=($?)
	#echo ${ARRAY[$count]}
	#echo ${ARRAY[@]}
	count=$[count+1]
	sleep 1
done

