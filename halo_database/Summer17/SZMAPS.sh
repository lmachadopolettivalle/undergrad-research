#!/bin/bash
#SBATCH --job-name=SZMapsRomC
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem=100000
#SBATCH --output="./SZMapsRomC"
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luisfernando.machadopolettivalle@yale.edu

source /home/fas/nagai/lm643/yt-conda/bin/activate
source ~/halo_database/environment.sh

count=1
ARRAY=("TSZ" "kSZ")
until python /home/fas/nagai/lm643/halo_database/Summer17/SZmaps.py ${ARRAY[$count]}; do
	#ARRAY+=($?)
	#echo ${ARRAY[$count]}
	#echo ${ARRAY[@]}
	count=$[count+1]
	sleep 1
done

