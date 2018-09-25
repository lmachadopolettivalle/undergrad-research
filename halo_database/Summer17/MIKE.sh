#!/bin/bash
#SBATCH --job-name=Vik_halfcut
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem=100000
#SBATCH --output="./Vik_halfcut"
#SBATCH --time=1:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luisfernando.machadopolettivalle@yale.edu

source /home/fas/nagai/lm643/yt-conda/bin/activate
source ~/halo_database/environment.sh

haloid=25
#CUT=1000 - on Oct16, 2017 Nagai chose to stick with CUT=1000
#hence, this is not a needed input anymore
Function='Vik'
#Function='Luis'
#Function='Mike'
#half='False'
half='True'
until python /home/fas/nagai/lm643/halo_database/Summer17/compare_temperatures.py $haloid $Function $half; do
	haloid=$[haloid+1]
	sleep 1
done
