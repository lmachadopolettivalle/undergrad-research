#!/bin/bash
#SBATCH --job-name=Scaling
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem=100000
#SBATCH --output="./Scaling"
#SBATCH --time=0:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luisfernando.machadopolettivalle@yale.edu

source /home/fas/nagai/lm643/yt-conda/bin/activate
source ~/halo_database/environment.sh

until python /home/fas/nagai/lm643/halo_database/Summer17/scaling.py 'C'; do
	sleep 1
done
