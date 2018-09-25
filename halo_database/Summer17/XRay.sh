#!/bin/bash
#SBATCH --job-name=XRayMap
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem=120000
#SBATCH --output="./XRayMap"
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luisfernando.machadopolettivalle@yale.edu

source /home/fas/nagai/lm643/yt-conda/bin/activate
source ~/halo_database/environment.sh


python /home/fas/nagai/lm643/halo_database/Summer17/maps.py "XRay"

