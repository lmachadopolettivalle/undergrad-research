#!/bin/bash
#SBATCH --job-name=300ksNOCORE5
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem=160000
#SBATCH --output="./300ksNOCORE5"
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luisfernando.machadopolettivalle@yale.edu

source /home/fas/nagai/lm643/yt-conda/bin/activate

haloid='5'
core='NO'
python /home/fas/nagai/lm643/halo_database/Summer17/MockXray/nocore_photonslist.py $haloid $core
