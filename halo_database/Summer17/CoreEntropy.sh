#!/bin/bash
## job name
#PBS -N Core01_1
## queue type
#PBS -q astro_prod
## total cpu resource allocation
#PBS -l nodes=2:ppn=8,mem=70gb
## run time allocation
#PBS -l walltime=6:00:00
## combine stderr and stdout into one file
#PBS -j oe
## where to put run time information
#PBS -o /home/fas/nagai/lm643/halo_database/outputs/$PBS_JOBNAME.txt
## email results: a = aborted, b = begin, e = end
#PBS -m abe
#PBS -M luisfernando.machadopolettivalle@yale.edu
#PBS -V
cd $PBS_O_WORKDIR

source /home/fas/nagai/lm643/yt-x86_64/bin/activate
source ~/halo_database/environment.sh

ARRAY=(8192 7936 7869 7779 7680 7552 7424 7394 7241 7212 7168 6937 6912 6656 6640 6400 6390 6350 6144 6069 5888 5795 5632 5529 5376 5271 5120 5107 5022 4864 4781 4608 4549 4352 4326 4173 4111 4096 3905 3840 3707 3584 3517 3478 3336 3328 3163 3072 2998 2840 2816 2690 2560 2547 2536 2411 2304 2281 2159 2048 2042 1945 1931 1826 1792 1726 1632 1550 1543 1536 1458 1378 1302 1280 1270)
count=0
until python /home/fas/nagai/lm643/halo_database/Summer17/core_entropy.py 1 ${ARRAY[$count]} $count 0.1; do
	#ARRAY+=($?)
    #echo ${ARRAY[$count]}
	#echo ${ARRAY[@]}
    count=$[count+1]
    sleep 1
done

