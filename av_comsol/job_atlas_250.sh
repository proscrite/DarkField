## PBS OPTIONS
#PBS -N garfieldpp
#PBS -q parallel
#PBS -l nodes=1:ppn=1
#PBS -l mem=192gb
#PBS -l cput=50:00:00
#PBS -M herrero@dipc.org
#PBS -m bae

module load ROOT/6.14.06-foss-2018b-Python-2.7.15
cd av_comsol
./av_comsol > results/250_out.dat
