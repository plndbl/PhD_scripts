#!/bin/bash

#SBATCH --partition=veryshort
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=6:00:00
#SBATCH --mem=10000M


module load apps/rosetta/2018.33
module load languages/anaconda3/3.6.5
module load languages/anaconda2/5.0.1

export cluster_Rosetta=/mnt/storage/software/apps/Rosetta-2018.33/rosetta_bin_linux_2018.33.60351_bundle
export AB_Rosetta=/mnt/storage/scratch/hb18661/rosetta_bin_linux_2019.35.60890_bundle

export PATH=${PATH}:/mnt/storage/software/apps/Rosetta-2018.33/rosetta_bin_linux_2018.33.60351_bundle/tools/protein_tools/scripts

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

for i in *pdb;
do
echo "NATRO" >> EOY4D2.resfile
echo "start" >> EOY4D2.resfile
pos_array=(16 18 19 20 21 22 23 24 25 26 27 28 30 32 33 34 36 37 40 73 74 76 77 78 79 80 81 82 83 84 85 86 87 88 90 91 92 93 94 95 96 97 98)

variable=$(echo $i |awk -F'[__]' '{print $3}')
variable_temp=$(echo $variable | sed -e 's/[A-Z]\(.*\)[A-Z]/\1/')
first_position=$(echo "${variable_temp}" | cut -c1-2)
second_position=$(echo "${variable_temp}" | cut -c3-4)


echo $first_position "A NATAA" >> EOY4D2.resfile
echo $second_position "A NATAA" >> EOY4D2.resfile
delete=($first_position)
delete1=($second_position)
pos_array_1=$(echo ${pos_array[@]/$delete})
pos_array_2=$(echo ${pos_array_1[@]/$delete1})

echo $pos_array_1

for t in ${pos_array_2[@]}; do
  echo $t
  echo $t "A AUTO NOTAA CH " >> EOY4D2.resfile
done

srun -u rosetta_scripts.static.linuxgccrelease -parser:protocol Design.xml @design.flags -s $i

rm EOY4D2.resfile
mv *output_design.pdb output_design
done 



