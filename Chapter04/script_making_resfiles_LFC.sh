#this script makes resfiles from standard Rosetta named Match outputs. 
for i in *.pdb;
do
echo "NATRO" >> $i.resfile
echo "start" >> $i.resfile
pos_array=(16 18 19 20 21 22 23 24 25 26 27 28 30 32 33 34 36 37 40 73 74 76 77 78 79 80 81 82 83 84 85 86 87 88 90 91 92 93 94 95 96 97 98)

variable=$(echo $i |awk -F'[__]' '{print $3}')
echo $variable
variable_temp=$(echo $variable | sed -e 's/[A-Z]\(.*\)[A-Z]/\1/')
echo $variable_temp
first_position=$(echo "${variable_temp}" | cut -c2-3)
echo $first_position
#second_position=$(echo "${variable_temp}" | cut -c3-4)


echo $first_position "A NATAA" >> $i.resfile
#echo $second_position "A NATAA" >> $i.resfile
delete=($first_position)
#delete1=($second_position)
pos_array_1=$(echo ${pos_array[@]/$delete})
#pos_array_2=$(echo ${pos_array_1[@]/$delete1})

echo $pos_array_1

for t in ${pos_array_1[@]}; do
  #echo $t
  echo $t "A AUTO NOTAA CH " >> $i.resfile
done
done
