#!/bin/bash
echo "0=$0"
echo "1=$1"
echo "2=$2"
echo "3=$3"
echo "4=$4"
echo "5=$5"

user=$3
for FILE in $1
do
	waiting=true
	while "$waiting" -eq "true"
	do
		squeue -u $user -o %j_%i_%n | egrep 'kallisto|STAR'  &> currentJobs.txt 
		nodeJobs=$(grep "$nodeName" currentJobs.txt | grep "kallisto" | wc -l)
		if [[ "$nodeJobs" -lt 200 ]] # queue a maximum of 200 jobs
		then
			echo "starting sbatch ${FILE} "
			sbatch "${FILE}"
			waiting="false";
			break;
		fi
		sleep 0.1 # sleep 0.1 seconds
	done
done

sleep 10

waiting=true
while "$waiting" -eq "true"
do
	squeue -u $user -o %j_%i_%n | egrep 'kallisto|STAR' &> currentJobs.txt 
	nodeJobs=$(grep "$nodeName" currentJobs.txt | egrep 'kallisto|STAR' | wc -l)
	if [[ "$nodeJobs" -lt 1 ]] # as long as there are jobs keep waiting
	then
		waiting="false";
		break;
	fi
	sleep 5 # sleep 5 seconds
done
echo "5="
echo $5
echo "2=$2"
find $2 -name '*.err' -exec grep -H "pseudoaligned" {} \; | tr ' ' '\t' | cut -f1,3,5 | sed -e 's/\[quant\]//g' > $4
mail -s "done!" $5 <<< "Finished:\n find $2 -name '*.err' -exec grep -H "pseudoaligned" {} \; | tr ' ' '\t' | cut -f1,3,5 | sed -e 's/\[quant\]//g' > $4/mappingPerSample.txt"
