#!/bin/bash
echo "0=$0"
echo "1=$1"
echo "2=$2"
echo "3=$3"
echo "4=$4"
echo "5=$5"
nodenames=(umcg-node011 umcg-node012 umcg-node013 umcg-node014 umcg-node015 umcg-node016 umcg-node017 umcg-node018)

user=$3
for FILE in $1
do
	waiting=true
	while "$waiting" -eq "true"
	do
		squeue -u $user -o %j_%i_%n_%t_%M_%N_%Y | egrep 'Kallisto|STAR'  &> currentJobs.txt
		nodeName=${nodenames[$a]}
		
		nodeJobs=$(egrep "$nodeName|__\(null\)" currentJobs.txt | wc -l)
		if [[ "$nodeJobs" -lt 1 ]] # queue a maximum of 200 jobs
		then
			echo "starting --nodelist $nodeName sbatch ${FILE} nodeJobs=$nodeJobs nodeName=$nodeName"
			sbatch --nodelist $nodeName "${FILE}" 
			waiting="false";
			break;
		fi
		a=$(($a+1));
		if [ $a -eq ${#nodenames[*]} ] ; then a=0; fi
		sleep 0.1 # sleep 0.1 seconds
	done
	sleep 0.1 # sleep 1 seconds
done

sleep 10

waiting=true
while "$waiting" -eq "true"
do
	squeue -u $user -o %j_%i_%n | egrep 'Kallisto|STAR' &> currentJobs.txt 
	nodeJobs=$(grep "$nodeName" currentJobs.txt | egrep 'Kallisto|STAR' | wc -l)
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
find $2 -name '*.err.gz' -exec zcat {} \;| grep -e '1:\|processe' | tr ' ' '\t' | sed -e ':a;N;$!ba;s/_1.fq.gz\n/\//g' | cut -f6,8,10| sed -e 's/\[quant\]//g' > $4
mail -s "done!" $5 <<< "Finished:\n find $2 -name '*.err.gz' -exec zcat {} \;| grep -e '1:\|processe' | tr ' ' '\t' | sed -e ':a;N;$!ba;s/_1.fq.gz\n/\//g' | cut -f6,8,10| sed -e 's/\[quant\]//g' > $4"
