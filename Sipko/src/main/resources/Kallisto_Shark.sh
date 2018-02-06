#!/bin/bash
echo "0=$0"
echo "1=$1"
echo "2=$2"
echo "3=$3"
echo "4=$4"
echo "5=$5"
echo "6=$6"

user=$3
#get the location of the file containing all the scripts that have already been starten (creates a new file here if it does not already exist)
startedJobs=$7
userJobFile="currentJobs_$user.txt"
#loop through all the .sh files
for FILE in $1
do
	#Submit the job to the queue
	echo "starting qsub ${FILE}" 
	qsub "${FILE}" 
done

#sleep a bit extra so all jobs will be listed when querying the current queue
sleep 10

#wait for all jobs to finish
waiting=true
while "$waiting" -eq "true"
do
	
	nodeJobs=$(qstat | grep $user | egrep -i 'kallistoJob' | wc -l)
	if [[ "$nodeJobs" -lt 1 ]] # as long as there are jobs keep waiting
	then
		waiting="false";
		break;
	fi
	sleep 5 # sleep 5 seconds
done

#all jobs are finished, send an email to the user
mail -s "done!" $5 <<< "Finished' > $4"
