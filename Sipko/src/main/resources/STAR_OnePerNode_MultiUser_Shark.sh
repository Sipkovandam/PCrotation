#!/bin/bash
echo "0=$0" #Location of this script
echo "1=$1" #Folder where all the .sh scripts are located
echo "2=$2" #Folder where the STAR results are put
echo "3=$3" #SLURM username
echo "4=$4" #Filename of the file that is going to contain all the mapping percentages (only works for Kallisto, so does nothing here)
echo "5=$5" #Email address used to mail the user when the script is done
echo "6=$6" #Command used to execute STAR
echo "7=$7" #Filename of the file containing which scripts are already ran

#get the username (e.g. umcg-svandam)
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
	
	nodeJobs=$(qstat | grep $user | egrep 'sh_STAR' | wc -l)
	if [[ "$nodeJobs" -lt 1 ]] # as long as there are jobs keep waiting
	then
		waiting="false";
		break;
	fi
	sleep 5 # sleep 5 seconds
done

#all jobs are finished, send an email to the user
mail -s "done!" $5 <<< "Finished' > $4"
