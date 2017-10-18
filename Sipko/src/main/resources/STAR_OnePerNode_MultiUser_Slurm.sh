#!/bin/bash
echo "0=$0" #Location of this script
echo "1=$1" #Folder where all the .sh scripts are located
echo "2=$2" #Folder where the STAR results are put
echo "3=$3" #SLURM username
echo "4=$4" #Filename of the file that is going to contain all the mapping percentages (only works for Kallisto, so does nothing here)
echo "5=$5" #Email address used to mail the user when the script is done
echo "6=$6" #Command used to execute STAR
echo "7=$7" #Filename of the file containing which scripts are already ran

#list of nodes that can be used. This script will fail if nodes are listed that do not exist or can't be used
nodenames=(umcg-node011 umcg-node012 umcg-node013 umcg-node014 umcg-node015 umcg-node017 umcg-node018)

#get the username (e.g. umcg-svandam)
user=$3
#get the location of the file containing all the scripts that have already been starten (creates a new file here if it does not already exist)
startedJobs=$7
userJobFile="currentJobs_$user.txt"
#loop through all the .sh files
for FILE in $1
do
	waiting=true
	while "$waiting" -eq "true"
	do
		#wait untill the jobs for this user are actually running before submitting new ones
			#This is to prevent multiple jobs occuplying nodes, blocking other users from helping you run your script (only one STAR job is scheduled per node since STAR becomes slower if you run multiple at the same time on the same node)
		jobsWaiting=$(squeue -u $user | egrep 'Priority|QOSGrpMemLi|QOSMaxMemoryPerUser' | wc -l)
		while [[ "$jobsWaiting" -gt 2 ]]
		do
			dateStamp=$(date)
			echo $dateStamp'	jobswaiting='$jobsWaiting
			sleep 10
			jobsWaiting=$(squeue -u $user | egrep 'Priority|QOSGrpMemLi|QOSMaxMemoryPerUser' | wc -l)
		done
		#get all the jobs running and on which node they are running/scheduled

		squeue -o %j_%i_%n_%t_%M_%N_%Y | egrep 'sh_STAR'  &> $userJobFile
		nodeName=${nodenames[$a]}
		
		#check if there is a job scheduled for that node already
		nodeJobs=$(egrep "$nodeName" $userJobFile | wc -l)
		if [[ "$nodeJobs" -lt 1 ]] # queue a maximum of 1 job per node
		then
			#check if the job was already started by another user or in a previous run
				#if so go to the next file
			isrunning=$(grep "$FILE" "$startedJobs" | wc -l)
			if [[ "$isrunning" -gt 0  ]]
			then
				echo "job already started;"$FILE" skipping"
				waiting="false";
				break;
			fi
			

			#Submit the job to the slurm queue
			echo "starting --nodelist $nodeName sbatch ${FILE} nodeJobs=$nodeJobs nodeName=$nodeName" 
			sbatch --nodelist $nodeName "${FILE}" 
			
			#Add the file to the list of scripts that has already been executed
			echo $FILE >> $startedJobs
			waiting="false";
			
			#Get the next nodename
			a=$(($a+1));
			if [ $a -eq ${#nodenames[*]} ] ; then a=0; fi
			
			#go to the next file
			break;
		fi
		#Get the next nodename
		a=$(($a+1));
		if [ $a -eq ${#nodenames[*]} ] ; then a=0; fi
		sleep 0.1 # sleep 0.1 seconds
	done
	#sleep for a bit after submitting a new job
	if [[ "$isrunning" -lt 1  ]]
		then
		sleep 5 # sleep 5 seconds
	fi
done

#sleep a bit extra so all jobs will be listed when querying the current queue
sleep 10

#wait for all jobs to finish
waiting=true
while "$waiting" -eq "true"
do
	squeue -u $user -o %j_%i_%n &> $userJobFile 
	nodeJobs=$(grep "$nodeName" $userJobFile | egrep 'sh_STAR' | wc -l)
	if [[ "$nodeJobs" -lt 1 ]] # as long as there are jobs keep waiting
	then
		waiting="false";
		break;
	fi
	sleep 5 # sleep 5 seconds
done

#all jobs are finished, send an email to the user
echo "5="
echo $5
echo "2=$2"
find $2 -name '*.err.gz' -exec zcat {} \;| grep -e '1:\|processe' | tr ' ' '\t' | sed -e ':a;N;$!ba;s/_1.fq.gz\n/\//g' | cut -f6,8,10| sed -e 's/\[quant\]//g' > $4
mail -s "done!" $5 <<< "Finished:\n find $2 -name '*.err.gz' -exec zcat {} \;| grep -e '1:\|processe' | tr ' ' '\t' | sed -e ':a;N;$!ba;s/_1.fq.gz\n/\//g' | cut -f6,8,10| sed -e 's/\[quant\]//g' > $4"
