###
#code to make runner work on both osx and linux. Essentially, it creates a file path to the script directory and saves this path as $0. In detail: 'if the operating system type is darwin (a mac), then use the greadlink function when the readlink function is called. Then use the greadlink function to find the script file named. In doing so, find the path of the script files directory and save as 'scriptdir'. This change is made for Macs because readlink doesn't run properly, but greadlink does. If the OS is not mac (eg. Linux), find the script file using the readlink function and save the path to the script file directory as 'scriptdir.' By using readlink to find where the scripts are, it means if this pipeline is copied onto another computer, the files can still be found. 

if
[[ $OSTYPE == darwin* ]]
then
readlink=$(which greadlink)
scriptdir="$(dirname $($readlink -f $0))"
else
scriptdir="$(dirname $(readlink -f $0))"
fi
###

#user defined variables that can be changed if needed:
#script directory is set to the directory where this script is stored on network/computer/github. This will obviously change between runners. By not hardcoding this variable, the script can be used and opened on any computer by any person with access.
workingdir=./
script=$scriptdir/01-fastqc.sh

#output directory is set as reads_fastqc. This changes between runners.
outdir=reads_fastqc
###

#Defines the function 'findSamples' into the commandline, but does not run the function. The function 'findSamples' looks in the reads directory to find the sample directory names.
function findSamples () {
find reads/ -mindepth 1 -maxdepth 1 -type d  -exec basename {} \;| tr ' ' '\n'
}

#make a file called 'reads_fastqc.' This will be created in the currently open working directory. 
mkdir ${outdir}

#Create a timestamp including the date, with the order year, month, day, hrs, mins, seconds.
timestamp=$(date +%Y%m%d-%H%M%S)

#Within the current working directory, make a directory called 'logs'
mkdir logs

#Define the variable logdir, which is like making a shortcut to the path workingdir/logs/outdir and is timestamped.
logdir="./logs/${outdir}.${timestamp}"

#Create a log directory (located in the logdir path). This log folder will contain the logs for the fastqc run.
mkdir $logdir

#Writes the fastqc and the runner scripts (called script.log) to a text file in the logs directory, so you know what version of the script you used during your run.
cat $script > "$logdir/script.log"
cat $0 > "$logdir/runner.log"

#Prints the fastqc script on the screen.
cat $script

#This says 'carry out the findSamples function. Then run the bash script 01-fastqc.sh on each of these files in parallel (ie. using multiple cores on computer). Store a script of the run to the logs directory within the reads_fastqc folder. 
findSamples | parallel bash $script {} \>logs/${outdir}.${timestamp}/{}.log 2\>\&1

#To run, go to the reads directory and call:
#bash ~/path_to/01-runner.sh
