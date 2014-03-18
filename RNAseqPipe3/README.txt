#pipe3 README

pipeline to process RNAseq data

Reads must be in a folder called reads and the scripts must be called from the parent directory of the reads folder

Each step has to be called in order and run in bash eg bash ~/path_to_01.runner.sh

To run all the main pathway script call the 0*-runner.sh at once, eg:
for script in ~/<path_to_pipe3>/0*runner.sh ; do echo $script ; bash $script ; done

Steps 02-scythe.sh requires the file "truseq_adapters.fasta" to be in the script directory or a symbolic link to it

Steps 05-subread.sh and 06-featureCounts.sh require a symbolic link called "subread_refdir" to the folder containing the index files, those files must have the prefix "TAIR10_gen_chrc"

