#RNAseqPipe3 README

---

For detailed descriptions of the tools implemented in this pipeline see [documentation by Bec](https://github.com/pedrocrisp/NGS-pipelines/blob/master/Docs/RNAseq.md)

---

Pipeline to process RNAseq data

Reads must be in a folder called reads, with subdirectories for each sample.  Each sample directory should contain one fasatq.gz; if there are multiple fastq files cat them together using 10-runner.sh (10-cat.sh) called as below.

Scripts are called from the parent directory of the reads folder

Each script and std error is printed to the logs directory

Each step is called in order and run in bash eg 

```
bash ~/path_to_01.runner.sh
```

To run all the main pathway script call the 0*-runner.sh at once, eg:

```
for script in ~/<path_to_pipe3>/0*runner.sh ; do echo $script ; bash $script ; done
```

---

Steps 02-scythe.sh requires the file "truseq_adapters.fasta" to be in the script directory or a symbolic link to it

Steps 05-subread.sh and 06-featureCounts.sh require a symbolic link called "subread\_refdir" to the folder containing the index files, those files must have the prefix "TAIR10\_gen\_chrc"

On the server:

```
ln -s /home/pete/workspace/refseqs/TAIR10_gen/ subread_refdir

ln -s /home/pete/workspace/refseqs/TAIR10_gen/TAIR10_gen_chrc.chrom.sizes TAIR10_gen_chrc.chrom.sizes

ln -s /home/pete/workspace/Illumina/truseq_adapters.fasta truseq_adapters.fasta

```

The default 06-featureCounts.sh assumes stranded RNAseq.  To run without stranded-ness call 17-featureCountsNosense.sh or to align to align 'anti-sense' (switch the strand of each feature) call 16-featureCountsAntisense.sh.
