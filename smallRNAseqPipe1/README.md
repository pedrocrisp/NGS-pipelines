#RNAseqPipe3 README

---

For detailed descriptions of some of the tools implemented in this pipeline see [documentation by Bec](https://github.com/pedrocrisp/NGS-pipelines/blob/master/Docs/RNAseq.md)

---

Pipeline to process sRNAseq data

Reads must be in a folder called reads, with subdirectories for each sample 

Scripts are called from the parent directory of the reads folder

Each script and std error is printed to the logs dircetory

Each step is called in order and run in bash eg 

```
bash ~/path_to_01.runner.sh
```

To run all the main pathway script call the 0*-runner.sh at once, eg:

```
for script in ~/<path_to_pipeline_folder>/0*runner.sh ; do echo $script ; bash $script ; done
```

---

Steps 02-scythe.sh requires the file "truseq_adapters.fasta" to be in the script directory or a symbolic link to it

Steps 03-cutadapt.sh requires the file "adapter.txt" to be in the working directory or a symbolic link to it

Steps 05-bowtie2.sh requires a symbolic link called "bowtie2_refdir" to the folder containing the index files, those files must have the prefix "TAIR10_allchr"

On the server:

```
ln -s /home/pete/workspace/refseqs/TAIR10_gen/ bowtie2_refdir

ln -s /home/pete/workspace/Illumina/truseq_sRNA_adapters.fasta truseq_adapters.fasta

```
