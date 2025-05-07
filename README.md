# RNA seq pipeline with Nextflow, Apptainer, and PBS 

## Tools and their pararmeters configuration 

This is a list of all tools used in this project and explanation for their usage: 

1. Fastqc

```
fastqc -t 2 -o fastqc_${sample_id}_logs -f fastq -q ${reads}
```
__-t__: works on file basis i.e one _thread_ per file. 

2. [Cutadapt](https://cutadapt.readthedocs.io/en/v4.8/guide.html#basic-usage)

```
cutadapt -q 20 -m 20 -j 8 -o trimmed_${sample_id}_1.fastq -p trimmed_${sample_id}_2.fastq $reads
```

Whether an input file needs to be decompressed or an output file needs to be compressed is detected automatically by inspecting the file name. Because all output files are short-lived intermediate files, so they are not compressed to speed up the process (output file not ending with .gz). 

__-q__ (or --quality-cutoff): trim low-quality ends from reads. If you specify a single cutoff value, the 3’ end of each read is trimmed. For Illumina reads, this is sufficient as their quality is high at the beginning, but degrades towards the 3’ end.

__-m__: Discard processed reads that are shorter than LENGTH.

__-j 8__: Increasing the number of _cores_ with -j will increase the number of reads per minute at near-linear rate.

__-p__: By default, all processed reads, no matter whether they were trimmed or not, are written to the output file specified by the -o option (or to standard output if -o was not provided). For paired-end reads, the second read in a pair is always written to the file specified by the -p option.

3. [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

```
bowtie2 -p 8 --very-sensitive -x bowtie2_index/ggal -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} -S mapped_${sample_id}.sam
```
__-p 8__: _Threads_ will run on separate processors/cores and synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing -p increases Bowtie 2's memory footprint. E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint by a few hundred megabytes. This option is only available if bowtie is linked with the pthreads library (i.e. if BOWTIE_PTHREADS=0 is not specified at build time).

4. [Samtools](https://www.htslib.org/doc/samtools-sort.html)

```
samtools sort -n -o mapped_${sample_id}_sorted.bam ${bam} 
```

__-n__: Sort by read names (i.e., the QNAME field) using an alpha-numeric ordering, rather than by chromosomal coordinates. 

5. [Htseq](https://htseq.readthedocs.io/en/release_0.11.1/count.html)

```
htseq-count -n 2 -f bam -r name -m union -s reverse -t exon -i gene_id mapped_${sample_id}_sorted.bam ${gff} > counts_${sample_id}.txt
```

__-r name__: For paired-end data, the alignment have to be sorted either by read name or by alignment position. 

__-s reverse__: The second read has to be on the same strand and the first read on the opposite strand. 

6. [Deseq2](https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html)

## Resources request

There is a need to know approximately what resources our jobs will need to configure _cpus_ and _memory_ for each PROCESS in Nextflow.

### Method 1: Submit a seperate job script

1. Run the command(s) in a PROCESS on one sample's pair-end read as a seperate job. For example, this is a job script to map reads (I ran code under /scratch because it is a fast local SSD disks)
```
#!/bin/bash

# Run this as "qsub submit_mapping.sh"

# Set a name for this run and the resource requirements,
# 8 CPU, 4 GB memory and 40 minutes wall time.
#PBS -N test_mapping
#PBS -l ncpus=8
#PBS -l mem=4GB
#PBS -l walltime=00:40:00

# Create a unique /scratch directory.
SCRATCH="/scratch/${USER}_${PBS_JOBID%.*}"
mkdir ${SCRATCH}

# Change to the input files path
cd ${PBS_O_WORKDIR}

# Copy your input data to this scratch directory.
cp <reads 1> <reads 2> <reference genome> ${SCRATCH}

# Change directory to the scratch directory and run your program.
cd ${SCRATCH}

# Run your program.
mkdir bowtie2_index
cd bowtie2_index
bowtie2-build ../<reference genome> ggal 
cd .. 
bowtie2 -p 8 --very-sensitive -x bowtie2_index/ggal -1 <reads 1> -2 <reads 2> -S mapped.sam
samtools view -S -b mapped.sam > mapped.bam

# Copy output results back to your working directory. 
mv ${SCRATCH}/mapped.bam ${PBS_O_WORKDIR}/output/

# Clean up
cd ${PBS_O_WORKDIR}
rm -rf ${SCRATCH}
```

2. Submit your job and save the job id

```
qsub submit_fastqc.sh
```

3. Run the following to see job stats from which we know how much RAM and how many cpus are actually used. 

```
qstat -fx <job id>
```
If the resources_used.mem < Resource_List.mem or resources_used.ncpus < Resource_List.ncpus, you should request less memory or cpus. 

### Method 2: Run Nextflow

Alternatively, we can run the Nextflow pipeline on one sample's pair-end reads, then repeat step 2 and 3 in Method 1. 

## Tips for effiency 

- Nextflow queueSize directive: The queueSize directive is part of the executor configuration in the nextflow.config file, and defines how many processes are queued at a given time. By default, Nextflow will submit up to 100 jobs at a time for execution. Increase or decrease this setting depending your HPC system quota and throughput. 

- Nextflow scratch directive: A common recommendation is to use the node's local scratch storage as the job working directory to avoid unnecessary use of the network shared file system and achieve better performance.

- Job submission limit: It is recommended that you configure your workflow to specify small, short jobs as using the local executor, leaving the larger and longer running jobs to the slurm executor.

## Run the program

1. Download the sample datasets here, which includes 4 pair-ended datasets - 2 treated and 2 untreated and save to the project directory.

2. Navigate terminal to the project directory.

3. Build apptainer image 

```
cd apptainerdef
apptainer build speedx-rnaseq.sif speedx-rnaseq.def
```

__(Optional)__ To check the if the apptainer works by starting an instance:

```
apptainer instance start speedx-rnaseq.sif speedx-rnaseq
```

Then, check if all software are installed by instance shell:

```
apptainer shell instance://speedx-rnaseq
```

Once everthing is verified, stop the instance: 
```
apptainer instance stop speedx-rnaseq
```

4. Run the Nextflow pipeline

```
nextflow run main.nf -profile apptainer
```

5. Check jobs submitted 

## Pipeline explanation 

With _process.container = '$baseDir/apptainerdef/speedx-rnaseq.sif', process.executor = 'pbspro', apptainer.enabled = true, process.scratch = /scratch/u162557, Nextflow: 

1. Submits a job to a computing note. 
2. In the computing note, creates a unique directory whose path is assigned to process.scratch. Here, I set process.scratch=/scratch/u162557.
3. Inside this scratch directory, creates a symlink for each input file required by the job execution. 
4. Mounts the computing node's scratch and current working directories to corresponding directories in the Apptainer container. 
5. Starts an Apptainer container and runs the commands in scratch directory of container. Hence, the outputs will appear in scratch of host as well. 
6. Move the output files from scratch to the shared working directory.
7. Deletes the scratch directory in the computing node.

## To do 
1. Extract, annotate, visualise DE genes
2. GO analysis
3. KEGG path analysis

## Resources: 

https://seqera.io/blog/5-more-tips-for-nextflow-user-on-hpc/
https://docs.ycrc.yale.edu/clusters-at-yale/guides/nextflow/
https://seqera.io/blog/5_tips_for_hpc_users/

