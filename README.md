# 2021_HW4_Gene-Expression

# Homework #3A – Aligning & counting reads for gene expression analysis

**DUE: 11:59PM, 6 November** 

## Description:
You’ll be conducting a relatively simple comparison of gene expression between two treatment groups, using short-read mRNA-seq data and a (provided) *de novo* transcriptome. This involves several overall steps: cleaning raw sequence reads, aligning reads to the transcriptome, counting how many reads align to each contig, and testing for differential expression (≈ different numbers of mapped reads) between treatment groups.

For each of these steps, there are a number of different programs / approaches that can be strung together into an analysis “pipeline”. All have their advantages and drawbacks, but there’s no single objective “best” way to analyze these data. For simplicity, you'll all be using the "wicked fast" pseudo-aligner `salmon`, but we will assign you a mix of DE analysis programs. For this homework, please read some of the documentation for your assigned program and run the data through them. Record your commands and note the metrics and output of your pipeline at each stage of the process. Please annotate this extensively with your thoughts on this pipeline, guided (but not limited!) by the questions below.

Once homework has been submitted, we will collate and compare the results of your various analyses and present them in a later class. We’ll discuss the results as a class, so as part of this assignment be prepared to discuss the particulars of your program, pros / cons, etc.

Start this homework early! There are several different steps, and this will give you time to troubleshoot and ask questions as needed.

This homework has quite a few steps. First, you will start with raw reads and a provided transcriptome assembly and prepare these data for differential expression analysis. This involves retrieving reads from the SRA, cleaning reads, aligning reads to the transcriptome, counting the number of reads aligned to each contig, and reformatting these count data to input into a differential analysis program. Next, you'll use your DE analysis program to compare expression between the two groups, do a little functional enrichment analysis using a provided files of GO terms, and answer some questions about your results and what they mean.

Note: If you're not formally enrolled but want to do this homework, please pick whichever analysis program strikes your fancy.

Please copy this document and change the name to `hw3_answers_[LASTNAME].md`, and reply in the document as prompted. We'll also ask you to create and save several scripts to be pushed with your HW repo. Everything we want you to submit is listed at the bottom of the document.

## The Data:

The data you'll be working with are from DeLeo & Bracken-Grissom (paper included in this repo):
* Species: Vertically-migrating shrimp, *Systellaspis debilis*
* 'Treatment' groups: day / deep and night / shallow
* Experiment details (abridged from paper):
Briefly, *S. debilis* were collected in the Florida Strait in July 2017 with "a 9-metre Tucker Trawl fitted with a light-tight, thermally insulated cod-end that could be opened and closed at depth. This method enabled specimen collection from specific depth intervals and maintenance at in situ temperatures prior to preservation. At the surface, species were identified under dim red light to avoid any damage to photosensitive tissues." "Day samples were collected in the morning/afternoon (presunset) from ~450 to  750 m, and night samples were collected around midnight (predawn) from  ~150 to 330 m." "Eye tissues were carefully dissected under a dissecting scope while submerged in RNAlater from five biological replicates corresponding to each sampling condition, day (n = 5) and night (n = 5)."
* Sampling structure (abridged from paper):
Five individuals captured at each time of day / depth, for a total of 10 individuals sequenced.
* Sequencing: mRNA sequencing, 150 bp paired-end reads, Illumina HiSeq

## Provided files:

Downsampled raw transcriptome sequence files from 10 samples, with two files per sample (forward and reverse reads). NOTE: These samples are NOT in your GitHub repo, since they're rather large. More on this in Step 0.

A _de novo_ transcriptome assembled by the authors for this project.

A file linking each contig to a series of GO terms for functional enrichment analysis.

An empty slurm script called `hw3_slurm_wrapper.txt` with placeholder info. Feel free to use this wrapper as the starting point for the (several!) slurm submissions you will do in the course of this homework.

## STEP 0: Get raw reads

Because of their size, the sequence files you'll be using for this project aren't in the GitHub repo. While they are technically small enough to be allowed, they take up enough space that we don't want each of you making a separate copy on the HPC, where space is precious. Instead, they're already on the HPC in a directory called `HW3_sequences` in the `collaboration` directory of the class directory (same place all your project directories are located; you should all have full access here). 

You are going to access these samples using a "symbolic link", a neat trick to let everybody work with the **same** reference samples in a central location accessible to all. This is especially handy for raw data or databases: files that are large and that will be accessed repeatedly and/or by multiple users.

To create a symbolic link to the raw samples, do this **from your HW3 directory on your personal space**:

`ln -s [FULL_PATH_TO_SEQUENCES] HW3_raw_seqs`

Double-check that the path to the sequences (in the `collaboration` directory) is correct! Use an absolute, not a relative path here.

You should see a new "file" named `HW3_raw_seqs` in your HW3 directory. Try `ls HW3_raw_seqs`. What happens? You should see a list of 20 gzipped samples. If you do - awesome! If not, double-check your path and try again (and maybe google around for more guidance on creating and using symbolic links).

Now you can use `HW3_raw_seqs` as if it were a directory with all the files in them. Basically, it points right to the single centralized repository and lets you use those samples without having to duplicate them in your own space.

## STEP 1: Clean raw reads

Raw reads off a sequencer are not always perfect – they can have some low-quality bases, and sometimes the adaptor sequence hasn’t been completely removed. As a first step, we’re going to clean up our raw sequence reads so we are working only with high-quality base calls. (This is usually your first step no matter what you’re doing downstream.)

There are a number of ways to do this, but we’re going to use Trim Galore! (The ! is part of the official name.) It is a wrapper script for cutadapt and fastqc, which it requires to run. 

Set up a conda environment called `HW3` for this homework, and add **all** necessary programs **except those that run in R** to this environment as you go through the homework. You can run R programs in conda, but it can be complicated to get everything installed poperly and it is awkward to create and view figures, so we're going to let you run R on your own computers if you prefer.

Here is the command for cleaning a single file:

`trim_galore -q 20 --phred33 --illumina -stringency 3 --paired --length 100 --fastqc [SEQFILE_FORWARD] [SEQFILE_REVERSE]`

What do the flags on this command mean? (Hint: the program output might be helpful here.)

> Answer: 

Write a bash loop to trim and clip all samples, embed it in a slurm script called `hw3_clean-reads_[LASTNAME].txt`, and run it on the HPC. Use the following parameters:

```
partition=compute
cpus-per-task=1
mem=5000
time=2:00:00
```

Push your `hw3_clean-reads_[LASTNAME].txt` slurm script to GitHub as part part of your finished homework.

## STEP 3: Alignment

Next, you need to align the individual reads for each sample to the transcriptome. Use `salmon` - it is `wicked fast` and is commonly used for gene expression analysis. Please check out the help files and other online resources to figure out how to use it.

Does salmon do “proper” alignment, or pseudo-alignment?  What's the difference?
> Answer: 

Could you use `salmon` to align transcriptome reads back to a **genome**?
> Answer: 

What resources did you use to figure out how to run `salmon` (please give URLs, etc)?
> Answer: 


Typically, you have to index your reference transcriptome (or genome). Use as a reference transcriptome `hw3_copepod_txm.fasta`, which is included as a gzip file in this repo. What command(s), if any, did you use to index the transcriptome for alignment?

```
```

Set up command(s) to align your trimmed reads to the indexed transcriptome, setting the following parameters (if possible):
`kmer = 25`

Write a script using Slurm called `hw3_align-[ALIGNER]\_[LASTNAME].txt` to run your alignment on Poseidon. Use the following slurm parameters:

```
partition=compute
ntasks=1
mem=5000
time=2:00:00
```

Push your `hw3_align-[ALIGNER]\_[LASTNAME].txt` slurm script to GitHub as part part of your finished homework. 

Copy the alignment command(s) for a single sample below, and explain each flag:
```
```

How much time did it take to run your alignment script? (Hint: check your slurm-generated email.)
> Answer: 

For sample SRR5004080, what % of reads aligned / mapped?
> Answer: 

Hint: This information will be automatically calculated and is lurking somewhere in the files you’ve already produced.

## STEP 3.5: Count aligned reads (**BWA & bowtie2 only**)

After reads are aligned, it’s time to extract quantification data. Many pseudo-aligners will do this automatically; “proper” aligners do not. Instead, they output sam (Sequence Alignment/Map) files with alignment info (sorry, folks with proper aligners!).

If you’re using `bowtie2` or `BWA`, modify and run the `hw3_count.txt` slurm script we’ve provided; call this modified script `hw3_count_[ALIGNER]\_[LASTNAME].txt`. You will need to add `Picard tools` and `samtools` to your conda environment. Push your personalized read-counting slurm script to GitHub as part part of your finished homework.

Note: In your logfile, you may get a series of warnings for the “CleanSam” command saying, “Error parsing text SAM file. MAPQ should be 0 for unmapped read.” Go ahead and ignore these, this has to do with how read mapping is coded by some aligners and doesn’t indicate any downstream issues.

## STEP 4: Looking at data & prepping input for differential expression analysis:

You should now have one counts file per sample giving the number of reads in that sample which mapped back to each contig in the reference transcriptome. (You may have other information as well, including TPM or contig length.) Using these count files as input, write a bash script to count how many contigs in each sample are covered by `> 10000` reads. The output should be the name of the file / sample, followed by the number of contigs with readcount `> 10000`. (It is fine if the sample name is on a different line from the count of highly-expressed contigs.) Paste this script below:
```
```

Paste the script output below:
```
```

Similarly, write a bash script to find the contig in each sample with the greatest number of reads mapping to it. The output should be the name of the file / sample, followed by contig name and the number of reads. (It is fine if the sample name is on a different line from the contig info.) Paste this script below:
```
```

Paste the script output below:
```
```

Use the `hw3_counts-to-table.py` python script (provided) to create a “master” tab-delimited table of read counts for all samples, where each row is a contig and each column is a sample. Each cell contains the raw (or imputed) read counts for that sample x contig combination. This will be the input file used in downstream differential expression analysis.

The top left corner of your  `.tsv` file should look something like this:
```
Contig  Sample_1    Sample_2    Sample_3
contig_1    13          27          4                
contig_2    456         1003        755
contig_3    0           1           3
```

Rename this file `hw3_counts_[ALIGNER]\_[LASTNAME].tsv`, and push it to GitHub as part of your finished homework.

Once you're finished with this homework, create a yaml file of your genex conda environment (while that environment is activated) like this:
`conda env export > hw3a_conda_genex_[LASTNAME].yml`

Congratulations! You are ready to proceed to differential expression analysis in the next homework.

About how much time did you spend on this homework?

For your homework, please push to GitHub:

1. `hw3a_answers_[LASTNAME].md`: An annotated copy of this readme file including your answers.
2. `hw3a_conda_genex_[LASTNAME].yml`: yaml file of your genex conda environment
3. `hw3_sra-retrieve_[LASTNAME].txt`: slurm script for retrieving reads from the SRA
4. `hw3_clean-reads_[LASTNAME].txt`: slurm script for cleaning raw reads
5. `hw3_align-[ALIGNER]\_[LASTNAME].txt`: slurm script for aligning clean reads to your transcriptome
6. `hw3_count_[ALIGNER]\_[LASTNAME].txt`: slurm script for translating alignments into a table of read counts (BWA / bowtie2 only) 
7. `hw3_counts_[ALIGNER]\_[LASTNAME].csv`: "master" contig-by-sample read count table

