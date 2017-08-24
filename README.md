# Geovars Pipeline Overview
This is an overview of the pipeline used to produce the structure graphs that showed clear populations

##### Before beginning
#
You'll need to do the following before begining
1. Open an instance of your preffered shell
2. Determine which hpc machine you want to use by 
    * Remotely log into each machine with *ssh* by typing `ssh netid@hpc_machine.extension` where `netid` is your ISU net id, `hpc_machine` is the name of whatever hpc resource you wish to use (i.e. speedy or rit2, etc), and `extension` is the name of the department that the machine is associated with (las for rit2 or ent for speedy). An example of this is `jtwalker@rit2.las.iastate.edu`
    * Run *htop* to determine which hpc machine is currently being used the least by typing `htop`. To quit *htop*, type `q`. [Guide to understanding htop output](https://peteris.rocks/blog/htop/)
3. Make a working directory
    * This can be created with `mkdir new_directory_name`
4. Copy all geovars sample files nescessary for analysis to the newly created working directory
    * `cp /path_to_bams/*.bam /home/netid/working_directory/`
#
It is also very useful to run things on a `screen`. The command to initiate a `screen` is simply `screen`. If you get kicked off of a remote machine or don't want to run things on the head node (please don't run things on the head node) use a screen. 
* To detatch from a screen so that you can come back to it later hold the `CONTROL` key and press `a` and `d`.
* To completely exit a `screen` session type `exit`

Another useful command to remember is `CONTROL-c` which will kill whatever process you are currently running
#
##### Steps completed by Novogene
1. Samples are sequenced, raw `.fastq`, clean `.fastq`, and reference genome aligned `.bam` files are returned
#
#
The following are the steps that I took in order to produce the *STRUCTURE* graphs which showed nice distinction between populations
##### Our pipeline
2. `.bam` files are merged into single '.bam'
3. Variants are called from merged `.bam`, producing a variant call file `.vcf`
4. `.vcf` is filtered to remove sites with low coverage, sites with low number of genotyped individuals, etc. 
5. Structure is run to evaluate population structure
6. Best *k* is evaluated with *Evanno method*
#
#
The `.bam` files that Novogene produced (using the *Burrows-Wheeler aligner*, or *BWA-mem*, are a good starting place, as they have supposedly had bad reads, over-represented sequences, etc. filtered from them. It might still be worth your time to look at the QC reports related to your samples in order to decide whether or not you want to start from initial `.fastq` files. 
#
#
##### 2. `.bam` files are merged into single '.bam'
#
#
```
module load picard
picard MergeSamFiles I=inputFile1.bam I=inputFile2.bam ... O=outputFile.bam
```
It seems that `Picard` will not take input in the form of a list of files to merge, so I have written the following *bash* script to create an input list and execute a `Picard` command which merges all `.bam` files within a given folder
```
#!/bin/bash # the shebang - tells the shell to use bash

module load picard # loads the picard module if you are on the ISU computing resources which use the linux program "modules"

ls *.bam > temp.sh # writes the names of all bam files in current directory to a list that is also a bash executable file

sed -i -e 's/^/I\=/' temp.sh # appends "I=" to the beginning of each line

sed -i '1ipicard MergeSamFiles' temp.sh # appends the program name and command to the top of the file

echo "O=output_all_bams.bam" >> temp.sh to run # appends an output command to the bottom of the temp file

sed -i -e ':a;N;$!ba;s/\n/ /g' temp.sh # removes all newline characters and replaces them with spaces so that the file now consists of only one line

bash temp.sh # runs the temporary file which we've created with this series of commands

rm temp.sh temp.she # cleans up mess
```
A file, `large_picard_job_creator.sh` in the associated Github repo contains this script without the possibly program breaking comments. To run this script you will need to place this file in your working directory and type `bash large_picard_job_creator.sh` (remember to use `TAB` autocomplete, as well as the `screen` program).

**THIS STEP WILL TAKE SEVERAL HOURS - POSSIBLY EVEN A DAY OR SO AND WILL TAKE A LARGE AMOUNT OF RAM  - BEST TO USE SPEEDY WHEN IT IS NOT BEING FULLY USED**
#
#
##### 3. Variants are called from merged `.bam`, producing a variant call file `.vcf`
#
#
Calling variants is a multistep process and can be a headache. Hopefully, numbering the steps in variant calling will make it a bit simpler for you
1. Load the program used for sorting and calling
`module load samtools`
2. Sort the merged `.bam` files. **This will take a LONG time, but is nescessary for steps down the pipe. This step may also error out if you don't have enough RAM available**
`samtools sort bam_from_previous_step.bam -o output_bam.bam`
3. Index the reference genome for use with *Samtools*. **This can be done at the same time as step 2, and also takes several hours to complete**
`samtools faidx /PATH_TO_FILE/GCF.fna`
4. Use *mpileup* to mark variants
`samtools mpileup -C 50 -q 1 -d 1000 -v -f /PATH_TO_REFERENCE/GCF.fna sorted_bam_file.bam > raw_vcf.vcf`
where `-C` adjusts mapping quality based on phred score
`-q` specifies minimum mapping quality
`-d` specifies the max depth that a variant position can have while still being included
`-v` specifies `.vcf` not binary `.bcf` output
`-f` points to fasta reference file
5. Load *BCFTools*
`module load bcftools`
6. Call variants (the previous step only marks variants, I guess?)
`bcftools call -mv raw_vcf.vcf -Ov -o final_vcf.vcf`
where `-Ov` specifies vcf output
`-o` points to output file
#
#
##### 4. `.vcf` is filtered to remove sites with low coverage, sites with low number of genotyped individuals, etc. 
#
#
This is one of the most important steps, but is, in my opinion, the most difficult. You'll need to remove variants that are errors and otherwise don't represent the biology of the samples, without removing variants that might alter your output in important ways. To filter use *Plink*. The general usage of plink is 
```
module load plink
plink --allow-extra-chr --input_file_type input_file.extension --filter_name parameter --desired_output_file_type -out output_filename
```
To produce the structure graphs that you've seen, I used `plink --vcf geovars.vcf --allow-extra-chr --maf 0.4 --geno 0.1 --recode -out output_file`
Where `--maf` specifies minor allele frequency to filter on
`--geno` specifies the maximum percent of missing genotypes allowed (I think I phrased this poorly, sorry)
`--recode` outputs `.ped` and `.fam` files respectively, which are what are fed into *PGDSpider* 
#
With a larger set of samples you'll need to play around with this more. Sorry that I can't offer more help
#
#
##### 5. *STRUCTURE* is run to evaluate population structure
#
#
You've likely already done the following before, but if you haven't used the HPC resources to run *STRUCTURE*, this might be useful
#
#
You'll need to download *PGDSpider* to your local machine if you don't already have it
[PGDSpider](http://www.cmpg.unibe.ch/software/PGDSpider/)
#
There is a `.spid` file that I created to make running this slightly less time consuming. It will be available on the Github associated with this doc. 
#
You'll also need to download the `.ped` and `.fam` files that you've generated with *Plink* to your local machine. You can use *secure copy* or *scp* to do this with the following command:
`scp netid@hpc-machine.extension.iastate.edu:~/path_to_working_directory/filename.extension /path_to_local_machine_working_directory/`
As an example
`scp jtwalker@rit2.las.iastate.edu:~/directory/file.ped /home/josh/directory_to_copy_to/`
* for speedy the extension is `.ent` rather than `.las`
#
After converting your `.ped` file into a `.strct_in` file via *PGDSpider* you'll need to use *scp* once more to copy the new `.strct_in` file to your working directory on whatever hpc machine you're using. The syntax for this is familiar (just reverse the source and destination addresses):
`scp /path_to_local_machine_working_directory/filename.ext netid@hpc_machine.extension.iastate.edu:~/path_to_working_directory/`
