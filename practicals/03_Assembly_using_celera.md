#Assembly using celera assembler

Celera Assembler (CA) is the oldest assembly program we will use. It was developed during the time of Sanger sequencing by the company Celera Genomics. Celera Assembler was used to assemble the *Drosophila* genome, as well as the human genome.

During the later years, the developers have continuously updated the program to be able to work with 454 and Illumina data, and most recently wth PacBio reads.

This tutorial involves two assemblies, one using the same Illumina data as for Velvet and SPADES, the other using PacBio data (see below). Even if you only intend to assemble one of the datasets, *please read this entire document*.

###Celera assembler and Illumina data


####Preparing the data

Celera uses its own data format for input reads, which often only contains a pointer to an input file, and some metadata. The `fastqToCA` command will produce the right input:

* `-insertsize` gives the average insert length for the library, as well as the expected standard deviation of that average
*  `-type sanger` tells Celera that 'quality values are PHRED, offset=33'
*  the other parameters are self-explanatory

To prepare the Illumina paired end and mate pair reads, the commands are:

```
fastqToCA -insertsize 300 30 \
-libraryname ecoli_frag -technology illumina -type sanger \
-mates /data/assembly/MiSeq_Ecoli_MG1655_50x_R1.fastq,/data/assembly/MiSeq_Ecoli_MG1655_50x_R2.fastq \
> ecoli_frag.frg

fastqToCA -insertsize 3000 600 -libraryname ecoli_jump \
-technology illumina -type sanger\
-mates /data/assembly/Nextera_MP_R1_50x.fastq,/data/assembly/Nextera_MP_R2_50x.fastq \
>ecoli_jump.frg
```

####Assembly

The `runCA` command starts the process. It can take many parameters to fine-tune assembly, set up memory and CPU use etc.  These parameters are usually added in a separate `.spec` file. Here, we do not need to do that because we restrict ourselves to:

* `-p` name of assembly
* `-d` name of output folder (tip: use the same name as for `-p`)
* `unitigger=bogart` the optimal algorithm for the so-called unitig-stage, see the Celera website for details

We will save the output from the command using `>runCA.out` in a file to be able to follow progress. `2>&1` makes sure any error-messages are written to the same file.

Run the assembly as follows:

**NOTE** the assembly will take many hours, so use the `screen` command! See [https://wiki.uio.no/projects/clsi/index.php/Tip:using_screen](https://wiki.uio.no/projects/clsi/index.php/Tip:using_screen)

```
runCA -p <asm_name> -d <folder_name> unitigger=bogart *.frg >runCA.out 2>&1
```

###Celera assembler and PacBio data

The Pacific Biosciences instrument, the PacBio RS II, can produce very long reads, albeit at much lower per-base qualities than other, short-read platforms. The raw per-base qualitites are too low to be able to assemble the reads directly (at least, as of today). To overcome this problem, different error-correction methods have been developed. For bacterial genomes, the most promising method needs a high coverage dataset (60-100 x coverag), and uses the shortest raw PacBio reads in the data, to error-correct the longest ones. The result is 20-30x coverage in high-quality very long reads. Currently, only MIRA and Celera Assemlber are abel to use such long reads for assembly.  
Unfortunately, error-correction takes a lot of memory and compute time. For this tutorial, you will start with reads that already have been error-corrected.

The metrics for the *input* read data are:

	Count	18354
	Sum	 bases 112 Mbp
	Average length	6093 bp
	N50 length	6472
	Largest	14592

####Preparing the data
As with the illumina data (see above) we will use fastqToCA. `-technology` is set to `sanger` for this kind of reads as there is no separate flag for PacBio reads and they can be treated as regular sanger-type reads.


```
fastqToCA -libraryname pacbioReads \
-technology sanger -type sanger \
-reads /data/assembly/PB_corrected.fastq \
 >PB_corrected.frg
```

####Assembly

In this case, we will provide a `spec` with parameters. This file is full of settings which are beyond this tutorial to explain. Please see the manual. 

First, copy the spec file:

```
cp /data/assembly/celera_pb.spec .
```

**NOTE** the `.` at the end of the command, this indicates 'the current directory'.

Feel free to have a look at this file.

Now we are ready to run celera.  
**NOTE** the assembly will take many hours, so use the `screen` command! See [https://wiki.uio.no/projects/clsi/index.php/Tip:using_screen](https://wiki.uio.no/projects/clsi/index.php/Tip:using_screen)

```
runCA -p <asm_name> -d <folder_name> -s celera_pb.spec PB_corrected.frg >runCA.out 2>&1
```

###Celera output
Celera generates a lot of files and folders, more than 20 GB for the illumina assembly. The most important folder is `<asm_name>/9-terminator/`. Here you will find:

* `<asm_name>.scf.fasta` --> scaffold sequences
* `<asm_name>.ctg.fasta` --> contig sequences
* many other files, see the manual.
