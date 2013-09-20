#Assembly using SPADES

Spades was written as an assembly program for bacterial genomes, from regular, as well as from whole-genome amplified samples. It performed very well in the GAGE-B competition, see [http://ccb.jhu.edu/gage_b/](http://ccb.jhu.edu/gage_b/).

Before assembly, SPADES will error-correct the reads.

NOTE that SPADES is not optimised for scaffolding using mate pairs. 

###Using SPADES

Spades can be used with paired end and mate pair data:

* The `--careful` flag is used to reduce the number of mismatches and short indels. 
* For each read file, a flag is used to indicate whether it is from a paired end (`--pe`) or mate (`--mp`) pair dataset, followed by a number for the dataset, and a number for read1 or read2. For example: `--pe1-1` and `--pe1-2` indicate pared end data set 1, read1 and read2, respectively.
* Similarly, use `--mp-1-1` and `--mp1-2` for the mate pair files. 
* Spades assumes mate pairs are in the orientation as they are in the original files coming from the Illumina instrument: <-- and --> ('outie' orientation, or 'rf' for reverse-forward). Our reads are in the --> and <-- ('innie', 'fr' for forward-reverse) orientation, so we add the `--mp1-fr` flag to let SPADES know about this
  
Other parameters:

* `-t` number of threads (CPUs) to use for calculations
* `-k` k-mers to use (this gives room for experimenting!)
* `-o` name of the output folder

####Setting up the assembly

First, create a new folder called `/home/<your_username>/assembly/spades` and `cd` into it.  
We will save the output from the command using `>spades.out` in a file to be able to follow progress. `2>&1` makes sure any error-messages are written to the same file.
Run the assembly as follows:

**NOTE** the assembly will take around two hours, so use the `screen` command! See [https://wiki.uio.no/projects/clsi/index.php/Tip:using_screen](https://wiki.uio.no/projects/clsi/index.php/Tip:using_screen)

```
spades.py -t 2 -k 21,33,55,77 --careful \
--pe1-1 /data/assembly/MiSeq_Ecoli_MG1655_50x_R1.fastq \
--pe1-2 /data/assembly/MiSeq_Ecoli_MG1655_50x_R2.fastq \
--mp1-1 /data/assembly/Nextera_MP_R1_50x.fastq \
--mp1-2 /data/assembly/Nextera_MP_R2_50x.fastq \
--mp1-fr -o <asm_name> >spades.out 2>&1
```

If the assembly is running in a 'screen', you can follow the output by checking the `spades.out` file.  

**TIP**: use this command to track the output as it is added to the file. Use ctrl-c to cancel.

```
tail -f spades.out
```

####SPADES output
* error-corrected reads
* contigs for each individual k-mer assembly
* final `contigs.fasta` and `scaffolds.fasta`

####Re-using error-corrected reads

Once you have run SPADES, you will have files with the error-corrected reads in `spades_folder/corrected/`. There will be one file for each input file, and one additional one for unpaired reads (where during correction, one of the pairs was removed from the dataset). Instead of running the full SPADES pipeline for your next assembly, you could add the error-corrected reads from the previous assembly. This will save time by skipping the error-correction step. I suggest to not include the files with unpaired reads.

Error-corrected read files are compressed, but SPADES will accept them as such (no need to uncompress).

Changes to the command line when using error-corrected reads:

* point to the error-corrected read files instead of the raw read files
* add the `--only-assembler` flag to skip correction


**Questions:**

* How does the assembly compare to the velvet assemblies using the (exact) same input data
* Which assembler is 'best' for this data - or can't you tell?

**Bonus exercise**: map the *error-corrected* paired end reads to the same assembly as the uncorrected reads, and check the results in the genome browser.

###Next steps
As for the previous assemblies, you could map reads back to the assembly, calls SNPs and Indels, run FRC and Quast, and visualise in the browser.
