#Assembly using SPADES

Spades: [http://bioinf.spbau.ru/spades](http://bioinf.spbau.ru/spades)  
Manual: [http://spades.bioinf.spbau.ru/release2.5.1/manual.html](http://spades.bioinf.spbau.ru/release2.5.1/manual.html)

Spades was written as an assembly program for bacterial genomes, from regular, as well as from whole-genome amplified samples. It performed very well in the GAGE-B competition, see [http://ccb.jhu.edu/gage_b/](http://ccb.jhu.edu/gage_b/).

Before assembly, SPADES will error-correct the reads.

NOTE that SPADES is not optimised for scaffolding using mate pairs. 

###Using SPADES

Spades can be used with paired end data only, or with paired end and mate pair data. The `--careful`flag is used to reduce the number of mismatches and short indels. For each read file, a flag is used to indicate whether it is from a paired end (`--pe`) or mate (`--mp`) pair dataset, followed by a number for the dataset, and a number for read1 or read2, for example: `--pe1-1` and `--pe1-2`.  
Other parameters:

* `-t` number of threads to use for calculations
* `-k` k-mers to use (room for experimenting!)
* `-o` name of the output folder

####Paired end reads only

First, create a new folder called `/home/<your_username>/assembly/spades` and `cd`into it.
The, run SPADES:

```
spades.py -t 2 -k 21,33,55,77 --careful \
--pe1-1 /data/assembly/MiSeq_Ecoli_MG1655_50x_R1.fastq \
--pe1-2 /data/assembly/MiSeq_Ecoli_MG1655_50x_R2.fastq \
-o <asm_name>

```

####Paired end and mate pairs

Add the paired ends as above, and use `--mp-1-1` and `--mp1-2` for the mate pair files. Spades assumes mate pairs are in the orientation as they are in the original files coming from the Illumina instrument: <-- and --> ('outie' orientation, or 'rf' for reverse-forward). Our reads are in the --> and <-- ('innie', 'fr' for forward-reverse) orientation, so we add the `--mp1-fr` flag to let SPADES know about this:

```
spades.py -t 2 -k 21,33,55,77 --careful \
--pe1-1 /data/assembly/MiSeq_Ecoli_MG1655_50x_R1.fastq \
--pe1-2 /data/assembly/MiSeq_Ecoli_MG1655_50x_R2.fastq \
--mp1-1 /data/assembly/Nextera_MP_R1_50x.fastq \
--mp1-2 /data/assembly/Nextera_MP_R2_50x.fastq \
--mp1-fr -o <asm_name2>
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
* add the `--only-assembler`flag to skip correction


**Bonus exercise**: map the error-corrected paired end reads to the same assembly as the uncorrected reads (calling SNPs and indels), and check the results in the genome browser.