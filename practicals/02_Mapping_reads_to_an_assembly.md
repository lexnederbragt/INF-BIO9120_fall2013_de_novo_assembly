#Mapping reads to an assembly and visualising the results

We will use `bwa` for mapping. this is the same program you used for the variant calling module. However, we will use a different version of `bwa`, one that is faster and more accurate with reads longer than 100bp. So instead of `bwa` we will use `bwa-0.7.5a`

####Indexing the assembly

Your new assembly now becomes the 'reference' for `bwa`. `bwa` needs an index of the sequences to make mapping go faster. For large genomes such as the human genome, this takes a long time (which is why you were given a indexed reference for the variant calling module). For the small bacterial genome we work with here this is very fast.

Move (using `cd`) to the folder with your final assembled sequenced, e.g. `contigs.fa` for velvet, or `scaffolds.fasta` for SPADES.  
Index the fasta file with:

```
bwa-0.7.5a index -a bwtsw <assembly.fasta>
```

Replace `<assembly.fasta>` with the name of your fasta file. Run `ls` to check the results, you should see a bunch of new files.


####Mapping paired end reads

Mapping the reads using `bwa-0.7.5a mem` yields SAM output. Instead of saving this output to disk, we will immediately convert it to a sorted BAM file by piping into the `samtools`program. 'Sorted' here means that the alignments of the mapped reads are in the order of the reference sequences, rather than random. Finally, we will generate an index of the sorted BAM file for faster searching later on.

First, create a new folder and `cd` into it:

```
mkdir bwa
cd bwa
```
Then do the mapping:

```
bwa-0.7.5a mem -t 2 ../<assembly.fasta> \
/data/assembly/MiSeq_Ecoli_MG1655_50x_R1.fastq \
/data/assembly/MiSeq_Ecoli_MG1655_50x_R2.fastq \
| samtools view -buS - | samtools sort - map_pe.sorted
```

Generate an index of the BAM file:

```
samtools index map_pe.sorted.bam
```

Explanation of some of the parameters:

* `../` means 'look in the folder one level up', i.e. where the fasta file is
* `-t 2`tells `bwa-0.7.5a mem` to use 2 threads (cpus)
* `-buS`tells `samtools view` that the input is in SAM format (`S`) and to output uncompressed (`u`) BAM format (`b`).
* the `-` for both `samtools` commands indicate that instead of using a file as input, the input comes from a pipe (technically, from 'standard in', or 'STDIN').
* `-map_pe.sorted` tells `samtools view` to call the outputfile `map_pe.sorted.bam`

If you would like to have a look at the alignments in the BAM file (which is in binary format), use `samtools view`again:

```
samtools view map_pe.sorted.bam |less
```

####Mapping mate pairs
Repeat the `bwa-0.7.5a mem` and `samtools` commands above, but:

* use the mate pair reads `Nextera_MP_R1_50x.fastq` and `Nextera_MP_R2_50x.fastq`
* change the output name to `map_mp.sorted`

####Plotting the insert size distribution
Since we know know where the pairs of reads map, we can obtaint he distance between them. That information is stored in the SAM/BAM outour in the 9th column, 'TLEN' (observed Template LENgth).

We will use python, and the python modules `pysam` and `matplotlib` to plot the distribution of insert sizes for a subset of the alignments. All this we will do in an IPython notebook, and interactive python web-based document with live coding and plotting of results. Don't worry, I'll *demonstrate the use of this notebook* before you try it all yourself.

* Copy the notbook file `/doc/assembly/practicals/Plot_insertsizes.ipynb` to the folder with the BAM files
* In the terminal, `cd` to the same folder
* In the terminal, write 

```
ipython notebook --pylab inline
```

* After a little bit, your webbrowser will start with a new tab labelled `IPython dashboard`, and the notbook `Plot_insertsizes` listed
* Click on the notebook name, it will open in a new tab
* Execute the cells as listed.
* For `infile`, use the name of the sorted BAM file for the mapping of the paired end or mate pair reads


####Visualising the assembly in a genome browser
For this part, we will use IGV again. 
Instead of using one of the build-in genomes, we will add the assembly as a new reference genome:

* Start the IGV program by typing `igv.sh`
* Choose `Genomes --> Load Genome from File…` (**NB** not File --> Load from File…)
* Select the fasta file with your assembly (**NB** the same file as you used for mapping the reads against!)

**Adding the mapped reads**  
Adding tracks to the browser is as simple as uploading a new file:

* Choose `File --> Load from File…`
* Choose the sorted BAM file of the paired end mapping 
* Repeat this for the BAM file of the mate pair mapping 
* Start browsing!
* Remember you can choose different sequences (contigs/scaffolds) from the drop-down menu at the top. Start browsing (one of) the longest scaffold(s)

**Question:**

* Do you see differences between some of the reads relative to the reference? What are these?
* Is coverage even? Are there gaps in the coverage, or peaks? Where?

####Adding the locations of gaps as another track
It would be convenient to be able to see the location of gaps in the browser. For this purpose run the following command (e.g., in the folder with the `bwa` results). we will use 10 bases as minimum gap length: `-m 10`

```
scaff2bed.py -i ../<assembly.fasta> -m 10 >gaps.bed
```

This will create a BED file with locations of the gaps. Add the BED file to the browser. Zoom in one gaps and look at the alignments.

**Question:**

* Are all gaps spanned by mate pairs? Tip: choose 'view as pairs' for the tracks

####Calling SNPs and INDELs
The alignments generated by `bwa` can be used to call SNPs and INDELs in the same way as you did for the cariant calling tutorial:

* Use the mapping of the paired end reads, not the mate pairs
* Run the approriate part of the variant calling tutorial, note that:
  * `picard FixMateInformation.jar` should be given the `map_pe.sorted.bam` file as input
  * Make sure you give the path to the reference fasta file, when needed, correctly
  * All `picard` commands should have the parameter `VALIDATION_STRINGENCY=SILENT` (**not** 'STRICT')
  * `picard ValidateSamFile.jar` will fail with too manyw arning, these can be ignored
  * All `gatk` commands should have the parameter `--validation_strictness SILENT` (**not** 'STRICT')
  
The commands yield two files: `snps_PE.vcf` and `indels.vcf. Add these to the browser. Have a look at these indels and the alignments around them and see whether you agree that there is a SNP or INDEL there.

###Next steps
Time to add more assemblies! See tutorials on Celera and SPADES
  