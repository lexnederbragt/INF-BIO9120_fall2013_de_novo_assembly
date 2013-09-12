#Assembly using velvet

###*De novo* assembly of Illumina reads using velvet

Learning points

* understand k-mers
* understand how to run Velvet
* understand the concept of k-mer coverage
* understand the effect of k-mer settings on assemblies
* understand how to use paired-end and mate pair information in Velvet
* understand the significance of some of the settings for Velvet 

####Assembling short-reads with Velvet

We will use Velvet to assemble Illumina reads on their own. Velvet uses the *de Bruijn graph* approach. 

We will assemble *E. coli K12* strain MG1655 which was sequenced on an Illumina MiSeq. The instrument read 150 bases from each direction.

We wil first use paired end reads only: 

`/data/assembly/MiSeq_Ecoli_MG1655_50x_R1.fastq`  
`/data/assembly/MiSeq_Ecoli_MG1655_50x_R2.fastq`

###Building the Velvet Index File

Velvet requires an index file to be built before the assembly takes place. We must choose a *k-* mer value for building the index. Longer *k-* mers result in a more stringent assembly, at the expense of coverage. There is no definitive value of *k* for any given project. However, there are several absolute rules:

* *k* must be less than the read length
* it should be an odd number. 

Firstly we are going to run Velvet in single-end mode, *ignoring the pairing information*. Later on we will incorporate this information.

First, 'go home':

```
cd /home/yourusername
```


or simply type

```
cd
```


Create the assembly folder:

```
mkdir assembly
cd assembly
mkdir velvet
cd velvet
```

Find a value of *k* (between 21 and 99) to start with, and record your choice in this google spreadsheet: [bit.ly/INFBIO1](bit.ly/INFBIO1). Run `velveth` to build the hash index (see below).

Program|Options|Explanation
-------|-------|-------------
velveth||Build the Velvet index file|
|foldername|use this name for the results folder
|value_of_k|use k-mers of this size
|-short|short reads (as opposed to long, Sanger-like reads)
|-separate|read1 and read2 are in separate files
|-fastq|read type is fastq


```
velveth <asm_name> <value_of_k> \  
-short -separate -fastq \  
/data/assembly/MiSeq_Ecoli_MG1655_50x_R1.fastq \  
/data/assembly/MiSeq_Ecoli_MG1655_50x_R2.fastq  
```

After it has finished, look in the folder `asm_name`. You should see the following files:

>`Log`  
>`Roadmaps`  
>`Sequences`  

`Log` is a useful file, this is a useful reminder of what commands you typed to get this assembly result, useful for reproducing results later on. Sequences contains the sequences we put in, and `Roadmaps` contains the index you just created.

Now we will run the assembly with default parameters:

```
velvetg <asm_name>
```

Velvet will end with a text like this:

`Final graph has ... nodes and n50 of ..., max ..., total ..., using .../... reads`

The number of nodes represents the number of nodes in the graph, which (more or less) is the number of contigs. Velvet reports its N50 (as well as everything else) in 'kmer' space. The conversion to 'basespace' is as simple as adding k-1 to the reported length.

Look again at the folder `asm_name`, you should see the following extra files:

`contigs.fa`  
`Graph`  
`LastGraph`  
`PreGraph`  
`stats.txt`

The important files are:

`contigs.fa` - the assembly itself  
`Graph` - a textual representation of the contig graph  
`stats.txt` - a file containing statistics on each contig

**Questions**

* What k-mer did you use?
* What is the N50 of an assembly with 7 contigs of sizes: 20, 9, 9, 6, 3, 2 and 1 long?
* What is the N50 of *your* assembly?
* What is the size of the largest contig?
* How many contigs are there in the `contigs.fa` file? Use `grep -c NODE contigs.fa`


Log your results in this google spreadsheet: `bit.ly/INFBIO1`


**We will discuss the results together and determine *the optimal* k-mer for this dataset.**

**Advanced tip:** You can also use VelvetOptimiser to automate this process of selecting appropriate *k*-mer values. VelvetOptimizer is included with the Velvet installation.

Now run `velveth` and `velvetg` for the kmer size determined by the whole class. Use this kmer from now on!

####Estimating and setting exp_cov

Much better assemblies are produced if Velvet understands the expected coverage for unique regions of your genome. This allows it to try and resolve repeats. The command `velvet-estimate-exp-cov.pl` is supplied with Velvet and will plot a histogram of k-mer frequency for each node in the graph, listing k-mer frequency, and the number of count of nodes with that frequency

`velvet-estimate-exp_cov.pl <asm_name>/stats.txt`

The output shows:

* k-mer coverage
* count of contigs (nodes in the *de Bruijn* graph) with that coverage
* series of *'s making up a histogram 

The *peak value* in this histogram can be used as a guide to the best k-mer value for `exp_cov`.

**Questions:**

* What k-mer frequency is the most frequent?
* Why?
* What do you think is the approximate expected k-mer coverage for your assembly?
* Convert this value from k-mer coverage into genome coverage (average number of times the genome was sequenced).

The formula is:

```
   Ck * L   
-------------  = C  
 (L - k + 1)
```

Where Ck = k-mer coverage, L = read length, k = k-mer size for your assembly

Now run velvet again, supplying the value for `exp_cov` (k-mer coverage, *not* genome coverage) corresponding to your answer:

```
velvetg <asm_name> -exp_cov <peak_k_mer_coverage>
```
**Question:**

* What improvements do you see in the assembly by setting a value for `exp_cov`?

####Setting *cov_cutoff*

You can also clean up the graph by removing low-frequency nodes from the *de Bruijn* graph using the `cov_cutoff` parameter. This will often result in better assemblies, but setting the cut-off too high will also result in losing bases. Using the histogram from previously, estimate a good value for `cov_cutoff`.

```
velvetg <asm_name> -exp_cov <your_value> -cov_cutoff <your_value>  
```

Try some different values for `cov_cutoff`, keeping `exp_cov` the same and record your assembly results.

####Asking velvet to determine the parameters

You can also ask Velvet to predict the values for you:

```
velvetg <asm_name> -exp_cov auto -cov_cutoff auto
```

**Questions:**

* What values of *exp_cov* and *cov_cutoff* did Velvet choose?
* Check the output to the screen. Is this assembly better than your best one?

####Incorporating paired-end information

Paired end information contributes additional information to the assembly, allowing contigs to be scaffolded. We will first re-index your reads telling Velvet to use paired-end information, by using `-shortPaired` instead of `-short` for `velveth`. Then, re-run velvetg using the best value of `k`, `exp_cov` and `cov_cutoff` from the previous step.

**!!! IMPORTANT Pick a new name for your assembly !!!**


```
velveth <asm_name2> <value_of_k> \  
-shortPaired -fastq -separate \  
/data/assembly/MiSeq_Ecoli_MG1655_50x_R1.fastq \  
/data/assembly/MiSeq_Ecoli_MG1655_50x_R2.fastq

velvetg <asm_name2> -exp_cov <value_of_exp_cov> \  
-cov_cutoff <value_of_cov_cutoff>  
```

**Questions:**

* How does doing this affect the assembly?
* what does velvet say about the insert size of the paired end library?

####Scaffold and contig metrics

The sequences in the `contigs.fa` file are actually scaffolds. Use the `assemblathon_stats.pl` script to generate metrics for this, and all following assemblies.


**The assemblathon stats script**  
The assemblathon [www.assemblathon.org](www.assemblathon.org) used a perl script to obtain standardized metrics for the assemblies that were submitted. Here we use (a slightly modified version of) this script. It takes the size of the genome, and one sequence fasta file as input. The script breaks the sequences into contigs when there are 20 or more N’s, and reports all sorts of metrics.


Program|Options|Explanation
-------|-------|-----------
assemblathon_stats.pl| |Provide basic assembly metrics
|-size|size (in Mbp, million basepairs) of target genome (optional)
|seq.fasta|fasta file of contigs or scaffolds to report on

Example, for a 3.2 Mbp genome:

```
assemblathon_stats.pl -s 3.2 scaffolds.fasta
```

OR, save the output to a file with

```
assemblathon_stats.pl -s 3.2 scaffolds.fasta > metrics.txt
```

Here, `>` (redirect) symbol used to ‘redirect’ what is written to the screen to a file.

**For this exercise**, use the the known length for this strain, 4.6 Mbp, for the genome size

Some of the metrics the script reports are:

* N50 is based on the total assembly size
* NG50 is based on the estimated/known genome size
* L50 (LG50) count: number of scaffolds/contigs at least N50 (NG50) bases

**Questions**

* How much of the estimated genome size is covered in the scaffolds
* how many gap bases ('N') are left in the scaffolds


####Looking for repeats

Have a look for contigs which are long and have a much higher coverage than the average for your genome. One tedious way to do this is to look into the `contigs.fa` file (with `less`). You will see the name of the contig ('NODE'), it's length and the kmer coverage. However, trying to find long contigs with high coverage this way is not very efficient.  

A faster was is to use the `stats.txt` file. The full description of this file is in the Velvet Manual, at [http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf](http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf). 

Relevant columns are:

1) ID --> sequence ID, same as 'NODE' number in the `contigs.fa` file  
2) lgth --> sequence 'length' (BUT see the velvet manual)  
6) short1_cov --> kmer coverage (column 6)  


Knowing this, we can use the `awk` command to select lines for contigs at least 1kb, with k-mer coverage greater than 60:

`awk '($2>=1000 && $6>=60)' stats.txt`

`awk` is an amazing program for tabular data. In this case, we ask it to check that column 2 ($2, the length) is at least 1000 and column 6 ($6, coverage) at least 60. If this is the case, awk will print the entire line. See [http://bit.ly/QjbWr7](http://bit.ly/QjbWr7) for more information on awk.

Find the contig with the highest coverage in the `contigs.fa` file. Perform a BLAST search using NCBI. 

**Question:**

* What is it?
* Is this surprising? Why, or why not?

####The effect of mate pair library reads

Long-range "mate-pair" libraries can also dramatically improve an assembly by scaffolding contigs. Typical sizes for Illumina are 2kb and 6kb, although any size is theoretically possible. You can supply a second library to Velvet. However, it is important that files are reverse-complemented first as Velvet expects a specific orientation. We have supplied a 3kb mate-pair library in the correct orientation.

**!!! IMPORTANT Pick a new name for your assembly !!!**

We will use `-shortPaired` for the paired end library reads as before, and add `-shortPaired2` for the mate pairs:

```
velveth <asm_name3> value_of_k \  
-shortPaired -separate -fastq \  
/data/assembly/MiSeq_Ecoli_MG1655_50x_R1.fastq \  
/data/assembly/MiSeq_Ecoli_MG1655_50x_R2.fastq \  
-shortPaired2 -separate -fastq \  
/data/assembly/Nextera_MP_R1_50x.fastq \  
/data/assembly/Nextera_MP_R2_50x.fastq  
```

We use auto values for velvetg because the addition of new reads will change the genome coverage:

```
velvetg <asm_name3> -cov_cutoff auto -exp_cov auto
```

**Questions:**

* What is the N50 of this assembly?
* How many scaffolds?
* How many bases are in gaps?
* What did velvet estimate for the insert length of the paired-end reads, and for the standard deviation? Use the last mention of this in the velvet output.
* And for the mate-pair library?


Make a copy of the contigs file and call it `velvet_pe+MP.fa`

Mate-pair data can contain significant paired-end contamination which generates misassemblies. To account for this, let Velvet know about the problem with the following flag:

```
velvetg <asm_name3> -cov_cutoff auto \  
-exp_cov auto -shortMatePaired2 yes
```

**Questions:**

* What is the N50 of this assembly?
* How many scaffolds?
* How many bases are in gaps?

Make a copy of the contigs file and call it `velvet_pe+mp2`


####Skipping the paired end reads
As both the mate pairs and the paired end *reads* are of the same length, and provide the same coverage, it could be intersting to try an assembly of the mate pair reads only. The read sequences would still be used to build the contigs, and the mate pair information to build scaffolds.

**!!! IMPORTANT Pick a new name for your assembly !!!**

The last assembly for this part then becomes:

```
velveth <asm_name4> <value_of_k> \  
-shortPaired -separate -fastq \  
/data/assembly/Nextera_MP_R1_50x.fastq \  
/data/assembly/Nextera_MP_R2_50x.fastq  

velvetg <asm_name4> -cov_cutoff auto \  
-exp_cov auto -shortMatePaired2 yes
```

**Questions:**

* What is the N50 of this assembly?
* How many scaffolds?
* How many bases are in gaps?
* How does this assembly compare to the previous ones?

Make a copy of the contigs file and call it `velvet_mp_only`
