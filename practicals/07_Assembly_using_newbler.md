#Assembly using SPADES

The newbler program was written by the company 454 Life Science that developed the first next-generation sequencing instrument. Even though it was not written for Illumina reads, it recently added that ability. Here, we will challende newbler with our Illumina datasets.

###Using newbler

Maker a new folder called `/home/<your_username>/assembly/newbler` and `cd` into it.

To run an assembly, use the `runAssembly` command described below. Please note:

* `runAssembly` is a program that sets up a newbler assembly and starts newbler
* `-o` specifies the name of the output folder
* `-cpu 2` tells newbler to use 2 cpus
* newbler will figure out the file format itself, and even which read files (paired end and mate pairs) belong together.

```
runAssembly -o newbler_pe+mp -cpu 2 \
/data/assembly/MiSeq_Ecoli_MG1655_50x_R1.fastq \
/data/assembly/MiSeq_Ecoli_MG1655_50x_R2.fastq \
/data/assembly/Nextera_MP_R1_50x.fastq \
/data/assembly/Nextera_MP_R2_50x.fastq \
>newbler.out 2>&1
```

**NOTE** Use the `screen` command! See [https://wiki.uio.no/projects/clsi/index.php/Tip:using_screen](https://wiki.uio.no/projects/clsi/index.php/Tip:using_screen).

The output that newbler sends to the screen/output file during assembly is explained in my blog at [http://contig.wordpress.com/2010/02/09/how-newbler-works/](http://contig.wordpress.com/2010/02/09/how-newbler-works/).
It can be also found in the file called `454NewblerProgress.txt` after the assembly is done.

####Newbler output

* `454Scaffolds.fna` contains the scaffolds.
* Have a look at the `454NewblerMetrics.txt` file.  For more details of what the different parts of this file mean, check [http://contig.wordpress.com/2010/03/11/newbler-output-i-the-454newblermetrics-txt-files](http://contig.wordpress.com/2010/03/11/newbler-output-i-the-454newblermetrics-txt-files).

**Questions:**

* How does the assembly compare to the other assemblies using the (exact) same input data?
* What did newbler determine as the insert size for the paired end and mate pair libraries?

###Next steps
As for the previous assemblies, you could map reads back to the assembly, calls SNPs and Indels, run FRC and Quast, and visualise in the browser.
**NOTE:** Some tools do not accept the `.fna` file extension. Solution: make a copy of the file and call it, for example, `newbler_pe+mp.fa`
