#Evaluating assemblies using FRCbam

Several tools exist for evaluating assemblies *in the absence of a reference genome*, by using the results of mapping reads back to the assembly. For this module, we will use FRCbam.

From the paper:

> […] based on the principle that the assembly precision can be predicted by identifying on each contig a set of suspicious regions (i.e., *features*): contigs are then sorted from the longest to the shortest, and for each feature threshold t only the longest contigs whose total sum of features is less than t are used to compute the genome coverage

These calculations are then used to provide 

>a new metric, Feature Response Curve (FRCurve), capable of capturing the trade-off between contigs’ contiguity and correctness.

FRCurve is very useful for comparing assemblies. For example, the assemblies produced by two different programs on the same dataset. Or, the same program with two different datasets (with and without mate pairs, for example).

To start this part of the tutorial, we will compare the 'best' velvet assembly with the 'best' SPADES assembly of both the paired end and mate pair data. After that, feel free to run FRCBam on any other assembly (once you have the reads mapped, it runs very fast).

The input for FRCBam is the bam files from the mapping of the reads. If you have not already done so, perform these mappings as described in the tutorial 'Mapping_reads_to_an_assembly'.

FRCBam calls feautures based on the observed distances between paired reads (and mates), and local read coverage. It will report for example regions where paired end distances are too short, or too long, or regions with very low of high coverage.

####Running FRC

First, create a folder called `/home<your_username>/assembly/FRC/` and move into that folder.

For each assembly you want to include:
* Check that there are BAM files of paired end and mate pairs mapped to the assembly (probably in a folder called `bwa`)
* For the command below, adjust the path to and name of the BAM files as needed
* Also make sure the name for <asm_name> is short, clear, and unique (otherwise any existing files from the last time the same name was used will be overwritten)
* run FRCBam as follows:

```
FRC --pe-sam <path/to/asm_folder>/bwa/map_pe.sorted.bam \
--pe-min-insert 250 --pe-max-insert 350 \
--mp-sam <path/to/asm_folder>/bwa/map_mp.sorted.bam \
--mp-min-insert 1650 --mp-max-insert 5000 \
--CEstats-PE-min -2 --CEstats-PE-max 3 \
--CEstats-MP-min -5 --CEstats-MP-max 5 \
--genome-size 4630000 --output <asm_name> > <asm_name_FRC.out>
```

* Inspect the output, pay attention to the filenames
* make a note of the total number of feautures
* Note that there is a file called `…Features.gff` that can be added to the browser, feel free to try this out.

**NOTE**

* The minimum and maximum insert sizes result from a few tries I did in preparation for the course. They are really library dependent. See also the results of the plotting of the insert size distributions
* The CEstats limits result from an exchange with the author of the program who was so kind as to advise me on which setteings were best for these libraries. 'CE stats' are related to the compression/expansion features (too wide or too narrow distribution of insert sizes)
 
The output of the `FRC` command is a set of files whose name all start with wat you entered for `--output`

####Plotting the FRC results

* While still in `home<your_username>/assembly/FRC`, copy the IPython notebook called `FRC_plotting.ipynb` from `/doc/assembly/practicals/` into the FRCplot
* Start the IPython notebook browser, open the notebook and follow the instructions to plot the FRcurves:

```
ipython notebook --pylab inline
```

NOTE:
The different types of feautures are explained on page 8 of the supplementary information of this article: [http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0052210](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0052210)
