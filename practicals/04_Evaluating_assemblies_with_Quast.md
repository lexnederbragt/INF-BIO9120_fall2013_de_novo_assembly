#Using Quast for assembly evaluation

The Quast program can be used to generate similar metrics as the assemblathon_stat.pl script, pluss some more and some visualisations.

Program|Options|Explanation
-------|-------|-------------
Quast||Evaluating genome assemblies
|-o|name of output folder
|-R|Reference genome
|-G|File with positions of genes in the reference (see manual)
|-T|number of threads (cpu's) to use
|sequences.fasta|one or more files with assembled sequences
|-l| comma-separates list of names for the assemblies, e.g. "assembly 1", "assembly 2" (in the same order as the sequence files)
|--scaffolds|input sequences are scaffolds, not contigs. They will be split at 10 N's or more to analyse contigs ('broken' assembly)
|--est-ref-size| estimated reference genome size (when not provided)
|--gene-finding| apply GenemarkS for gene finding

See the manual for information on the output of Quast:
[http://quast.bioinf.spbau.ru/manual.html#sec3](http://quast.bioinf.spbau.ru/manual.html#sec3)

**NOTE** Quast will produce a html report file `report.html` that you can open in your browser.

####Using Quast *without* a reference genome
Note that the `--scaffold` option is not used here for simplification

```
quast.py -T 2 \
--est-ref-size 4640000 \
--gene-finding \
-o out_folder_name \
../path/to/scaffolds1.fasta \
../path/to/scaffolds2.fasta \
-l "Assembly 1, Assembly 2"
```

An additional advantage of adding the `--gene-finding` flag is that Quast will provide a `gff` file with the predicted genes. This file can be added to your genome browser session as a separate track, feel free to try this out.

####Using Quast *with* the reference genome

```
quast.py -T 2 \
-o out_folder_name \
-R /data/assembly/ref/NC_000913_K12_MG1655.fasta \
-G /data/assembly/ref/e.coli_genes.gff \
../path/to/scaffolds1.fasta \
../path/to/scaffolds2.fasta \
-l "Assembly 1, Assembly 2"
```

