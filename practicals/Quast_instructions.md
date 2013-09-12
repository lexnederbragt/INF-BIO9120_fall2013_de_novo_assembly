#Using Quast for assembly evaluation

The Quast program can be used to generate similar metrics as the assemblathon_stat.pl script, pluss some more and some visualisations.

Program|Options|Explanation
-------|-------|-------------
Quast||Evaluating genome assemblies
|-o|name of output folder
|-R|Reference genome
|-G|File with positions of genes in the reference (see manual)
|-t|number of threads (cpu's) to use
|sequences.fasta|one or more files with assembled sequences
|-l| comma-separates list of names for the assemblies, e.g. "assembly 1", "assembly 2" (in the same order as the sequence files)
|--scaffolds|input sequences are scaffolds, not contigs. They will be split at 10 N's or more to analyse contigs
|--est-ref-size| estimated reference genome size (when not provided)

See the manual for information on the output of Quast:
[http://quast.bioinf.spbau.ru/manual.html#sec3](http://quast.bioinf.spbau.ru/manual.html#sec3)

**NOTE** Quast will produce a html report file `report.html` that you can open in your browser.

####Using Quast *without* a reference genome

```
quast.py -t 2 \
--scaffolds \
-o out_folder_name \
../path/to/scaffolds1.fasta \
../path/to/scaffolds2.fasta \
-l "Assembly 1, Assembly 2"
```

####Using Quast *with* a reference genome

```
quast.py -t 2 \
--scaffolds \
-o out_folder_name \
-R /data/assembly/ref/NC_000913_K12_MG1655.fasta \
-G /data/assembly/ref/e.coli_genes.gff \
../path/to/scaffolds1.fasta \
../path/to/scaffolds2.fasta \
-l "Assembly 1, Assembly 2"
```