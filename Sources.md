
##Sources of programs, scripts, datafiles etc

###Datafiles

* Miseq 2 x 150 paired end reads
  * from [http://www.illumina.com/science/data_library.ilmn](http://www.illumina.com/science/data_library.ilmn)
  * random subsampling using seqtk [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)
* Nextera mate pair reads
  * from Illumina basespace [https://basespace.illumina.com/‎](https://basespace.illumina.com/‎), look for "Nextera Mate Pair (E. Coli)" [https://basespace.illumina.com/project/294296/Nextera-Mate-Pair-E-Coli](https://basespace.illumina.com/project/294296/Nextera-Mate-Pair-E-Coli)
* PacBio reads
  * from [http://www.pacificbiosciences.com/devnet/files/software/hgap/HGAp_BAS_H5_DATA/HGAp_BAS_H5_DATA/ecoli_MG1655](http://www.pacificbiosciences.com/devnet/files/software/hgap/HGAp_BAS_H5_DATA/HGAp_BAS_H5_DATA/ecoli_MG1655)
  * pre-assembled using smrtanalysis 2.0.1 from PacBio, see [http://pacbiodevnet.com/](http://pacbiodevnet.com/)

###Read QC

* FastQC v0.10.1 from [http://www.bioinformatics.babraham.ac.uk/projects/fastqc/](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

###Assembly programs

* Velvet version 1.2.09 from [http://www.ebi.ac.uk/~zerbino/velvet/](http://www.ebi.ac.uk/~zerbino/velvet/)
* SPAdes genome assembler v.2.5.1 from [http://bioinf.spbau.ru/spades](http://bioinf.spbau.ru/spades)
* Celera Assembler from [http://sourceforge.net/apps/mediawiki/wgs-assembler/](http://sourceforge.net/apps/mediawiki/wgs-assembler/). We used the version from CVS tip ([http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=Check_out_and_Compile](http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=Check_out_and_Compile)) with IDs 4371 (2013-08-01) for all parts, except 4393 (2013-08-24)for AS_CGW_main.C 
* Newbler was version 2.8, and can be requested from [http://454.com/contact-us/software-request.asp](http://454.com/contact-us/software-request.asp)

###Other programs

* bwa version 0.7.5a-r405 [http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/)
* samtools version: 0.1.18 from [http://samtools.sourceforge.net/)(http://samtools.sourceforge.net/])
* IPython and the IPython notebook from [http://ipython.org/](http://ipython.org/)
* IGV version 2.3 from [http://www.broadinstitute.org/igv/](http://www.broadinstitute.org/igv/)
* FRC_align (FRCbam) version 1.0 from [https://github.com/vezzi/FRC_align](https://github.com/vezzi/FRC_align)
* REAPR version: 1.0.16 from [http://www.sanger.ac.uk/resources/software/reapr/](http://www.sanger.ac.uk/resources/software/reapr/)
* QUAST from [http://bioinf.spbau.ru/quast](http://bioinf.spbau.ru/quast)

###Scripts

* velvet-estimate-exp_cov.pl is included in the velvet distribution
* assemblathon_stats.pl See [https://github.com/lexnederbragt/sequencetools](https://github.com/lexnederbragt/sequencetools). Modified from [https://github.com/ucdavis-bioinformatics/assemblathon2-analysis](https://github.com/ucdavis-bioinformatics/assemblathon2-analysis)
* scaff2bed.py can be found as scaffoldgap2bed.py at [https://github.com/lexnederbragt/sequencetools](https://github.com/lexnederbragt/sequencetools)
