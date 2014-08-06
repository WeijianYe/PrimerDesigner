PrimerDesigner
==============

##Overview
A quick tool for automating primer3 design and primer specificity check workflow. This repository contains several perl scripts(PrimerDesigner.pl, fastaDeal.pl, ConvertPrimer2Fq.pl, PrimerStat.pl, PrimerFinalTable.pl), you have to put them all in the identical directory and apply the option "--thirdD" to specify the source.

##Software Requirements
The following dependencies have been confirmed to work for running the 'PrimerDesigner.pl' pipeline, make sure that you use the correct version.
* **[primer3](http://primer3.sourceforge.net/releases.php)**: v2.3.5+
* **[Burrows-Wheeler Aligner](http://sourceforge.net/projects/bio-bwa/files/£©**: bwa-0.7.10+
* **[SAMtools](http://sourceforge.net/projects/samtools/files/)**: v0.1.16+

##USAGE:
There are two ways for the primer design, they are fairly indenpendent.

#For region input
	perl PrimerDesigner.pl -region region.txt --refD refDir --bwa bwa --primer primer3 --thirdD demoDir --out outdir primerDesign.I primerDesign.O
Note:
region.file: Name, chrN (do not use Chr/chromosome, please.) start end, tab separated.
refDir: you have to specify the directory contained your reference files ( must be splited by chromosomes and have a wholegenome fa file named hg19.fa or hg18.fa )

#For fa format input
	perl PrimerDesigner.pl -fa region.fa -refD refDir --bwa bwa --primer primer3 --thirdD demoDir --out outdir primerDesign.I primerDesign.O
Note:
For the region.fa, the title of each reads have to be like "chr_pos1_pos2" (eg:>chr10_43677329_43677329), which means the position of each variants.In the sequence body, you can use "[]" or "<>" to require the included regions or excluded regions. (eg: ATCT<CC>C[CTC]ATTATG, PrimerDesigner will design the primers flank the "CTC" region and forbids the primers in the central "CC".

Parameters
========
Some parameters should be noticed.
--rep: You can use a mispriming library, fa-format, primer3 should use to screen for interspersed repeats or for other sequence to avoid as a location for primers. Suggestion library below,
`HUMAN (contains microsatellites)<http://bioinfo.ut.ee/cgi-bin/primer3-0.4.0/cat_humrep_and_simple.cgi>`_
RODENT_AND_SIMPLE (contains microsatellites)
RODENT (does not contain microsatellites)
DROSOPHILA






