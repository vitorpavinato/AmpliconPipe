
## **PypeAmplicon v1.0: Python pipeline for analysis of amplicon data**

**Authors:** 

Saranga Wijeratne & Vitor A. C. Pavinato


**CITATION**

Wijeratne S & Pavinato VAC. PypeAmplicon v1.0: Python pipeline for analysis of amplicon data.
[![DOI](https://zenodo.org/badge/114859807.svg)](https://zenodo.org/badge/latestdoi/114859807)


**INSTALLATION**

You should install the softwares listed below and the required python packages (see main script) 
in order to run the pipeline. For the softwares installation, follow their instructions provided. 
For the python packages, install them either by pip or anaconda depending how you managed your 
python installation you have.

**Third part softwares:**
- [usearch](https://www.drive5.com/usearch/);
- [bbmap and bbsplit](https://sourceforge.net/projects/bbmap/);
- [seqtk](https://github.com/lh3/seqtk)
- [picard](https://broadinstitute.github.io/picard/);
- [freebayes](https://github.com/ekg/freebayes).


**USAGE**
```
python pypeamplicon.py pypeamplicon.cfg
```

**CONFIGURATION FILE**

In order to run the pipeline you should speficy the configuration file. This file contains paths 
and necessary information to run the pipeline. Below are detailed information in how to set each 
configuration section.

__[SOFTWARES]__

You should install and specify their path.

__[REF]__

You should specify the path for the file containing the FASTA sequence you are using as a reference 
sequence. It can be either the reference genome/scaffold OR the fragment homologous to the amplicon.

__[FQDIR]__

You should specify where the .fastq files were stored.

__[WORKINGDIR]__

You should specify your working directory. We suggest the following organization:
- <Project_folder> as the working directory;
	- <reference_file> Where you store your reference file;
	- <raw_files> Where you store your raw fastq files;

__[COMMANSUFFIX]__

In this section you can specify the type of the sequencing data you have
- **PE:** for paired-end data;
- **SE:** for single-end data;

For **PE data** you should provide the sulfix for the first-read files (R1) and for the second-read 
files (R2). However, in order to work, your raw fastq files should have this sulfix at the end of the 
file name (right before .fastq).

__[CLUSTERFILTER]__

Clustering filter section allows you to control the behavior of two parameters that define the filtering 
of your raw data:
1. filtercutoff: the maximum percentage of divergence allowed between reads within a cluster of reads;
2. idcutoff: the minimum proportion of similarities that clusters of reads should have to be considered 
the same cluster;


**POSTPROCESSING SCRIPTS**

At the final stage of the pipeline, the variant calling is performed using freebayes. However if, for some 
reason, you would like to have more control at this stage, we provide a sh script that runs freebayes. The 
parameters set in this script were tested in three independet datasets. You can modify it to fit your needs.

Script __freebayes-pypeamplicon.sh__

```
sh freebayes-pypeamplicon.sh <input.bam> <output.vcf> <populations-file>
```

The **population-file** is a tab limited file containing in each row:

__<sample_id>	tab <pop_id>__

If you are working with amplicon sequencing and have a list of targeted SNPs (already tested in other experiments), 
you probably are interested in process the raw vcf of the amplicon sequencing experiment to remove the non-targeted snps. 
We have two other sh scripts that allow you to:
1. Keep only the targeted SNPs in the filtered .vcf file;
2. Keep only the non-target SNPs in another .vcf file.

To keep only the targeted SNPs

```
sh vcfintersect-targets.sh <input.vcf> <output.vcf>
```

To keep only the non-targeted SNPs
```
sh vcfintersect-non-targets.sh <input.vcf> <output.vcf>
```

In both scripts you need to specify a BED file that contains the position of the targeted SNPs.
