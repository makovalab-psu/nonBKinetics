# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* This repository accompanies the paper *Non-B DNA affects polymerization speed and error rate in sequencers and living cells*
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

#FILE FORMATS
##.mf

The .mf file describes the interval of interest and thus overlaps with non-B DNA annotations. This file format will either contain feature-restricted intervals or full 100bp windows.
This tab delimited file contains 4 required columns:

**chromosome**

**start** in 1-based format

**end** in 1-based format

**feature length** number of nucleotides in feature, note that this field is present in both feature-restricted and 100bp windows. Occassionally, feature can be longer than 100bp. In such cases, we trim it down to 100bp and restrict our analysis to trimmed windows.   

**IPD values** variable number of columns containing IPD values from ipdSummary, each column represents one nucleotide from feature or 100bp window, thus at most 100 IPD values are present
.collapsed

#MODULES
##Non-B DB and STR annotations.

parse_repeats.py
##Interval-Wise Testing for differences in IPDs. 

###bash_IPD.sh###
Input : unaligned cmp.h5 from PacBio sequencing project, reference genome

Output : IWT-ready files

###IWT###
load_IPD_data.r

IWT_IPD.r

feature_length_IPD.r


##Effect of sequence composition on IPD. 

sequence_composition_IPD.r

##Experimental characterization of G-quadruplexes. 

thermostability.r

#Data formatting
generateEmptyTrack.py script, employing the Linux function shuf

##SMRT sequencing errors
runErrorStatistics_optimized.R, generateErrorsFullWindow.sh and generateErrors.sh scripts

##
Illumina sequencing errors.
bash_illumina.sh script and its dependencies

##Variants from human-orangutan divergence. 
bash_divergence.sh script and its dependencies

##Variants from the 1000 Genomes project.
bash_diversity.sh script and its dependencies

##Somatic mutations from The Cancer Genome Atlas.
bash_somatic.sh script and its dependencies

##Comparison of errors and variants between motifs and motif-free regions. 
error_plots.Rnw

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact