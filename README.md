# README #

This repository accompanies the paper *Non-B DNA affects polymerization speed and error rate in sequencers and living cells*

* 1.0.0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

#Sofware versions

Samtools : Version: 1.3.1 (using htslib 1.3.1)
Bedtools : Version:   v2.26.0
Bedops   : Version:  2.4.20

PacificBiosciences SMRT-Link :  smrtanalysis-2.3.0

Python 2.7.3 (default, Mar 13 2014, 11:03:55)

Naive Variant Caller, from modified version at tools-blankenberg/tools/naive_variant_caller/tools/naive_variant_caller.py

#FILE FORMATS
##.mf

The .mf file describes the interval of interest and thus overlaps with non-B DNA annotations. This file format will either contain feature-restricted intervals or full 100bp windows.
This tab delimited file contains 4 required columns:

#DATA SOURCE
##Genome in a Bottle Son - PacBio

Download files from url:
https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/hdf5/

##Genome in a Bottle Son - Illumina

Monika please complete


##non-B DNA annotation

Url for portal:
https://nonb-abcc.ncifcrf.gov/apps/Query-GFF/feature/

Enter following options:
Species: Human
Classes: NonBDAS
Search by Type: Chromosome
Chromosome: 1 to 22 (one chromosome at a time only)
Query Type: all features for the region
Features Types to Retrieve: choose all but Short Tandem Repeats



**chromosome**

**start** in 1-based format

**end** in 1-based format

**feature length** number of nucleotides in feature, note that this field is present in both feature-restricted and 100bp windows. Occassionally, feature can be longer than 100bp. In such cases, we trim it down to 100bp and restrict our analysis to trimmed windows.   

**IPD values** variable number of columns containing IPD values from ipdSummary, each column represents one nucleotide from feature or 100bp window, thus at most 100 IPD values are present
.collapsed

##.collapsed
The .collapsed file describes the final error rates in given intervals.
**chromosome**

**start** in 1-based format

**end** in 1-based format

**total rate** as fraction <0,1>

**mismatch rate** as fraction <0,1>

**insertion rate** as fraction <0,1>

**deletion rate** as fraction <0,1>

#MODULES
##Non-B DB and STR annotations.

Non-B DB annotation files downloaded at: https://nonb-abcc.ncifcrf.gov/apps/Query-GFF/feature/

Microsatellites annotation files obtained through STR-FM (see Fungtammasan, A. et al. Accurate typing of short tandem repeats from genome-wide sequencing data and its applications. Genome Res.2015.)

Microsatellites alignment and collapsing done in: *MicrosatAnnotation/bash_microsat.sh*


##Interval-Wise Testing for differences in IPDs. 

###bash_IPD.sh###
Input : unaligned cmp.h5 from PacBio sequencing project, reference genome

Output : IWT-ready files

###IWT###
Loading data in R in IWTomics format and creating subsamples for test: load_IPD_data.r

Perform IWT test and plot results: IWT_IPD.r

IPD analysis and IWT test for motifs with different length: feature_length_IPD.r


##Effect of sequence composition on IPD. 

Sequence composition obtained in: *CompoRegress/bash_comporegress.sh*

Perform sequence composition analysis: sequence_composition_IPD.r


##Experimental characterization of G-quadruplexes. 

Analyze Tm and delta-epsilon experimental results: thermostability.r


#Data formatting

generateFeatures.py


This script restricts the intervals to features only.

Input: .mf file

Output: .mf file with FeatureOnly suffix

generateEmptyTrack.py script, employing the Linux function shuf


This script generates matching controls for each motif by subsampling from motif-free windows.

Input: .mf file, optionally output directory

Output: .mf file

generateControls.sh

Script that generates 10 sets of control datasets.
Input: folder with .mf files; we will generate 10 controls for these files
Output: 10 folders with matching controls for each motif

##SMRT sequencing errors
runErrorStatistics_optimized.R, generateErrorsFullWindow.sh and generateErrors.sh scripts

##Illumina sequencing errors.
#bash_illumina.sh script and its dependencies

Now with Naive Variant Caller:

bash runIlluminaErrorDiscoveryTest.sh Joint/ bam/0 var_output collapsed_output

##Variants from human-orangutan divergence. 
bash_divergence.sh script and its dependencies

##Variants from the 1000 Genomes project.
bash_diversity.sh script and its dependencies

##Somatic mutations from The Cancer Genome Atlas.
bash_somatic.sh script and its dependencies

##Estimation of falsely reported errors due to misalignment.
misalignment_simulation.sh and its dependencies

##Impact of different aligners on calling sequencing errors.
five_aligners.sh and its dependencies

##Comparison of errors and variants between motifs and motif-free regions. 
error_plots.Rnw

#Filtering

The accuracy analysis for Illumina and PacBio removed both repeatmasked regions and regions with variants in HG002 compared to hg19.
Following high quality HG002 calls were used:
HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz

joint_filtering.sh

filter_out_repeatmasked.sh

filter_out_true_variants.sh

filter_out_repeatmasked.sh and filter_out_true_variants.sh can be run separately or together in the joint_filtering.sh. 

#Data plotting
plot_heatmap

#Data processing
Simple scripts for converting between formats used in our study:

mf_into_gff_and_bed.sh

mf2gff.sh

convert_gff_to_mf.sh

gff2mf.py

#Preparing the motif and motif-free windows

inputs : .mf files, RepeatMasked file, True Variants file

1. mf2gff.sh
2. filter_out_repeatmasked.sh (Diversity, Divergence, TCGA) & joint_filtering.sh (Genome in a Bottle - Pacbio + Illumina)
	Create sepearte folders for RepeatMasked and Joint
3. convert_gff_to_mf.sh (use on filtered files. Don't forget to remove the orginal GFFs. Only the output of step 2 will work with that script - needs gff2mf.py)
4. generateFeatures.py (remove the flanking regions of motifs in their 100bp windows. When original motifs were bigger than 100bp, only kept the 100bp at the center of the motifs)
5. generateControls.sh (generates 10 random controls for each motif. Needs generateEmptyTrack.py)