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


##Microsat annotation:

Microsatellites annotation files obtained through STR-FM (see Fungtammasan, A. et al. Accurate typing of short tandem repeats from genome-wide sequencing data and its applications. Genome Res.2015.)


##1000Genomes Project VCF files :

Download files from url:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/


##Multi-Z Alignment:

Download files from url (chromosomes 1 to 22):

http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/multiz100way/maf/


##TCGA:

Url for portal:

https://portal.gdc.cancer.gov/

Proper authentification required, data is not public access. 

##DeCode Trios:

Url for portal:

https://www.ebi.ac.uk/ena/data/view/PRJEB21300


#FILE FORMATS

##.mf

The .mf file describes the interval of interest and thus overlaps with non-B DNA annotations. This file format will either contain feature-restricted intervals or full 100bp windows.
This tab delimited file contains 4 required columns:

1. **chromosome**
2. **start** in 1-based format
3. **end** in 1-based format
4. **feature length** number of nucleotides in feature, note that this field is present in both feature-restricted and 100bp windows. Occassionally, feature can be longer than 100bp. In such cases, we trim it down to 100bp and restrict our analysis to trimmed windows.
5. **IPD values** variable number of columns containing IPD values from ipdSummary, each column represents one nucleotide from feature or 100bp window, thus at most 100 IPD values are present

##.collapsed
The .collapsed file describes the final error rates in given intervals.

1. **chromosome**
2. **start** in 1-based format
3. **end** in 1-based format
4. **total rate** as fraction <0,1>
5. **mismatch rate** as fraction <0,1>
6. **insertion rate** as fraction <0,1>
7. **deletion rate** as fraction <0,1>

#MODULES

##Steps to prepare IPD and Depth data in PacBio:

1. Download and extract the data
2. Recover *.bax.h5
3. Aligning:

	`pbalign --nproc 8 $inp ./hg19.fa ${inp}.cmp.h5 --forQuiver --metrics IPD`
	
	nproc: number of processor used - depends on your machine
	
4. Merge alignment into one file:

	`find ./ -name "*cmp.h5" -type f -exec cmph5tools.py merge --outFile out_all.cmp.h5 {} +`
	
5. Sanity checks: 

	use coverage.py (see explaination inside script) ?
	
	`cmph5tools.py summarize out_all.cmp.h5`
	
6. Deep sort alignment (optional?):

	`cmph5tools.py sort --deep out_all.cmp.h5`
	
7. Separate alignment in chromosomes:

	`cmph5tools.py select --where "(Reference == 'chr${inp}')" --outFile chr${inp}.cmp.h5 out_all.cmp.h5`
	
8. Compute IPDs:

	`ipdSummary.py $inp --reference ../hg19.fa --outfile ${inp}.idp`
	
	`cat *pickle > 52X.pickle`
	
	`python cleanIPD.py 52X.pickle > 52XIPD`
9. Compute Depth:
	`python cleanDepth.py 52X.pickle > 52XDepth`


##Steps to prepare IWT input:

1. Download non-B DNA annotation
2. Filter G4 into G4Plus and G4Minus and format:

	`grep '+' G_Quadruplex_Motifs | cut -f 1-6 > G4Plus`
	
3. Format unstranded motifs:

	`cut -f 1-6 Direct_Repeats > DirectRepeats`
	
4. Format Microsat Annotation:

	`python format.py all.per1.TRs.bed.10 Mono`
	
	`python line_up_microsats.py Mono > Mono.aligned`
	
	`python parse_repeats.py Mono.aligned`
	
	`for a in *n; do cut -f1,2,3,4,5 $a > ../GFF/$a; done;`
	
	Repeat for Di, Tri and Tetra
	
5. Prepare windows for study:	#we still have an issue of features bigger than windows overlapping with surronding windows
	
	`python ../prepare_windows.py ../GFF/ > Windows_Ready`
	
6. Sort windows:	#/!\ This step is error prone /!\ Make sure sorting was properly performed /!\
	
	`cat Windows_Ready | env LC_ALL=C sort -k 1,1d -k 2,2n > Windows_Sorted`
	
7. Collect values in windows:	#adjust collect_values_in_windows.py for F or R
	
	`python ../collect_values_in_windows.py Windows_Sorted 52XIPD > Windows_Collected_F`
	
	`python ../collect_values_in_windows.py Windows_Sorted 52XIPD > Windows_Collected_R`
	
8. Split by Featurs:

	`python ../../split_by_feature.py Windows_Collected_F`
	
	`python ../../split_by_feature.py Windows_Collected_R`
	
9. IPD data now ready for IWT. ReDo step 7 and 8 with 52XDepth for IWT on Depth.


###IWT###

1. Loading data in R in IWTomics format and creating subsamples for test: 

	`load_IPD_data.r`

2. Perform IWT test and plot results: 

	`IWT_IPD.r`

3. IPD analysis and IWT test for motifs with different length: 

	`feature_length_IPD.r`


##Effect of sequence composition on IPD.

1. Obtain coordinates of windows: 

	`python bedformatting.py Windows_Collected_F`

2. Get sequence inside each window:

	`bedtools getfasta -s -fi hg19.fa -bed Windows_Collected_F.bed > getFastaF`
	
3. Get sequence composition of each window:

	`python compo.py Windows_Collected_F.bed getFastaF > Windows_Collected_F.compo`
	
4. Split by Features:

	`python split_by_feature.py Windows_Collected_F.compo`
	
5. ReDo steps 1 to 4 with Windows_Collected_R




Sequence composition obtained in: *CompoRegress/bash_comporegress.sh*

6. Perform sequence composition analysis: 

	`sequence_composition_IPD.r`


##Experimental characterization of G-quadruplexes. 

1. Analyze Tm and delta-epsilon experimental results: 
	
	`thermostability.r`


#Data formatting

Monika: please complete

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


#Filtering

Monika: please complete

The accuracy analysis for Illumina and PacBio removed both repeatmasked regions and regions with variants in HG002 compared to hg19.
Following high quality HG002 calls were used:
HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz

joint_filtering.sh

filter_out_repeatmasked.sh

filter_out_true_variants.sh

filter_out_repeatmasked.sh and filter_out_true_variants.sh can be run separately or together in the joint_filtering.sh. 

#Data processing

Monika: please complete

Simple scripts for converting between formats used in our study:

mf_into_gff_and_bed.sh

mf2gff.sh

convert_gff_to_mf.sh

gff2mf.py

#Preparing the motif and motif-free windows

Monika: please complete

inputs : .mf files, RepeatMasked file, True Variants file

1. mf2gff.sh
2. filter_out_repeatmasked.sh (Diversity, Divergence, TCGA) & joint_filtering.sh (Genome in a Bottle - Pacbio + Illumina)
	Create sepearte folders for RepeatMasked and Joint
3. convert_gff_to_mf.sh (use on filtered files. Don't forget to remove the orginal GFFs. Only the output of step 2 will work with that script - needs gff2mf.py)
4. generateFeatures.py (remove the flanking regions of motifs in their 100bp windows. When original motifs were bigger than 100bp, only kept the 100bp at the center of the motifs)
5. generateControls.sh (generates 10 random controls for each motif. Needs generateEmptyTrack.py)



##SMRT sequencing errors

Monika: please complete

runErrorStatistics_optimized.R, generateErrorsFullWindow.sh and generateErrors.sh scripts

##Illumina sequencing errors.

Monika: please complete

#bash_illumina.sh script and its dependencies

Now with Naive Variant Caller:

bash runIlluminaErrorDiscoveryTest.sh Joint/ bam/0 var_output collapsed_output

##Variants from human-orangutan divergence.

1. Concatenate all chomosomes multiple alignment:

	`cat chr*.maf > WG.maf`
	
2. Filter by species of interest (homo - pongo - rhesus outgroup):

	`python filter_out_maf.py > WG.filtered.maf`
	
3. Parse for SNPs:

	`mafparser.py WG.filtered.maf > WGSNP.gff`
	
4. With use.galaxy.org, parse for INDELs:
	1. use .mf files for features and controls
	2. use Extract MAF Blocks on .mf files
	3. use Fetch Indels on extracted MAF blocks
	
5. Format Indels in GFF and concatenate:

	`python parse_indels.py indels.fetched > WGINDELS.gff`
	
	`cat WGSNP.gff WGINDELS.gff > Divergence.gff`
	
6. Intersect variants with features coordinates:

	`bedtools intersect -wa -wb -b Divergence.gff -a ${inp}.gff -loj > ${inp}.intersect`
	
	`python parse_intersect.py ${inp}.intersect > ${inp}.collapsed`
	
	`python reorder.py ${inp}  ${inp}.collapsed`
	
7. Compute number of variants inside features:

	`python rates.py ${inp}.intersect`

##Variants from the 1000 Genomes project.

1. Extract polymorphisms from 1000G VCF files, split by type and frequency range:

	`python filterVCFforMAF.py ALL.chr${X}.* chr${X}`
	
2. Concatenate chromosomes:

	`cat highfreqsnpchr* > highfreqsnp.intervals`
	
	`cat highfreqindelchr* > highfreqindel.intervals`
	
	`cat lowfreqsnpchr* > lowfreqsnp.intervals`
	
	`cat lowfreqindelchr* > lowfreqindel.intervals`
	
3. With use.galaxy.org, get multiple alignments (pan, gorilla, pongo, nomascus around indels:

	1. extract MAF blocks
	
	2. maf to intervals
	
4. Polarize indels:

	`python Join.py highfreqindel.maf.intervals highfreqindel.intervals`
	
	`python Polaryze.py highfreqindel.maf.intervals.joined`
	
5. Format output in GFF:

	`python Intervals_to_gff.py highfreqsnp.intervals > highfreqsnp.gff`
	
	`python Polarize_to_gff.py highfreqindel.polarized > highfreqindel.gff`
	
	`cat highfreqsnp.gff highfreqindel.gff > highfreq1kG.gff`
	
6. Intersect polymorphisms with features coordinates:

	`bedtools intersect -wa -wb -b highfreq.gff -a ${inp}.gff -loj > ${inp}.intersect`
	
	`python parse_intersect.py ${inp}.intersect > ${inp}.collapsed`
	
	`python reorder.py ${inp}  ${inp}.collapsed`

7. Compute number of polymorphisms inside features:

	`python rates.py ${inp}.intersect`
	
6. Redo steps 3 to 7 for lowfreq

##Somatic mutations from The Cancer Genome Atlas.

1. Obtain somatic mutations as described in Material and Methods and format into GFF
2. Count each segregation site once (a coordinate can happen only once)
3. Intersect variants with features coordinates:

	`bedtools intersect -wa -wb -b Divergence.gff -a ${inp}.gff -loj > ${inp}.intersect`
	
	`python parse_intersect.py ${inp}.intersect > ${inp}.collapsed`
	
	`python reorder.py ${inp}  ${inp}.collapsed`
	
4. Compute number of variants inside features:

	`python rates.py ${inp}.intersect`
	
	
##DeCode Trios:

1. Extract De Novo Mutations from DeCode VCF files and format to GFF:

	`python DeCode2GFF.py`
	
2. With use.galaxy.org, use liftOver to change coordinates from hg38 to hg19
	
3. Intersect variants with features coordinates:

	`bedtools intersect -wa -wb -b DecodeDNMhg19.gff -a ${inp}.gff -loj > ${inp}.intersect`
	
	`python parse_intersect.py ${inp}.intersect > ${inp}.collapsed`
	
	`python reorder.py ${inp}  ${inp}.collapsed`
	
4. Compute number of variants inside features:

	`python rates.py ${inp}.intersect`

##Estimation of falsely reported errors due to misalignment.

misalignment_simulation.sh and its dependencies

##Impact of different aligners on calling sequencing errors.

five_aligners.sh and its dependencies

##Comparison of errors and variants between motifs and motif-free regions.

error_plots.Rnw



#Data plotting
plot_heatmap