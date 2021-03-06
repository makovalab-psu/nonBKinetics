# README #

This repository accompanies the paper *Long-read sequencing technology indicates genome-wide effects of non-B DNA on polymerization speed and error rate *

All figures and tables were created with scripts present in that directory. Supplementary Figures/Tables deriving from main Figures/Tables are not described in details in this README as they were generated with the same scripts, usually by filtering the inputs.

Please check that your version of this repository is up to date (https://bitbucket.org/makova-lab/kinetics_wmm)

* 1.0.0
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

# Sofware versions

Samtools : Version: 1.3.1 (using htslib 1.3.1)
Bedtools : Version: v2.26.0
Bedops   : Version: 2.4.20

PacificBiosciences SMRT-Link :  smrtanalysis-2.3.0

Python 2.7.3 (default, Mar 13 2014, 11:03:55)

Naive Variant Caller, from modified version at tools-blankenberg/tools/naive_variant_caller/tools/naive_variant_caller.py


# DATA SOURCE

---

## Genome in a Bottle Son - PacBio

Download files from url:

https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/hdf5/

https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/MtSinai_blasr_bam_GRCh37/

## non-B DNA annotation

Url for portal:

https://nonb-abcc.ncifcrf.gov/apps/Query-GFF/feature/

Enter following options:

Species: Human

Classes: NonBDAS

Search by Type: Chromosome

Chromosome: 1 to 22 (one chromosome at a time only)

Query Type: all features for the region

Features Types to Retrieve: choose all but Short Tandem Repeats


## Microsatellites annotation:

Microsatellites annotation files obtained through STR-FM (see Fungtammasan, A. et al. Accurate typing of short tandem repeats from genome-wide sequencing data and its applications. Genome Res.2015.)


## 1000Genomes Project VCF files :

Download files from url:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/


## Multi-Z Alignment:

Download files from url (chromosomes 1 to 22):

http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/multiz100way/maf/


# CUSTOM FILE FORMATS

---

## .mf

The .mf file describes the interval of interest and thus overlaps with non-B DNA annotations. This file format will either contain feature-restricted intervals or full 100bp windows.
This tab delimited file contains 4 required columns:

1. **chromosome**
2. **start** in 1-based format
3. **end** in 1-based format
4. **feature length** number of nucleotides in feature, note that this field is present in both feature-restricted and 100bp windows. Occassionally, feature can be longer than 100bp. In such cases, we trim it down to 100bp and restrict our analysis to trimmed windows.
5. **IPD values** variable number of columns containing IPD values from ipdSummary, each column represents one nucleotide from feature or 100bp window, thus at most 100 IPD values are present

## .collapsed
The .collapsed file describes the final error rates in given intervals.

1. **chromosome**
2. **start** in 1-based format
3. **end** in 1-based format
4. **total rate** as fraction <0,1>
5. **mismatch rate** as fraction <0,1>
6. **insertion rate** as fraction <0,1>
7. **deletion rate** as fraction <0,1>


# FIGURE 2

---

## Steps to prepare IPD and Depth data in PacBio.

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


## Steps to prepare IWT input.

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
	
7. Collect IPD values in windows:	#adjust collect_values_in_windows.py for F or R
	
	`python ../collect_values_in_windows.py Windows_Sorted 52XIPD > Windows_Collected_F`
	
	`python ../collect_values_in_windows.py Windows_Sorted 52XIPD > Windows_Collected_R`
	
8. Split by Features: #/!\ Files created are .mf files. The extension requires to be added manually. /!\ Windows with no IPD not filtered out in this step.

	`python ../../split_by_feature.py Windows_Collected_F`
	
	`python ../../split_by_feature.py Windows_Collected_R`
	
9. Filter out windows with no IPD and create feature length files (repeat on all files generated at step 8):

	`FILES=./*`
	
	`for file in $FILES`
	
	`do`
	
	`	awk 'BEGIN {OFS=FS="\t"} {$1="chr"$1; if(gsub(/nan/,"NA")<100) print}' $file > $file"_filtered"`
	
	`	awk 'BEGIN {OFS=FS="\t"} {print $1, $2, $3, $4}' $file"_filtered" > $file"_filtered_length"`
	
	`	sed -i -r 's/(\s+)?\S+//4' $file"_filtered"`
	
	`done`

10. IPD data now ready for IWT. ReDo step 7, 8 and 9 with 52XDepth for IWT on Depth.


## Run IWT.

1. Loading data in R in IWTomics format and creating subsamples for test (see detailed comments inside R script): 

	`load_IPD_data.r`

2. Perform IWTomics test and plot results (see detailed comments inside R script): 

	`IWT_IPD.r`

3. IPD analysis and IWT test for motifs with different length (see detailed comments inside R script): 

	`feature_length_IPD.r`


## Effect of sequence composition on IPD.

1. Obtain coordinates of windows: 

	`python bedformatting.py Windows_Collected_F`

2. Get sequence inside each window:

	`bedtools getfasta -s -fi hg19.fa -bed Windows_Collected_F.bed > getFastaF`
	
3. Get sequence composition of each window:

	`python compo.py Windows_Collected_F.bed getFastaF > Windows_Collected_F.compo`
	
4. Split by Features:

	`python split_by_feature.py Windows_Collected_F.compo`
	
5. Filter out windows with no IPD and create composition files (repeat on all files generated at step 4):

	`FILES=./*`
	
	`for file in $FILES`
	
	`do`
	
	`	awk 'BEGIN {OFS=FS="\t"} {gsub(/[A]/,"",$1); gsub(/[TCG]/,"\t",$1); $3=$3+1; print $2, $3, $4, $1}' $file > $file"_composition"`
	
	`	intersectBed -wa -f 1 -a $file"_composition" -b $file"_filtered" > "../forward/"$file"_filtered_composition"`
	
	`	awk 'BEGIN {OFS=FS="\t"} {gsub(/[A]/,"",$1); gsub(/[TCG][TCG]/,"\t",$1); gsub(/[TCG]/,"\t",$1); $3=$3+1; print $2, $3, $4, $1}' $file > $file"_composition_di"`
	
	`	intersectBed -wa -f 1 -a $file"_composition_di" -b $file"_filtered" > "../forward/"$file"_filtered_composition_di"`
	
	`done`

	
6. ReDo steps 1 to 5 with Windows_Collected_R

7. Perform sequence composition analysis (see detailed comments inside R script): 

	`sequence_composition_IPD.r`



# FIGURE 3

---

## Experimental characterization of G-quadruplexes. 

1. Get the 10 most common G4 and collect IPD values. Use the following line to retrieve the sequences of interest.

	`awk '{if ($4 == "GGGTGGAGGGTGGGAGGAGGG") print $0}' Windows_Collected_F >> intermolecular`

2. Analyze delta-epsilon and Tm experimental results (see detailed comments inside R script):
	
	`thermostability.r`



# TABLE 1

---


## Preparing the motif and motif-free windows

Requires: .mf files, RepeatMasked file with coordinates of repeats in a genome, and a file with high quality calls that contrasts sequenced individual with a reference. 


1. **Generate corresponding .gff file for each .mf file.** The reason for this conversion is that .gff files are a suitable format for the intersection with other datasets.

	`mf2gff.sh folder_with_mf_files`
	This script requires the name of the folder with .mf files as a single parameter. It will generate a corresponding .gff file for each .mf file.
  
 
2. **Filter out windows that meet certain criteria.** For example, one might wish to remove windows overlapping with the known repeats in a human genome as identified by RepeatMasker. For analysis of SMRT errors, it might be also desirable to additionally remove the known variants in an analyzed human individual  (HG002). This allows one to analyze sequencing errors without confounding them with the real variants in which the sequenced individual differs from a reference genome. 

	Create separate folders for *RM* (repeatmasked), *trueVariants* (personal variants) and *joint* (both).

	`joint_filtering.sh` (SMRT errors, removes repeats and variants)
	
	This step generates new files with suffixes *interspersed_containing.gff and *interspersed_filtered.gff (overlaps with repeats or not) and suffixes *varContaining.gff and *varFiltered.gff (overlapping with variants or not), as well as prefixes joint_ (both filters applied).
	
	Following high-quality HG002 calls were used for the filtering step:
	HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz
	
	The repeat track can be downloaded from http://genome.ucsc.edu/cgi-bin/hgTables as a bed file. We used human interspersed repeats for filtering and named the file Human_interspersed.
	
	joint_filtering.sh calls both filter_out_repeatmasked.sh and filter_out_true_variants.sh; these can be also run independently.
	
	`filter_out_repeatmasked.sh` (removes only repeats, suitable for Diversity, Divergence) 
	`filter_out_true_variants.sh` (removes only variants)
	
3. **Use new coordinates in order to re-create .mf files that are filtered.** Use on filtered .gff files. Only the output of step 2 will work with the script. 
	Requires gff2mf.py (or gff2mf_regVar.py). These scripts take two input parameters: original .mf file (before filtering) and .gff file with coordinates of a new subset of windows 
	
	To automatically create filtered .mf files for folders trueVariants, RM and joint that should contain filtered .gff files, one can run:
	`convert_gff_to_mf.sh original_folder_with_mf_files`
	
	
4.	**Keep only motifs and remove the flanking regions around motifs inside their 100bp windows.**
	This step is optional, depending on whether you plan to analyze only features or full 100bp windows.
	`generateFeatures.py mf_file` will generate a new file that will only include features
	(Note that when original motifs are bigger than 100bp, only 100bp centered around the motifs will be kept; note that this is only relevant for DirectRepeats that can be very long.)

	
5.  **Generates 10 random controls for each motif.** 
	`generateControls.sh folder_wit_mf files`
	Script that generates 10 sets of control datasets for the .mf files inside the provided folder.
	Output: 10 folders with matching controls for all the motifs
	
	The script *requires a control file to be present*. Typically, a control file contains a large number of windows. From this set of windows, as many will be used as are present in a motif file. For example, if a motif file contains 2,500 windows, the subset of randomly selected 2,500 control windows will be used. This control file is expected to be named *Empty*.mf and therefore it's important to make sure that only one control file is present in a folder before running this script. Since small variation in sampling is to be expected, as many as 10 independent folders with controls will be created.
	
	
	Requires generateEmptyTrack.py (or generateEmptyTrack_regVar.py).
	Note that these scripts can be run independently by running:
	`python generateEmptyTrack.py mf_file control_file output_directory` where
	*mf_file* is the .mf file for which the control should be generated
	*control_file* is a large file with all the control windows (these would be regions of the genomes that do not contain non-B DNA)
	*output_directory* is an optional parameter for the location of an output
	
Important note: the generated controls can be either from anywhere in the genome or localized in a close proximity of a motif (such as 0.5 Mb upstream or downstream) in order to account for the regional variation. In such case, scripts in a folder controls_regVar can be used. These scripts differ slightly from the scripts described above; for example, generateEmptyTrackRegVar.py instead of generateEmptyTrack.py, or gff2mf_regVar.py instead of gff2mf.py

## Data formatting

**generateFeatures.py**
This script restricts the intervals to features only.

Input: .mf file

Output: .mf file with FeatureOnly suffix

**generateEmptyTrack.py** script
This script generates matching controls for each motif by subsampling from motif-free windows. Sampling is implemented with the Linux function shuf.

Input: .mf file, optionally output directory

Output: .mf file


## SMRT sequencing errors

Dependency: runErrorStatistics_optimized.R, reorder_pacbio.py


1. **Compute errors.**

	`generateErrors.sh`

2. **Merge and reorder.** Since many windows are run in parallel, the results need to be merged and re-ordered to match the order in the original .mf files.

	`merge_and_reorder.sh`



# FIGURE 4

---

## Variants from human-orangutan divergence.

1. Concatenate all chomosomes multiple alignment:

	`cat chr*.maf > WG.maf`
	
2. Filter by species of interest (homo - pongo - rhesus outgroup):

	`python filter_out_maf.py > WG.filtered.maf`
	
3. Parse for SNPs:

	`mafparser.py WG.filtered.maf > WGSNP.gff`
	
4. With use.galaxy.org, parse for INDELs:
	1. upload .mf files for features and controls
	2. use Extract MAF Blocks on .mf files
	3. use Fetch Indels on extracted MAF blocks
	
5. Format Indels in GFF and concatenate:

	`python parse_indels.py indels.fetched > WGINDELS.gff`
	
	`cat WGSNP.gff WGINDELS.gff > Divergence.gff`
	
6. Intersect variants with features coordinates - use .gff created in TABLE 1:

	`bedtools intersect -wa -wb -b Divergence.gff -a ${inp}.gff -loj > ${inp}.intersect`
	
	`python parse_intersect.py ${inp}.intersect > ${inp}.collapsed`
	
	`python reorder.py ${inp}  ${inp}.collapsed`
	
7. Compute number of variants inside features:

	`python rates.py ${inp}.intersect`

## Variants from the 1000 Genomes project.

1. Extract polymorphisms from 1000G VCF files, split by type and frequency range:

	`python filterVCFforMAF.py ALL.chr${X}.* chr${X}`
	
2. Concatenate chromosomes:

	`cat highfreqsnpchr* > highfreqsnp.intervals`
	
	`cat highfreqindelchr* > highfreqindel.intervals`
	
	`cat lowfreqsnpchr* > lowfreqsnp.intervals`
	
	`cat lowfreqindelchr* > lowfreqindel.intervals`
	
3. With use.galaxy.org, get multiple alignments (pan, gorilla, pongo, nomascus around indels):

	1. Upload the indel.intervals

	2. extract MAF blocks
	
	3. maf to intervals
	
4. Polarize indels:

	`python Join.py highfreqindel.maf.intervals highfreqindel.intervals`
	
	`python Polaryze.py highfreqindel.maf.intervals.joined`
	
5. Format output in GFF:

	`python Intervals_to_gff.py highfreqsnp.intervals > highfreqsnp.gff`
	
	`python Polarize_to_gff.py highfreqindel.polarized > highfreqindel.gff`
	
	`cat highfreqsnp.gff highfreqindel.gff > highfreq1kG.gff`
	
6. Intersect polymorphisms with features coordinates - use .gff created in TABLE 1:

	`bedtools intersect -wa -wb -b highfreq.gff -a ${inp}.gff -loj > ${inp}.intersect`
	
	`python parse_intersect.py ${inp}.intersect > ${inp}.collapsed`
	
	`python reorder.py ${inp}  ${inp}.collapsed`

7. Compute number of polymorphisms inside features:

	`python rates.py ${inp}.intersect`
	
6. Redo steps 3 to 7 for lowfreq (not presented in manuscript)

## Create the figure

Input needed: regression_composition_results_log.RData (from FIGURE 2)

1. Create the dataframe:

	`python MergeIPDvsErrors.py errors.collapsed .mf divergence.collapsed diversity.collapsed`
	
2. Paste with composition file:

	`paste .compo _IPD_PbError_Div_1kG > _Compo_IPD_PbError_Div_1kG`
	
3. Select compositions (repeat on all files generated at step 2):

	`FILES=./*`
	
	`for file in $FILES`
	
	`do`
	
	`	tr ' ' \\t < $file > $file"_tab"`
	
	`	awk 'BEGIN {OFS=FS="\t"} {gsub(/[A]/,"",$2); gsub(/[TCG]/,"\t",$2); gsub(/[N]/,"",$1); print $4, $5, $6, $7, $1, $2}' $file"_tab" > $file"_composition"`
	
	`	awk 'BEGIN {OFS=FS="\t"} {gsub(/[A]/,"",$3); gsub(/[TCG][TCG]/,"\t",$3); gsub(/[TCG]/,"\t",$3); gsub(/[N]/,"",$1); print $4, $5, $6, $7, $1, $3}' $file"_tab" > $file"_composition_di"`
	
	`done`

3. Compute residuals from composition regression:

	`compute_residuals.r`

4. Analyze errors vs residual IPD:

	`mismatches_vs_residualIPD.r`





# SUPPLEMENTARY NOTE 1

---


Input needed: GQuadPlus (from FIGURE 2) , Genome in a Bottle bam file (only chr21 reported in manuscript)

1. Create the flanks:

	`python flanking.py`
	
2. Find the position of read ends:

	`python readends.py`

3. Intersect and return as tab-delimited file:

	`bedtools intersect -wa -wb -a GQuadPlus2kflanks.gff -b chr21Ends.gff > GQuadPlus2kflanks_chr21Ends.intersect`

	`python Intersect2csv.py`

4. Build the histograms:

	`Rscript plothistograms.r`





# SUPPLEMENTARY NOTE 2

---


Input needed: IPD_forward.RData (from FIGURE 2)

1. Create autocorrelation plots and fast fourier transform plots using script:
	`IPD_periodicity.r`




# SUPPLEMENTARY NOTE 4 & TABLE S9

---


Input from Figure 4.

1. All quantile calculations found in Supplementary Note 4 & Table S9 can be found in:
	`quantiles.r`




# FIGURE S5

---
For this analysis, we retrieved the information about the individual molecules during PacBio sequencing. Therefore we extracted the information about MoleculeId, TemplateStart, TemplateEnd, ReadStart, ReadEnd, and Strand using cmph5tools. Next, we intersected this information with the coordinates of our motifs. This was achieved using the script

`runPasses.sh bam_folder original_mf_gff_folder`

where a first argument (bam_folder) stores the alignments to human reference in a cmp.h5 file format, as one file per each chromosome. The second argument is a directory that contains .gff files to be used for an intersection.

Next, *each pass of a molecule that intersects either the G-Quadruplexes or the control windows* is analyzed using the following R script:
`passes.Rnw`
Note that this script uses R libraries h5r and pbh5. While pbh5 library is very useful for accessing IPD information from within an R script, it might be difficult to install due to its dependencies.

Finally, statistical analysis is performed using the script
`IPD_different_passes.r`
