#!/bin/bash

featuresDir=$1  # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/RM"
readLength=$2   # e.g. 100

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Create a mock genome consisting of motifs with 99-bp flanks, separated by
# runs of 100 Ns

# some prep work ...

twoBitInfo ${hg19}/seq/hg19.2bit /dev/stdout \
  | awk '{ print $1,"0",$2 }' \
  > temp.features5.hg19.clippers

separatorLen=${readLength}
separator=`echo ${readLength} \
             | awk '{ sep="";
                      for (i=1;i<=$0;i++) sep=sep"N";
                      print sep; }'`

# main loop to create the mock genome

:> targets/features5.fa
:> temp.features5.target.map.xxx
fakeChrNum=0
echo ${featureList} | tr " " "\n" \
  | while read featureSpec ; do
      fakeChrNum=$((fakeChrNum+1))
      fakeChr=chr${fakeChrNum}fake
      #
      featureFile=${featuresDir}/${feature}_filtered.dat
      numIntervals=`wc -l ${featureFile} | awk '{ print $1 }'`
      echo "=== ${fakeChr} ${feature}  ${numIntervals} intervals ==="
      #
      cat ${featureFile} \
        | dilate_intervals $((readLength-1)) \
        | clip_intervals --sorted temp.features5.hg19.clippers \
        | awk '{ print $1":"$2"-"$3 }' \
        > temp.features5.${feature}.ucsc
      twoBitToFa ${hg19}/seq/hg19.2bit /dev/stdout \
            -noMask -seqList=temp.features5.${feature}.ucsc \
        > temp.features5.${feature}.fa.xxx
      #
      fakeStart=0
      :> temp.features5.${feature}.map.txt
      echo ">${fakeChr} ${feature}+$((readLength-1))bp_flanks" \
        > temp.features5.${feature}.fa
      cat temp.features5.${feature}.fa.xxx \
        | fasta_to_one_line --nosep --noarrow \
        | tr ":-" "  " \
        | while read chr start end nucs ; do
            #echo "=== ${feature} ${chr} ${start} ${end} to ${fakeChr} ${fakeStart} ==="
            echo ${fakeChr} ${fakeStart} ${chr} ${start} \
              >> temp.features5.${feature}.map.txt
            echo ${nucs}      >> temp.features5.${feature}.fa
            echo ${separator} >> temp.features5.${feature}.fa
            nucsLen=`echo ${nucs} | awk '{ print length($0) }'`
            fakeStart=$((fakeStart+nucsLen+separatorLen))
            done
      #
      cat temp.features5.${feature}.fa \
        | fasta_set_line_length 100 \
        >> targets/features5.fa
      #
      cat temp.features5.${feature}.map.txt \
        >> temp.features5.target.map.xxx
      #
      rm temp.features5.${feature}.ucsc
      rm temp.features5.${feature}.fa.xxx
      rm temp.features5.${feature}.fa
      rm temp.features5.${feature}.map.txt
      done

# finish the interval mapping files

cat temp.features5.target.map.xxx \
  > targets/features5.map.txt
rm temp.features5.target.map.xxx  

cat targets/features5.map.txt \
  | awk '{ print $3,$4,$1,$2 }' \
  > targets/features5.reverse_map.txt

# build the index(es) ...

bwa index -a bwtsw targets/features5.fa

samtools faidx targets/features5.fa

# convert to 2bit so I can do easy lookups ...

faToTwoBit targets/features5.fa targets/features5.2bit

# make a list of sam @SQ header records ...

cat targets/features5.map.txt \
  | awk '{ print $3 }' \
  | sort -u \
  > temp.features5.hg19.chroms

twoBitInfo ${hg19}/seq/hg19.2bit /dev/stdout \
  | select_lines_by_name temp.features5.hg19.chroms \
  | awk '{ print "@SQ SN:"$1" LN:"$2 }' \
  | tr " " "\t" \
  > targets/features5.map.sq_headers

