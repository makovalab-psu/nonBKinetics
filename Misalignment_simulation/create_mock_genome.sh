#!/bin/bash

mockName=$1     # e.g. "plum"
featuresDir=$2  # e.g. "/nfs/brubeck.bx.psu.edu/scratch6/monika/true_variants/RM"
readLength=$3   # e.g. 100

featureList="APhasedRepeats DirectRepeats GQuadPlus GQuadMinus InvertedRepeats MirrorRepeats ZDNAMotifs"

# Create a mock genome consisting of motifs with 99-bp flanks, separated by
# runs of 100 Ns

# some prep work ...

twoBitInfo ${hg19}/seq/hg19.2bit /dev/stdout \
  | awk '{ print $1,"0",$2 }' \
  > temp.${mockName}.hg19.clippers

separatorLen=${readLength}
separator=`echo ${readLength} \
             | awk '{ sep="";
                      for (i=1;i<=$0;i++) sep=sep"N";
                      print sep; }'`

# main loop to create the mock genome

:> targets/${mockName}.fa
:> temp.${mockName}.target.map.xxx
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
        | clip_intervals --sorted temp.${mockName}.hg19.clippers \
        | awk '{ print $1":"$2"-"$3 }' \
        > temp.${mockName}.${feature}.ucsc
      twoBitToFa ${hg19}/seq/hg19.2bit /dev/stdout \
            -noMask -seqList=temp.${mockName}.${feature}.ucsc \
        > temp.${mockName}.${feature}.fa.xxx
      #
      fakeStart=0
      :> temp.${mockName}.${feature}.map.txt
      echo ">${fakeChr} ${feature}+$((readLength-1))bp_flanks" \
        > temp.${mockName}.${feature}.fa
      cat temp.${mockName}.${feature}.fa.xxx \
        | fasta_to_one_line --nosep --noarrow \
        | tr ":-" "  " \
        | while read chr start end nucs ; do
            #echo "=== ${feature} ${chr} ${start} ${end} to ${fakeChr} ${fakeStart} ==="
            echo ${fakeChr} ${fakeStart} ${chr} ${start} \
              >> temp.${mockName}.${feature}.map.txt
            echo ${nucs}      >> temp.${mockName}.${feature}.fa
            echo ${separator} >> temp.${mockName}.${feature}.fa
            nucsLen=`echo ${nucs} | awk '{ print length($0) }'`
            fakeStart=$((fakeStart+nucsLen+separatorLen))
            done
      #
      cat temp.${mockName}.${feature}.fa \
        >> targets/${mockName}.fa
      #
      cat temp.${mockName}.${feature}.map.txt \
        >> temp.${mockName}.target.map.xxx
      #
      rm temp.${mockName}.${feature}.ucsc
      rm temp.${mockName}.${feature}.fa.xxx
      rm temp.${mockName}.${feature}.fa
      rm temp.${mockName}.${feature}.map.txt
      done

# finish the interval mapping files

cat temp.${mockName}.target.map.xxx \
  > targets/${mockName}.map.txt
rm temp.${mockName}.target.map.xxx  

cat targets/${mockName}.map.txt \
  | awk '{ print $3,$4,$1,$2 }' \
  > targets/${mockName}.reverse_map.txt

# build the index(es) ...

bwa index -a bwtsw targets/${mockName}.fa

samtools faidx targets/${mockName}.fa

# convert to 2bit so I can do easy lookups ...

faToTwoBit targets/${mockName}.fa targets/${mockName}.2bit

# make a list of sam @SQ header records ...

cat targets/${mockName}.map.txt \
  | awk '{ print $3 }' \
  | sort -u \
  > temp.${mockName}.hg19.chroms

twoBitInfo ${hg19}/seq/hg19.2bit /dev/stdout \
  | select_lines_by_name temp.${mockName}.hg19.chroms \
  | awk '{ print "@SQ SN:"$1" LN:"$2 }' \
  | tr " " "\t" \
  > targets/${mockName}.map.sq_headers

