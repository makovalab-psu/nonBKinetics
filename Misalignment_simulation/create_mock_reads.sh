#!/bin/bash

mockName=$1         # e.g. "plum"
sequencingDepth=$2  # e.g. 60
readLength=$3       # e.g. 100
subsPerMillion=$4   # e.g. 2000
indelsPerMillion=$5 # e.g. 10

# Sample reads from the mock genome

# prep-work; determine how many possible read-start positions there are on each
# chromosomes; note that each mock chromosome contains one feature

cat targets/${mockName}.fa \
  | fasta_to_one_line --nosep --noarrow \
  | split_by_nth_column --column=1 temp.{name}.zzz

ls temp.*.zzz \
  | sed "s/temp\.//" \
  | sed "s/\.zzz//" \
  | while read chrom ; do
      echo "=== ${chrom} ==="
      cat temp.${chrom}.zzz \
        | awk '{ print ">"$1,$2; print $3; }' \
        | fasta_set_line_length 100 \
        > targets/${mockName}.${chrom}.fa
      done

ls targets/${mockName}.chr*.fa \
  | sed "s/.*${mockName}\.//" \
  | sed "s/\.fa//" \
  | while read chrom ; do
      cat targets/${mockName}.${chrom}.fa \
        | fasta_N_intervals --complement \
        | awk '{ n++; t+=($4-$3)-(readLength-1) }
           END { print chrom,n,t }' \
           chrom=${chrom} readLength=${readLength}
      done \
  | awk 'BEGIN { print "#chrom intervals startPositions" } { print $0 }' \
  | line_up_columns \
  > targets/${mockName}.possible_reads

# sample reads from each chromosome

ls temp.*.zzz \
  | sed "s/temp\.//" \
  | sed "s/\.zzz//" \
  | while read chrom ; do
      echo "=== ${chrom} ==="
      genomeLenth=`cat targets/${mockName}.possible_reads \
                     | select_lines_by_name --name=${chrom} \
                     | awk '{ print $3 }'`
      #
      numReads=$((genomeLenth*sequencingDepth/readLength))
      echo "${chrom} ${numReads} reads" | commatize_column 2 --sep=space
      #      
      rootName=${mockName}.${chrom}.reads
      cat targets/${mockName}.${chrom}.fa \
        | awk '{ print $1 }' \
        | simulate_reads ${numReads}x${readLength} \
            --seed=${rootName} \
            --sam=reads/${rootName}.truth.sam.gz \
            --strand=forward --prohibitn --fastq --quality=I \
            --noise=${subsPerMillion}/1000000 \
            --indel=${indelsPerMillion}/1000000,0 \
            --name=R[9]_{chrom}_{zstart} \
            --progress=10K \
        | gzip \
        > reads/${rootName}.fastq.gz
    done

# convert the truth sams to sorted bams

ls temp.*.zzz \
  | sed "s/temp\.//" \
  | sed "s/\.zzz//" \
  | while read chrom ; do
      echo "=== ${chrom} ==="
      gzip -dc temp/${mockName}.${chrom}.reads.truth.sam.gz \
        | samtools sort \
            -T temp.samtools/${mockName}.${chrom}.reads.truth \
            -o alignments/${mockName}.${chrom}.reads.truth.bam
      samtools index \
        alignments/${mockName}.${chrom}.reads.truth.bam \
        alignments/${mockName}.${chrom}.reads.truth.bai
    done

