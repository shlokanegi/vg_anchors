#!/bin/bash

GENE=$1
THREADS=$2
FILE=$3
FILE_NAME="${FILE##*/}"
FILE_PREFIX="${FILE_NAME%%.*}"

echo "running minimap2 alignment of $FILE_PREFIX on $GENE"

if [ ! -d ../m2_alignments/$GENE/sam ]; then mkdir -p ../m2_alignments/$GENE/sam; fi
if [ ! -d ../m2_alignments/$GENE/bam ]; then mkdir -p ../m2_alignments/$GENE/bam; fi

if [ ! -f ../genes/$GENE.mmi ]; then minimap2 -d -k21 -w11 ../genes/$GENE.mmi ../genes/$GENE.mmi; fi
if [ ! -f "../m2_alignments/$GENE/sam/$FILE_PREFIX.sam" ]; then minimap2 -axsr --sam-hit-only ../genes/$GENE.mmi $FILE -t $THREADS > ../m2_alignments/$GENE/sam/$FILE_PREFIX.sam; fi

samtools view -b ../m2_alignments/$GENE/sam/$FILE_PREFIX.sam -o ../m2_alignments/$GENE/bam/$FILE_PREFIX.bam
samtools sort ../m2_alignments/$GENE/bam/$FILE_PREFIX.bam -o ../m2_alignments/$GENE/bam/$FILE_PREFIX.sorted.bam
samtools index ../m2_alignments/$GENE/bam/$FILE_PREFIX.sorted.bam


#!/bin/bash

THREADS=$1
GENE_FASTA=$2
echo "running minimap2 index"

GENE_MMI="${GENE_FASTA%.*}"
if [ ! -f "$GENE_MMI.mmi" ]; then echo "running minimap2 index on $GENE_FASTA" && minimap2 -d $GENE_MMI.mmi $GENE_FASTA -t $THREADS; fi
