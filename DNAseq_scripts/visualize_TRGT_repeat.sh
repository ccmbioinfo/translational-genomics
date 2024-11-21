#!/bin/bash

VCF=$1
BAM=$2
REPEAT_ID=$3
OUTPUT_DIR=`dirname $VCF`
SAMPLE=`basename $VCF | cut -d '.' -f1`
TRGTv1=/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/tools/TRGTv1.0.0/trgt
REF=/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_1/data/reference/human_GRCh38_no_alt_analysis_set.fasta
REPEATS=/hpf/largeprojects/ccmbio/ccmmarvin_shared/pacbio_longread_pilot_phase_2/tools/TRGTv1.0.0/repeats/pathogenic_repeats.hg38.bed


$TRGTv1 plot --genome $REF \
       --repeats $REPEATS \
       --vcf $VCF \
       --spanning-reads $BAM \
       --repeat-id ${REPEAT_ID} \
       --show meth \
       --image ${OUTPUT_DIR}/${SAMPLE}_${REPEAT_ID}.svg
