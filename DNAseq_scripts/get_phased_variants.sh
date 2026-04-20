#!/bin/bash 
# This script will query all variants in a specified phase block (haplotype block) that are in phase with a variant of interest
# To obtain the phase set ID and phased genotype, first query the variant coordinates, i.e. tabix $VCF chr:pos-end
# e.g. sh get_phased_variants.sh  /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/DECODER/DSK_DNA_130-cohort.joint.GRCh38.small_variants.phased.vcf.gz  61735909 DSK_DNA_130-03-A1 "1|0"

VCF=$1 # family joint-genotyped VCF
PS=$2 # PS is the phase set ID for the variant in question, e.g. a de novo variant in a proband. 
SAMPLE=$3
PHASED_GT=$4 # phased genotype for variant in question, e.g. "0|1"

SAMPLE_CRG2=`echo $SAMPLE | sed 's/_DNA_//' | sed 's/-A1//' | sed 's/-/_/g'`
OUTDIR="../../../results/phased_variants/"
echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGENOTYPE" > $OUTDIR/${SAMPLE}_variants.txt 
bcftools view $VCF -s $SAMPLE | bcftools filter -i "TYPE='snp' & FORMAT/PS = $PS" | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%SAMPLE=%GT]\n' | grep "$PHASED_GT" >> $OUTDIR/${SAMPLE_CRG2}_variants.txt
