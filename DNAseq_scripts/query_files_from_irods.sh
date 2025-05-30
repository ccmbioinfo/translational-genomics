#!/bin/bash
# Query the iRODS collection for files from PacBio humanwgs pipeline produced by TCAG

LIMSID=$1
PROJECT=$2
OUT=../../../files_from_irods/DECODER/iRods_paths/${LIMSID}_PacBio_files.txt

# retrieve deep variant joint-genotyped VCF
echo -e "iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like "%$LIMSID%" and DATA_NAME like '%cohort.joint.GRCh38.small_variants.phased.vcf.gz'" 
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.small_variants.phased.vcf.gz'" > $OUT
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.small_variants.phased.vcf.gz.tbi'" >> $OUT

# retrieve joint-genotyped SV VCF
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.structural_variants.phased.vcf.gz'" >> $OUT
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.structural_variants.phased.vcf.gz.tbi'" >> $OUT

# retrieve BAMs
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.haplotagged.bam'" >> $OUT
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.haplotagged.bam.bai'" >> $OUT

# retrieve CNV TSVs
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%_hg38_PACBIO_cnv.tagged.tsv'" >> $OUT

# retrieve TRGT VCFs and spanning BAMs
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%trgt.sorted.phased.vcf.gz'" >> $OUT
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%trgt.sorted.phased.vcf.gz.tbi'" >> $OUT   
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%.trgt.spanning.sorted.bam%'" >> $OUT

