#!/bin/bash
# Query the iRODS collection for files from PacBio humanwgs pipeline produced by TCAG

LIMSID=$1
PROJECT=$2
OUT=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${PROJECT}/iRods_paths/${LIMSID}_PacBio_files.txt
echo "Filepaths will be saved to $OUT"

# retrieve deep variant joint-genotyped VCF
{ iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.small_variants.phased.vcf.gz%'" | grep 'fam\|cohort' > $OUT; } || { iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.deepvariant.phased.vcf.gz%'" | grep 'fam\|cohort' > $OUT; }

# retrieve joint-genotyped SV VCF
{ iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.structural_variants.phased.vcf.gz%'" | grep 'fam\|cohort' >> $OUT; } || { iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.pbsv.phased.vcf.gz%'" | grep 'fam\|cohort' >> $OUT; }

# retrieve BAMs
{ iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.haplotagged.bam%'" >> $OUT; } || { iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%GRCh38.aligned.haplotagged.bam%'" >> $OUT; }

# retrieve CNV VCFs
{ iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%hificnv.vcf.gz%'" >> $OUT; } || { iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like 'hificnv%.vcf.gz%'" >> $OUT; }

# retrieve TRGT VCFs and spanning BAMs
{ iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%trgt.sorted.phased.vcf.gz%'" >> $OUT; } || { iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%trgt.sorted.vcf.gz%'" >> $OUT; }
iquest "%s/%s" "SELECT COLL_NAME, DATA_NAME where COLL_NAME like '%$LIMSID%' and DATA_NAME like '%.trgt.spanning.sorted.bam%'" >> $OUT

grep -v CAT_NO_ROWS_FOUND $OUT > $OUT.tmp
mv $OUT.tmp $OUT

