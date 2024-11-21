#!/bin/bash
#usage: sh PacBio_setup.sh <analyses> <iRods_query_date> <project>

set -euo pipefail

analyses=$1 # Path to sample sheet. Should be tab separated and three columns, Family_ID, Sample_ID, Decoder_ID. Family_ID and Sample_ID here refer to the IDs after which the files are named.
# For example family ID is 1745 and sample ID is 1741_SK0125, and the DECODER ID is DSK_007.01. If files are named with DECODER ID, just use the DECODER family and sample IDs.
date=$2 # Date that iRods files were queried/retrieved, e.g. 2024-08-27 
project=$3 # e.g genesteps or DECODER

# Input validation
if [ -z $analyses ]; then
	echo "Please provide a sample metadata TSV"
	exit 1
fi

# Check if input file exists
if [ ! -f "$analyses" ]; then
    echo "Error: Sample file $sample_file does not exist"
    exit 1
fi

if [ -z $date ]; then
	echo "Please the date that iRods was queried"
	exit 1
fi

if [ -z $project ]; then
	echo "Please provide a project ID, e.g. DECODER"
	exit 1
fi

CREDS="PT_credentials.csv"
METADATA_DIR="/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/sample_sheets"
ANALYSIS_DIR="/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/analyses/${project}"
CRG2_PACBIO=~/crg2-pacbio

if [ -z $analyses ]; then
        echo 'Must specify path to requested analysis TSV file, exiting'
        exit
fi

module load python/3.7.1


python3 get_HPO_pedigree_genome_clinic.py $analyses $CREDS

while IFS=$'\t' read -r family sample_id project_id_original
do
        project_id=`echo $project_id_original | tr -d '_' | tr '.' '_'` # e.g. convert DSK_018.03 to DSK018_03 so crg2 doesn't mess up wildcards
        project_family_original=`echo $project_id_original | cut -d '.' -f1`
        project_family=`echo $project_id_original | cut -d '.' -f1 | tr -d '_' `

    if [ "$family" != "Family_ID" ]; then
        # create analysis directory for family
        if [ ! -d ${ANALYSIS_DIR}/${project_family} ]; then
                echo "Setting up analysis directory for $family"
                mkdir -p ${ANALYSIS_DIR}/${project_family}/PacBio
                # copy pipeline files to analysis directory 
                cp ${CRG2_PACBIO}/config.yaml ${CRG2_PACBIO}/crg2-pacbio.sh ${CRG2_PACBIO}/slurm_profile/slurm-config.yaml ${ANALYSIS_DIR}/${project_family}/PacBio
                sed -i "s/NA12878/$project_family/" ${ANALYSIS_DIR}/${project_family}/PacBio/config.yaml

                # add HPO terms to config.yaml
                echo "Finding HPO terms"
                HPO=`ls /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/HPO/${project}/${project_family}*`
                echo $HPO
                if [ -z $HPO ]; then
                        echo "No HPO terms found"
                fi
                sed -i "s+hpo: \"\"+hpo: \"${HPO}\"+"  ${ANALYSIS_DIR}/${project_family}/PacBio/config.yaml

                # add pedigree to config.yaml
                echo "Finding pedigree"
                ped=`ls /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/pedigrees/${project}/${project_family}*`
                echo $ped
                if [ -z $ped ]; then
                        echo "No pedigree found"
                fi
                sed -i "s+ped: \"\"+ped: \"${ped}\"+"  ${ANALYSIS_DIR}/${project_family}/PacBio/config.yaml


		# not C4R, so remove C4R sample columns from reports
		c4r: True
		sed -i "s/c4r: True/c4r: False/"  ${ANALYSIS_DIR}/${project_family}/PacBio/config.yaml
                # add targets to job submission script
                sed -i "s+{SLURM}+{SLURM} -p sv/${project_family}.pbsv.csv small_variants/coding/${project_family} small_variants/panel/${project_family} small_variants/panel-flank/${project_family}   pathogenic_repeats/${project_family}.known.path.str.loci.csv repeat_outliers/${project_family}.repeat.outliers.annotated.csv TRGT_denovo/${project_family}.TRGT.denovo.annotated.csv+g" ${ANALYSIS_DIR}/${project_family}/PacBio/crg2-pacbio.sh

                # create samples.tsv 
                echo -e "sample\tBAM" > ${ANALYSIS_DIR}/${project_family}/PacBio/samples.tsv

                # create units.tsv 
                echo "Populating units.tsv with inputs"
                echo -e "family\tplatform\tsmall_variant_vcf\ttrgt_vcf_dir\tpbsv_vcf\ttrgt_pathogenic_vcf_dir" > ${ANALYSIS_DIR}/${project_family}/PacBio/units.tsv
                deepvariant=`grep $family ${METADATA_DIR}/deepvariant_${date}.txt`
                deepvariant=`basename $deepvariant`
                deepvariant=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${deepvariant}
                pbsv=`grep $family ${METADATA_DIR}/pbsv_${date}.txt`
                pbsv=`basename $pbsv`
                pbsv=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${pbsv}
                echo -e "$project_family\tPACBIO\t${deepvariant}\ttrgt\t${pbsv}\t" >> ${ANALYSIS_DIR}/${project_family}/PacBio/units.tsv

                # create CNV directory
                mkdir ${ANALYSIS_DIR}/${project_family}/PacBio/cnv

                # create trgt directory
                mkdir ${ANALYSIS_DIR}/${project_family}/PacBio/trgt
        fi

        # add sample and bams
        echo "Adding BAM file to samples.tsv for $sample_id"
        BAM=`grep $sample_id ${METADATA_DIR}/BAMs_${date}.txt`
        BAM=`basename $BAM`
        BAM=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${BAM}

        # add sample and BAM to samples.tsv
        project_sample=`echo $project_id_original | cut -d '.' -f2`
        echo -e "${project_sample}\t$BAM"  >> ${ANALYSIS_DIR}/${project_family}/PacBio/samples.tsv

        # copy TCAG annotated CNV files for merging
        echo "Copying CNV for $sample_id"
        CNV=`grep $sample_id /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/sample_sheets/CNV_${date}.txt`
        CNV_filename=`basename $CNV`
        CNV=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${CNV_filename}
        cp $CNV ${ANALYSIS_DIR}/${project_family}/PacBio/cnv
        # rename CNV file
        sample_id_no_period=`echo $sample_id | tr '.' '_'` # TCAG replaces period in DECODER IDs with underscore
        sample_id_no_period=`echo ${sample_id_no_period}_A1`
        sed -i "s/${sample_id_no_period}/${project_id}/g" ${ANALYSIS_DIR}/${project_family}/PacBio/cnv/${CNV_filename}

        # trgt
        trgt=`grep $sample_id_no_period ${METADATA_DIR}/trgt_${date}.txt`
        trgt=`basename $trgt`
        trgt=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${trgt}
        echo $trgt
        mkdir ${ANALYSIS_DIR}/${project_family}/PacBio/trgt
        cp $trgt ${ANALYSIS_DIR}/${project_family}/PacBio/trgt
        trgt=`basename $trgt`
        trgt=${ANALYSIS_DIR}/${project_family}/PacBio/trgt/$trgt

        # replace sample names in VCFs
        echo "Replacing sample IDs in VCFs"
        deepvariant=`grep $family ${METADATA_DIR}/deepvariant_${date}.txt`
        deepvariant=`basename $deepvariant`
        deepvariant=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${deepvariant}
        pbsv=`grep $family ${METADATA_DIR}/pbsv_${date}.txt`
        pbsv=`basename $pbsv`
        pbsv=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${pbsv}
        # ID map file
        echo $sample_id_no_period $project_id > sample_rename.txt
        # deepvariant
        bcftools reheader -s sample_rename.txt  -o $deepvariant.vcf.gz $deepvariant
        mv ${deepvariant}.vcf.gz $deepvariant
        # pbsv
        bcftools reheader -s sample_rename.txt  -o $pbsv.vcf.gz $pbsv
        mv $pbsv.vcf.gz $pbsv
        # trgt
        trgt_sample_id=`bcftools query -l $trgt`
        echo $trgt_sample_id $project_id > sample_rename.txt
        bcftools reheader -s sample_rename.txt  -o $trgt.vcf.gz $trgt
        mv $trgt.vcf.gz $trgt

        if [ "$sample_id" != "$project_id_original" ]; then
                echo "Sample ID and DECODER ID do not match, renaming files"
                echo $sample_id $project_id_original
                # # rename sample in BAM
                # BAM_prefix=`echo $BAM | sed 's/.bam//g'`
                # samtools view -H $BAM  | sed "s/SM:[^\t]*/SM:${project_id}/g" | samtools reheader - $BAM > ${BAM_prefix}.rename.bam
                # samtools index ${BAM_prefix}.rename.bam
                # sed -i  "s+$BAM+${BAM_prefix}.rename.bam+" ${ANALYSIS_DIR}/${project_family}/PacBio/samples.tsv
                
        fi
        fi
done<$analyses

#submit one CNV merge job per family
for project_family in `awk '{print $3}' $analyses | cut -d '.' -f1 | tr -d '_' |  uniq`
do
         if [ "$project_family" != "ProjectID" ]; then
              cd ${ANALYSIS_DIR}/${project_family}/PacBio/cnv
              echo "Submitting CNV merge job"
              sbatch ${CRG2_PACBIO}/scripts/merge.cnv.reports.sh $project_family
              cd $ANALYSIS_DIR
         fi
done
