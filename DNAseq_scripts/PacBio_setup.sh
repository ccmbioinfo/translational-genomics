#!/bin/bash
# Usage: sh PacBio_setup.sh <analyses> <project>
# This file sets up crg2-pacbio analysis directories and downloads HPO terms and pedigrees from Phenotips. 
# It also replaces TCAG sequence IDs in VCFs with crg2-pacbio compatible IDs (ie family_sample), and checks that all samples in the analysis TSV are represented in the Phenotips pedigree.

set -euo pipefail


analyses=$1 # Path to sample sheet. Should be tab separated and three columns, Family_ID, sequence_id, Decoder_ID. Family_ID and sequence_id here refer to the IDs after which the files are named.
# For example family ID is 1745 and sample ID is 1741_SK0125, and the DECODER ID is DSK_007.01. 
project=$2 # e.g genesteps or DECODER

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

if [ -z $project ]; then
	echo "Please provide a project ID, e.g. DECODER"
	exit 1
fi

CREDS="PT_credentials.csv"
METADATA_DIR="/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}"
ANALYSIS_DIR="/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/analyses/${project}"
CRG2_PACBIO=~/crg2-pacbio
TODAY=`date +%Y-%m-%d`

if [ -z $analyses ]; then
        echo 'Must specify path to requested analysis TSV file, exiting'
        exit
fi

module load python/3.7.1
module load bcftools

# get HPO terms and pedigrees from Phenotips for all families in sample sheet
python3 get_HPO_pedigree_genome_clinic.py -sample_sheet $analyses -credentials $CREDS

for project_family in `awk '{print $3}' $analyses | cut -d '.' -f1 | tr -d '_' | uniq`
do
        echo $project_family
        FAMILY_DIR=${ANALYSIS_DIR}/${project_family}/PacBio_${TODAY}
        if [ -d ${FAMILY_DIR} ]; then
                echo "Analysis directory already exists for $project_family. Re-creating samples.tsv."
                # now re-initialize sample file
                echo -e "sample\tbam\tcase_or_control" > "${FAMILY_DIR}/samples.tsv"
        fi
done

while IFS=$'\t' read -r family sequence_id project_id sample_type status notes uploaded_to_DNAStack lab LIMS
do
        project_id=`echo $project_id | tr -d '_' | tr '.' '_' | tr -d '\r'` # e.g. convert DSK_018.03 to DSK018_03 so crg2-pacbio doesn't mess up wildcards
        project_family=`echo $project_id | cut -d '_' -f1`
	sequence_id=`echo $sequence_id | tr -d '\r'`
        echo $project_id  $project_family $sequence_id


    if [ "$family" != "Family_ID" ]; then
        # create analysis directory for family
        FAMILY_DIR=${ANALYSIS_DIR}/${project_family}/PacBio_${TODAY}
        if [ ! -d ${FAMILY_DIR} ]; then
                echo "Setting up analysis directory for $family"
                mkdir -p ${FAMILY_DIR}
                # copy pipeline files to analysis directory 
                cp ${CRG2_PACBIO}/config.yaml ${CRG2_PACBIO}/crg2-pacbio.sh ${CRG2_PACBIO}/slurm_profile/slurm-config.yaml ${FAMILY_DIR}
                sed -i "s/NA12878/$project_family/" ${FAMILY_DIR}/config.yaml

                # add HPO terms to config.yaml
                 echo "Finding HPO terms"
                 HPO=`ls -t /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/HPO/${project}/${project_family}* | head -n 1` # get most recent HPO file
                 echo $HPO
                 if [ -z $HPO ]; then
                         echo "No HPO terms found"
                 fi
                 sed -i "s+hpo: \"\"+hpo: \"${HPO}\"+"  ${FAMILY_DIR}/config.yaml

                # add pedigree to config.yaml
                 if [ "$project" = "DECODER" ]; then
                         echo "Finding pedigree"
                         ped=`ls -t /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/pedigrees/${project}/${project_family}* | head -n 1` # get most recent pedigree
                         echo $ped
                         if [ -z $ped ]; then
                                 echo "No pedigree found"
                         fi
                         sed -i "s+ped: \"\"+ped: \"${ped}\"+"  ${FAMILY_DIR}/config.yaml
                 fi


		# not C4R, so remove C4R sample columns from reports
		sed -i "s/c4r: True/c4r: False/"  ${FAMILY_DIR}/config.yaml

                # create samples.tsv 
                echo -e "sample\tBAM\tcase_or_control" > ${FAMILY_DIR}/samples.tsv

                # create units.tsv 
                echo "Populating units.tsv with inputs"
                echo -e "family\tplatform\tsmall_variant_vcf\tpbsv_vcf\tcnv_dir" > ${FAMILY_DIR}/units.tsv
                deepvariant=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${family}-cohort.joint.GRCh38.small_variants.phased.vcf.gz
                if [ ! -f $deepvariant ]; then
                        deepvariant=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${family}.joint.GRCh38.small_variants.phased.vcf.gz
                        if [ ! -f $deepvariant ]; then
                                tcag_family=`echo $family | sed 's/DSK_/DSK_DNA_/'`
                                deepvariant=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${tcag_family}-cohort.joint.GRCh38.small_variants.phased.vcf.gz
                                if [ ! -f $deepvariant ]; then
                                        deepvariant=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${sequence_id}.GRCh38.small_variants.phased.vcf.gz # singleton sample
                                fi
                        fi
                fi
                sv=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${family}-cohort.joint.GRCh38.structural_variants.phased.vcf.gz
                if [ ! -f $sv ]; then
                        sv=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${family}.joint.GRCh38.structural_variants.phased.vcf.gz
                        if [ ! -f $sv ]; then
                                tcag_family=`echo $family | sed 's/DSK_/DSK_DNA_/'`
                                sv=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${tcag_family}-cohort.joint.GRCh38.structural_variants.phased.vcf.gz
                                if [ ! -f $sv ]; then
                                        sv=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${sequence_id}.GRCh38.structural_variants.phased.vcf.gz  # singleton sample
                                fi
                        fi
                fi
                echo -e "$project_family\tPACBIO\t${deepvariant}\t${sv}\tcnv/vcfs" >> ${FAMILY_DIR}/units.tsv

                # create CNV directory
                mkdir ${FAMILY_DIR}/cnv

                # create trgt directory
                mkdir ${FAMILY_DIR}/trgt 
        fi

        # add sample and bams
        echo "Adding BAM file to samples.tsv for $sequence_id"
        BAM=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${sequence_id}.GRCh38.haplotagged.bam

        # add sample and BAM to samples.tsv
        project_sample=`echo $project_id | cut -d '_' -f2 | tr -d '\r' | tr -d ' '`
        echo -e "${project_sample}\t$BAM"  >> ${FAMILY_DIR}/samples.tsv

        # copy CNV VCFs
        echo "Copying CNV for $sequence_id"
        CNV=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${sequence_id}.GRCh38.hificnv.vcf.gz
        CNV_filename=`basename $CNV`
        if [ ! -d ${FAMILY_DIR}/cnv/vcfs ]; then
                mkdir -p ${FAMILY_DIR}/cnv/vcfs
        fi
        cp ${CNV}* ${FAMILY_DIR}/cnv/vcfs


        # trgt
        trgt=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${sequence_id}.GRCh38.trgt.sorted.phased.vcf.gz
        echo $trgt
        echo "Copying trgt for $sequence_id"
        if [ ! -d ${FAMILY_DIR}/trgt ]; then
                mkdir -p ${FAMILY_DIR}/trgt
        fi
        cp ${trgt}* ${FAMILY_DIR}/trgt
        trgt=`basename $trgt`
        trgt=${FAMILY_DIR}/trgt/$trgt

        # replace sample names in VCFs
        echo "Replacing sample IDs in VCFs"
        deepvariant=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${family}-cohort.joint.GRCh38.small_variants.phased.vcf.gz
        if [ ! -f $deepvariant ]; then
                deepvariant=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${family}.joint.GRCh38.small_variants.phased.vcf.gz
                if [ ! -f $deepvariant ]; then
                        tcag_family=`echo $family | sed 's/DSK_/DSK_DNA_/'`
			deepvariant=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${tcag_family}-cohort.joint.GRCh38.small_variants.phased.vcf.gz
			if [ ! -f $deepvariant ]; then
                                deepvariant=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${sequence_id}.GRCh38.small_variants.phased.vcf.gz # singleton sample
			fi
                fi
        fi
        sv=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${family}-cohort.joint.GRCh38.structural_variants.phased.vcf.gz
        if [ ! -f $sv ]; then
                sv=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${family}.joint.GRCh38.structural_variants.phased.vcf.gz
                if [ ! -f $sv ]; then
                        tcag_family=`echo $family | sed 's/DSK_/DSK_DNA_/'`
                        sv=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${tcag_family}-cohort.joint.GRCh38.structural_variants.phased.vcf.gz
			if [ ! -f $sv ]; then
                                sv=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${sequence_id}.GRCh38.structural_variants.phased.vcf.gz  # singleton sample
			fi
                fi
        fi
        # ID map file
        echo $sequence_id $project_id > sample_rename.txt
        # deepvariant
        bcftools reheader -s sample_rename.txt  $deepvariant |  bcftools view  -Oz -o $deepvariant.vcf.gz
        mv ${deepvariant}.vcf.gz $deepvariant
        # pbsv
        bcftools reheader -s sample_rename.txt  $sv |  bcftools view  -Oz -o $sv.vcf.gz
        mv $sv.vcf.gz $sv
        # trgt
        trgt_sequence_id=`bcftools query -l $trgt`
        echo $trgt_sequence_id $project_id > sample_rename.txt
        bcftools reheader -s sample_rename.txt  $trgt |  bcftools view  -Oz -o $trgt.vcf.gz
        mv $trgt.vcf.gz $trgt
        # hificnv
        for vcf in `ls ${FAMILY_DIR}/cnv/vcfs/*.vcf.gz`; do
            bcftools reheader -s sample_rename.txt  $vcf |  bcftools view  -Oz -o $vcf.vcf.gz
            mv $vcf.vcf.gz $vcf
            tabix $vcf # need to re-index otherwise bcftools merge complains
        done


        fi
done<$analyses

for project_family in `awk '{print $3}' $analyses | cut -d '.' -f1 |  uniq`
do
         echo $project_family
         if [ "$project_family" != "Decoder_ID" ]; then
              # check if all samples in analysis TSV are represented in the pedigree 
              ped_family=`echo $project_family | tr -d '_'`
              ped=`ls -t /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/pedigrees/${project}/${ped_family}* | head -n 1`
              # get all samples in the pedigree
              samples=`awk '{print $2}' $ped`
              for sample in `grep  $project_family $analyses | grep LRWGS |  awk '{print $3}' |  tr -d '_' | tr '.' '_'`; do
                  if ! echo $samples | grep -q $sample; then
                      echo "Sample $sample not found in pedigree for $project_family, exiting"
                      exit 1
                  fi
              done
         fi
done
