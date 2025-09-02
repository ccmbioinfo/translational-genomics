#!/bin/bash
#usage: sh PacBio_setup.sh <analyses> <project>

set -euo pipefail


analyses=$1 # Path to sample sheet. Should be tab separated and three columns, Family_ID, sequence_id, Decoder_ID. Family_ID and sequence_id here refer to the IDs after which the files are named.
# For example family ID is 1745 and sample ID is 1741_SK0125, and the DECODER ID is DSK_007.01. If files are named with DECODER ID, just use the DECODER family and sample IDs.
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
ANALYSIS_DIR="/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/analyses/test/${project}"
CRG2_PACBIO=~/crg2-pacbio

if [ -z $analyses ]; then
        echo 'Must specify path to requested analysis TSV file, exiting'
        exit
fi

module load python/3.7.1
module load bcftools


#python3 get_HPO_pedigree_genome_clinic.py -sample_sheet $analyses -credentials $CREDS

for project_family in `awk '{print $3}' $analyses | cut -d '.' -f1 | tr -d '_' | uniq`
do
        echo $project_family
        if [ -d ${ANALYSIS_DIR}/${project_family}/PacBio ]; then
                echo "Analysis directory already exists for $project_family. Re-creating samples.tsv and units.tsv files."
                # now re-initialize sample file
                echo "sample" > "${ANALYSIS_DIR}/${project_family}/PacBio/samples.tsv"
        fi
done

while IFS=$'\t' read -r family sequence_id project_id sample_type status notes uploaded_to_DNAStack lab LIMS
do
        project_id=`echo $project_id | tr -d '_' | tr '.' '_' | tr -d '\r'` # e.g. convert DSK_018.03 to DSK018_03 so crg2 doesn't mess up wildcards
        project_family=`echo $project_id | cut -d '_' -f1`
	sequence_id=`echo $sequence_id | tr -d '\r'`
        echo $project_id  $project_family $sequence_id


    if [ "$family" != "Family_ID" ]; then
        # create analysis directory for family
        if [ ! -d ${ANALYSIS_DIR}/${project_family}/PacBio ]; then
                echo "Setting up analysis directory for $family"
                mkdir -p ${ANALYSIS_DIR}/${project_family}/PacBio
                # copy pipeline files to analysis directory 
                cp ${CRG2_PACBIO}/config.yaml ${CRG2_PACBIO}/crg2-pacbio.sh ${CRG2_PACBIO}/slurm_profile/slurm-config.yaml ${ANALYSIS_DIR}/${project_family}/PacBio
                sed -i "s/NA12878/$project_family/" ${ANALYSIS_DIR}/${project_family}/PacBio/config.yaml

                # add HPO terms to config.yaml
                 echo "Finding HPO terms"
                 HPO=`ls -t /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/HPO/${project}/${project_family}* | head -n 1` # get most recent HPO file
                 echo $HPO
                 if [ -z $HPO ]; then
                         echo "No HPO terms found"
                 fi
                 sed -i "s+hpo: \"\"+hpo: \"${HPO}\"+"  ${ANALYSIS_DIR}/${project_family}/PacBio/config.yaml

                # add pedigree to config.yaml
                 if [ "$project" = "DECODER" ]; then
                         echo "Finding pedigree"
                         ped=`ls /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/pedigrees/${project}/${project_family}* | tail -n 1` # last file is the most recent
                         echo $ped
                         if [ -z $ped ]; then
                                 echo "No pedigree found"
                         fi
                         sed -i "s+ped: \"\"+ped: \"${ped}\"+"  ${ANALYSIS_DIR}/${project_family}/PacBio/config.yaml
                 fi


		# not C4R, so remove C4R sample columns from reports
		sed -i "s/c4r: True/c4r: False/"  ${ANALYSIS_DIR}/${project_family}/PacBio/config.yaml
                # add targets to job submission script
                if [ "$project" = "genesteps" ]; then
                        sed -i "s+{SLURM}+{SLURM} -p sv/${project_family}.pbsv.csv small_variants/coding/${project_family} small_variants/panel/${project_family} small_variants/panel-flank/${project_family}   pathogenic_repeats/${project_family}.known.path.str.loci.csv repeat_outliers/${project_family}.repeat.outliers.annotated.csv+g" ${ANALYSIS_DIR}/${project_family}/PacBio/crg2-pacbio.sh
                fi
                # create samples.tsv 
                echo -e "sample\tBAM\tcase_or_control" > ${ANALYSIS_DIR}/${project_family}/PacBio/samples.tsv

                # create units.tsv 
                echo "Populating units.tsv with inputs"
                echo -e "family\tplatform\tsmall_variant_vcf\ttrgt_vcf_dir\tpbsv_vcf\ttrgt_pathogenic_vcf_dir" > ${ANALYSIS_DIR}/${project_family}/PacBio/units.tsv
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
                echo -e "$project_family\tPACBIO\t${deepvariant}\ttrgt\t${sv}\t" >> ${ANALYSIS_DIR}/${project_family}/PacBio/units.tsv

                # create CNV directory
                mkdir ${ANALYSIS_DIR}/${project_family}/PacBio/cnv

                # create trgt directory
                mkdir ${ANALYSIS_DIR}/${project_family}/PacBio/trgt 
        fi

        # add sample and bams
        echo "Adding BAM file to samples.tsv for $sequence_id"
        BAM=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${sequence_id}.GRCh38.haplotagged.bam

        # add sample and BAM to samples.tsv
        project_sample=`echo $project_id | cut -d '_' -f2 | tr -d '\r' | tr -d ' '`
        echo -e "${project_sample}\t$BAM"  >> ${ANALYSIS_DIR}/${project_family}/PacBio/samples.tsv

        # copy TCAG annotated CNV files for merging
        echo "Copying CNV for $sequence_id"
        CNV=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${sequence_id}.GRCh38.hificnv_annPipelineRevv1.7.1_20250120_hg38_PACBIO_cnv.tagged.tsv
        CNV_filename=`basename $CNV`
        if [ ! -d ${ANALYSIS_DIR}/${project_family}/PacBio/cnv ]; then
                mkdir -p ${ANALYSIS_DIR}/${project_family}/PacBio/cnv
        fi
        cp $CNV ${ANALYSIS_DIR}/${project_family}/PacBio/cnv
        # rename CNV file
        sequence_id_no_period=`echo $sequence_id | tr '.' '_'` # TCAG replaces period in DECODER IDs with underscore
        project_id=`echo $project_id | tr -d '\r' | tr -d ' '`
        #sequence_id_no_period=`echo ${sequence_id_no_period}_A1`
        sed -i "s/${sequence_id_no_period}/${project_id}/g" ${ANALYSIS_DIR}/${project_family}/PacBio/cnv/${CNV_filename}

        # trgt
        trgt=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}/${sequence_id}.GRCh38.trgt.sorted.phased.vcf.gz
        echo $trgt
        echo "Copying trgt for $sequence_id"
        if [ ! -d ${ANALYSIS_DIR}/${project_family}/PacBio/trgt ]; then
                mkdir -p ${ANALYSIS_DIR}/${project_family}/PacBio/trgt
        fi
        cp ${trgt}* ${ANALYSIS_DIR}/${project_family}/PacBio/trgt
        trgt=`basename $trgt`
        trgt=${ANALYSIS_DIR}/${project_family}/PacBio/trgt/$trgt

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
        echo $sequence_id_no_period $project_id > sample_rename.txt
        # deepvariant
        bcftools reheader -s sample_rename.txt  -o $deepvariant.vcf.gz $deepvariant
        mv ${deepvariant}.vcf.gz $deepvariant
        # pbsv
        bcftools reheader -s sample_rename.txt  -o $sv.vcf.gz $sv
        mv $sv.vcf.gz $sv
        # trgt
        trgt_sequence_id=`bcftools query -l $trgt`
        echo $trgt_sequence_id $project_id > sample_rename.txt
        bcftools reheader -s sample_rename.txt  -o $trgt.vcf.gz $trgt
        mv $trgt.vcf.gz $trgt

        fi
done<$analyses

for project_family in `awk '{print $3}' $analyses | cut -d '.' -f1 |  uniq`
do
         echo $project_family
         if [ "$project_family" != "Decoder_ID" ]; then
              # check if all samples in analysis TSV are represented in the pedigree 
              ped_family=`echo $project_family | tr -d '_'`
              ped=`ls /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/pedigrees/${project}/${ped_family}* | tail -n 1`
              # get all samples in the pedigree
              samples=`awk '{print $2}' $ped`
              # check if all samples in analysis TSV are represented in the pedigree
              for sample in `grep  $project_family $analyses | grep LRWGS |  awk '{print $3}' |  tr -d '_' | tr '.' '_'`; do
                  if ! echo $samples | grep -q $sample; then
                      echo "Sample $sample not found in pedigree for $project_family, exiting"
                      exit 1
                  fi
              done
              # exclude TRGT de novo rule if family is not a trio 
              project_family=$(echo "$project_family" | tr -d '_' | tr -d '\r')
              if [ $(wc -l ${ANALYSIS_DIR}/${project_family}/PacBio/samples.tsv | awk '{print $1}') -lt 4 ]; then
                echo "Family $project_family is not a trio, excluding TRGT de novo rule"
                sed -i "s+{SLURM}+{SLURM} -p sv/${project_family}.pbsv.csv small_variants/coding/${project_family} small_variants/panel/${project_family} small_variants/panel-flank/${project_family}   pathogenic_repeats/${project_family}.known.path.str.loci.csv repeat_outliers/${project_family}.repeat.outliers.annotated.csv+g" ${ANALYSIS_DIR}/${project_family}/PacBio/crg2-pacbio.sh
              fi
              # submit one CNV merge job per family
              cd ${ANALYSIS_DIR}/${project_family}/PacBio/cnv
              echo "Submitting CNV merge job"
              sbatch ${CRG2_PACBIO}/scripts/merge.cnv.reports.sh $project_family
              cd $ANALYSIS_DIR
         fi
done
