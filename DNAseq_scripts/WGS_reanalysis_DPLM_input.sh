#!/bin/bash
# Set up crg2 pipeline run for sequences from SickKids DPLM clinical genetics lab
# Usage: WGS_reanalysis_DPLM_input.sh <sample TSV file> <project>
# Example:  sh WGS_reanalysis_DPLM_input.sh ../sample_sheets/DECODER/analysis_requests/analyses_DSK054_DSK045.tsv DECODER

set -euo pipefail

# Input validation
sample_file=$1
project=$2

if [ -z $sample_file ]; then
	echo "Please provide a sample metadata TSV"
	exit 1
fi

# Check if input file exists
if [ ! -f "$sample_file" ]; then
    echo "Error: Sample file $sample_file does not exist"
    exit 1
fi

if [ -z $project ]; then
	echo "Please provide a project ID, e.g. DECODER"
	exit 1
fi

CREDS="PT_credentials.csv"
METADATA_DIR="/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/sample_sheets"
ANALYSIS_DIR="/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/analyses/${project}"
CRG2=~/crg2
FASTQ_DIR="/hpf/largeprojects/tgnode/data/fastq_20241120"


while IFS=$'\t' read -r family sample_id decoder_id sample_type status notes dnastack lab 
do
	# Skip header
	[ "$decoder_id" = "DecoderID" ] && continue
	
	# Format IDs
	crg2_decoder_id=$(echo "$decoder_id" | tr -d '_' | tr '.' '_') # e.g. convert DSK_018.03 to DSK018_03 so crg2 doesn't mess up wildcards
	crg2_sample_id=$(echo "$crg2_decoder_id" | cut -d '_' -f2)
	family=$(echo "$decoder_id" | cut -d '.' -f1 | tr -d '_')
	
	ANALYSIS_DIR_FAM="${ANALYSIS_DIR}/${family}/${sample_type}_reanalysis"
	echo "Processing family: $family"
	echo "Analysis directory: $ANALYSIS_DIR_FAM"

	# Create and setup analysis directory if it doesn't exist
	if [ ! -d "$ANALYSIS_DIR_FAM" ]; then
		echo "Setting up analysis directory for $family"
		mkdir -p "$ANALYSIS_DIR_FAM"
		
		# Copy config files
		cp "${CRG2}/config_hpf.yaml" "${CRG2}/dnaseq_slurm_hpf.sh" "${CRG2}/slurm_profile/slurm-config.yaml" "$ANALYSIS_DIR_FAM"
		
		# Update config file
		sed -i "s/NA12878/$family/" "${ANALYSIS_DIR_FAM}/config_hpf.yaml"
		[ "$sample_type" = "WGS" ] && sed -i "s/wes/wgs/" "${ANALYSIS_DIR_FAM}/config_hpf.yaml"
		
		# Set HPO terms and pedigree
		HPO="../HPO/${project}/${family}*"
		pedigree="../pedigrees/${project}/${family}*"
		if [ ! -f $HPO ]; then
			echo "Warning: No HPO file found for family ${family}"
		else
			sed -i "s+hpo: \"\"+hpo: \"${HPO}\"+" "${ANALYSIS_DIR_FAM}/config_hpf.yaml"
		fi
		if [ ! -f $pedigree ]; then
			echo "Warning: No pedigree file found for family ${family}" 
		else
			sed -i "s+ped: \"\"+ped: \"${pedigree}\"+" "${ANALYSIS_DIR_FAM}/config_hpf.yaml"
		fi
		
		# Initialize sample files
		echo "sample" > "${ANALYSIS_DIR_FAM}/samples.tsv"
		echo -e "sample\tplatform\tfq1\tfq2\tbam\tcram" > "${ANALYSIS_DIR_FAM}/units.tsv"
	fi
	
	# Add sample to samples.tsv
	echo "$crg2_sample_id" >> "${ANALYSIS_DIR_FAM}/samples.tsv"

	# Find and format fastq paths
	R1=$(find "${FASTQ_DIR}" -name "*${sample_id}*R1*fastq.gz" | paste -sd "," -)
	R2=$(find "${FASTQ_DIR}" -name "*${sample_id}*R2*fastq.gz" | paste -sd "," -)
	
	if [ -z "$R1" ] || [ -z "$R2" ]; then
		echo "Warning: No fastq files found for sample ${sample_id}"
	fi
	
	echo -e "$crg2_sample_id\tILLUMINA\t$R1\t$R2\t\t" >> "${ANALYSIS_DIR_FAM}/units.tsv"
	
done < "${sample_file}"
