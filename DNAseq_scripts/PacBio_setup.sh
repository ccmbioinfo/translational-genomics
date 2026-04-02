#!/bin/bash

module load python/3.7.1
module load tabix
module load bcftools

set -euo pipefail

analyses=$1 # Path to sample sheet. Should be tab separated and three columns, Family_ID, sequence_id, Decoder_ID. Family_ID and sequence_id here refer to the IDs after which the files are named.
# For example family ID is 1745 and sample ID is 1741_SK0125, and the DECODER ID is DSK_007.01. 
project=$2 # e.g genesteps or DECODER

# Input validation
if [ -z $analyses ]; then
	echo "Please provide a sample metadata TSV"
	exit 1
fi

if [ -z $project ]; then
	echo "Please provide a project ID, e.g. DECODER"
	exit 1
fi

python3 PacBio_setup.py --analyses $analyses --project $project
