#!/bin/bash
#SBATCH --job-name=irsync-to-irods
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out
#SBATCH --mail-user=your.email@email.ca
#SBATCH --mail-type=END,FAIL

# This script archives directories on Isilon/PowerScale to iRODS.
# Note that it does not automatically delete the local files after syncing to iRODS. This can be done manually after the script finishes without error, but first run final_check.sh. 
# Usage: sh 01_irsync_to_iRODS.sh <local_path> <archival_path>
# Example usage: sbatch irsync_to_iRODS.sh CHU_SJ /resarchivezone/tgnode/external_datasets/CHUSJ

set -euo pipefail

local_path=$1 # local path on Isilon/PowerScale
archival_path=$2 # tgnode iRODS base path: /resarchivezone/tgnode/

# load iRODS module
module load irods_client/4.3.1

# set iRODS environment variables
export IRODS_ENVIRONMENT_FILE=~/.irods/irods_environment.json

if [ -z "$local_path" ] || [ -z "$archival_path" ]; then
    echo "Usage: $0 <local_path> <archival_path>"
    exit 1
fi

if [ ! -d "$local_path" ]; then
    echo "Local path $local_path does not exist"
    exit 1
fi

# sync the files to iRODS with -K flag for checksum verification
echo "syncing files from $local_path to $archival_path..."
irsync -rvK $local_path i:$archival_path
echo "irsync complete, syncing again to be safe..."

# just to be safe, sync the files to iRODS again (not sure if this is still an issue, https://github.com/irods/irods/issues/2822)
irsync -rvK $local_path i:$archival_path
echo "irsync complete!"
