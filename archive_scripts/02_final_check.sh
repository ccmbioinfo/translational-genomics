#!/bin/bash
# This script double checks that archiving was successful by comparing the number of files in the local directory to the number of files in the archive.
# Thanks Lynette Lau for the inspo! 
# If the check is "OK", you can delete the local directory that you archived. If there is a "WARNING", some files were not transferred or mirrored appropriately.
# Usage: sh 02_final_check.sh <local_path> <archival_path>
# Example usage: sh final_check.sh /hpf/largeprojects/tgnode/data/CHU_SJ /resarchivezone/tgnode/external_datasets/CHUSJ

set -euo pipefail

local_path=$1
archival_path=$2

module load irods_client/4.3.1

if [ -z "$local_path" ] || [ -z "$archival_path" ]; then
    echo "Usage: $0 <local_path> <archival_path>"
    exit 1
fi

if [ ! -d "$local_path" ]; then
    echo "Local path $local_path does not exist"
    exit 1
fi

n_archive=$(ils -lr ${archival_path} | grep " & " | wc -l) # lists all the files in the archive directory, including mirrored files
n_local=$(find $local_path -type f | wc -l)
n_local_2=$((${n_local} * 2)) # the number of local files multiplied by two should equal the number of files in the archive directory
if [[ $n_archive -ne $n_local_2 ]]; then
        WARN="WARNING" # if these are not equal, some files were not transferred or mirrored appropriately
else
        WARN="OK"
fi

echo -e "iRODS_count_with_mirror\tIsilon_count_x2\tcheck"
echo -e "$n_archive\t$n_local_2\t$WARN"

