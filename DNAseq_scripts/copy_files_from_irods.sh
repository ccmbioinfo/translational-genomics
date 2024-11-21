#!/bin/bash
# take a list of files in the iRods archive and copy them to Isilon

irods_files=$1
project=$2
dir=/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/files_from_irods/${project}

if [ ! -d $dir ]; then
	echo "Directory does not exist, please double check project name"
	exit
fi

while read -r line
do
    echo $line
    iget $line $dir
done<$irods_files
