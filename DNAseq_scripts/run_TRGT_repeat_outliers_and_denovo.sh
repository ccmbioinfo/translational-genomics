#!/bin/bash

family=$1
project=$2
ANALYSIS_DIR="/hpf/largeprojects/tgnode/sandbox/mcouse_analysis/analyses/${project}/${family}"

cp ~/crg2-pacbio/config.yaml ${ANALYSIS_DIR}/PacBio/
sed -i "s/NA12878/${family}/" ${ANALYSIS_DIR}/PacBio/config.yaml


HPO=`ls /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/HPO/${project}/${family}*`
echo $HPO
if [ -z $HPO ]; then
	echo "No HPO terms found"
fi
sed -i "s+hpo: \"\"+hpo: \"${HPO}\"+"  ${ANALYSIS_DIR}/PacBio/config.yaml

# add pedigree to config.yaml
echo "Finding pedigree"
ped=`ls /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/pedigrees/${project}/${family}*`
echo $ped
if [ -z $ped ]; then
	echo "No pedigree found"
fi
sed -i "s+ped: \"\"+ped: \"${ped}\"+"  ${ANALYSIS_DIR}/PacBio/config.yaml

# not C4R, so remove C4R sample columns from reports
              
sed -i "s/c4r: True/c4r: False/"  ${ANALYSIS_DIR}/PacBio/config.yaml

# add TRGT-denovo to job submission script
cp ~/crg2-pacbio/crg2-pacbio.sh ${ANALYSIS_DIR}/PacBio/
sed -i "s+{SLURM}+{SLURM} -p TRGT_denovo/${family}.TRGT.denovo.annotated.csv+g" ${ANALYSIS_DIR}/PacBio/crg2-pacbio.sh
