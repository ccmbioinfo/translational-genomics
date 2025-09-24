#!/bin/bash
# run from /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/analyses/
module load python/3.7.1

analysis=$1 # e.g. DSK002/WGS_reanalysis
ngs=$2

family=`dirname $analysis`
HPO=$(ls /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/HPO/*/${family}* | tail -n 1)


mkdir ${analysis}/report/all

if [[ "$ngs" == "WGS" ]]; then
    # MELT
    cp ${analysis}/report/MELT/* ${analysis}/report/all
    # coding report
    cp ${analysis}/report/coding/*/*wes*  ${analysis}/report/all
    # de novo report
    cp ${analysis}/report/denovo/*/*denovo*  ${analysis}/report/all
    # panel report
    panel=`ls ${analysis}/report/panel/*/*wgs* | grep -v clinical` 
    panel_date=`basename $panel | cut -d'.' -f3`
    cp $panel ${analysis}/report/all/${family}.wgs.panel.${panel_date}.csv
    # panel flank report
    panelflank=`ls ${analysis}/report/panel-flank/*/*wgs* | grep -v clinical` 
    panel_date=`basename $panelflank | cut -d'.' -f3`
    cp $panelflank ${analysis}/report/all/${family}.wgs.panel-flank.${panel_date}.csv
    # sv report
    cp ${analysis}/report/sv/* ${analysis}/report/all
    # MT report
    cp ${analysis}/report/mitochondrial/*202* ${analysis}/report/all
    # STR reports
    cp ${analysis}/report/str/*202* ${analysis}/report/all

    rm ${analysis}/report/all/*clinical*
else
    # coding report
    cp ${analysis}/report/coding/*/*wes*  ${analysis}/report/all
    # somatic/mosaic report
    cp ${analysis}/report/gatk_somatic/*/*wes*  ${analysis}/report/all
fi 

# add HPO terms to small variant reports
for f in ${analysis}/report/all/*w*csv
do
    if [ "$f" = *"MELT"* ]; then
        python3 ~/crg2-pacbio/utils/add_hpo_terms_to_wes.py $HPO $f
        rm $f
    fi
done

