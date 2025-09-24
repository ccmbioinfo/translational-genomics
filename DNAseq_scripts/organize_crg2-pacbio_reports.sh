#!/bin/bash
# run from /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/analyses/
module load python/3.7.1

analysis=$1 # e.g. DSK002/PacBio
family=`dirname $analysis`
HPO=`ls -t /hpf/largeprojects/tgnode/sandbox/mcouse_analysis/HPO/*/${family}*txt | head -n 1`

mkdir ${analysis}/reports

# coding report
cp ${analysis}/small_variants/coding/*/*wes*  ${analysis}/reports
# panel report
panel=${analysis}/small_variants/panel/*/*wgs* 
panel_date=`basename $panel | cut -d'.' -f3`
cp $panel ${analysis}/reports/${family}.wgs.panel.${panel_date}.csv
# panel flank report
panelflank=${analysis}/small_variants/panel-flank/*/*wgs* 
panel_date=`basename $panelflank | cut -d'.' -f3`
cp $panelflank ${analysis}/reports/${family}.wgs.panel-flank.${panel_date}.csv

# add HPO terms to small variant reports
for f in ${analysis}/reports/*wes*csv ${analysis}/reports/*wgs*csv
do
	echo $f
	python3 ~/crg2-pacbio/utils/add_hpo_terms_to_wes.py $HPO $f
    rm $f
done

# cnv report
cp ${analysis}/cnv/${family}.cnv.202*.csv ${analysis}/reports
# sv report
cp ${analysis}/sv/*pbsv.202*csv ${analysis}/reports
# tr report
cp ${analysis}/pathogenic_repeats/${family}.known.path.str.loci.202*csv ${analysis}/reports
# tr outlier report
cp ${analysis}/repeat_outliers/${family}.repeat.outliers.annotated.202*csv  ${analysis}/reports 
# TRGT de novo report
cp ${analysis}/TRGT_denovo/${family}*202*csv ${analysis}/reports

