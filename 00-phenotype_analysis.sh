#!/bin/bash

set -e        # stop the script if a command fails

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

workingdir="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/00-phenotype"

cat /data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/differential_splicing_leafcutter/amyg_groups_diffsplice-input_20190502.txt | cut -d" " -f1-2 | sed 's/_/ /g' | cut -d" " -f 1,6 | sed 's/ /\t/g' > $workingdir/amygdala_subjects_pheno.txt

cat $workingdir/amygdala_subjects_pheno.txt | sort -k2 > $workingdir/tmp
mv -f $workingdir/tmp $workingdir/amygdala_subjects_pheno.txt

cat /data/NHLBI_BCB/Fayaz/31-leafcutter/lieber/differential_splicing_leafcutter/sacc_cov_groups_file.txt | cut -d" " -f1-2 | sed 's/_/ /g' | cut -d" " -f 1,6 | sed 's/ /\t/g' > $workingdir/sacc_subjects_pheno.txt

cat $workingdir/sacc_subjects_pheno.txt | sort -k2 > $workingdir/tmp
mv -f $workingdir/tmp $workingdir/sacc_subjects_pheno.txt



