#!/bin/bash

set -e        # stop the script if a command fails

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

function do_featurecounts_merge () {
	date
	########################################################################################################################
        sids="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/00-phenotype/sacc_subjects_pheno.txt"
	feature="gene"
#	brainregion="amygdala"
	brainregion="sacc"
	out_dir="/data/NHLBI_BCB/Fayaz/34-lncRNAKB-BP-JHU/featureCounts"
#	out_dir="featurecounts_exon"

	for i in `cat $sids | cut -f1 | head -1`; do
           SAMPLE=$i
	   echo -e "$feature\t$i" > tmp1
	   cat $out_dir/$SAMPLE/${SAMPLE}.${feature}featureCounts.txt | sed '1,2d' | sort -k1,1 | cut -f 1,7 >> tmp1
	done

	for i in `cat $sids | cut -f1 | sed '1,1d' `; do
           SAMPLE=$i
	   echo -e "$i" > tmp2
	   cat $out_dir/$SAMPLE/${SAMPLE}.${feature}featureCounts.txt | sed '1,2d' | sort -k1,1 | cut -f 7 >> tmp2
	   paste tmp1 tmp2 > tmp3
	   mv -f tmp3 tmp1
	done

	mv -f tmp1  $brainregion.$feature.featurecount.txt
	rm -f tmp2
	rm -f tmp3
	########################################################################################################################
	date
}

do_featurecounts_merge

