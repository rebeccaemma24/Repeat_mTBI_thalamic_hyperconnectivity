#!/bin/bash
#
# For each BETA file, mask with participant's own GM mask and get average beta value

BETA_file='BETA_files.txt' #All subjects, 1st source only
GM_file='GM_cortex_files.txt'
outdir='/outdir'

# Loop through subjects first
for n in {1..43}; do
    b_fn=`sed "${n}q;d" ${BETA_file}` 
    BETA_fn=`echo ${b_fn} | sed s/......$//`
    GM_fn=`sed "${n}q;d" ${GM_file}`
    #GM_fn='tpl-MNI152NLin2009cAsym_res-02_label-GM_probseg_mask04.nii'

    # Loop through ROIs 1-16
    for i in {1..16}; do
        printf -v roi_no "%02d" ${i} 
        echo `fslstats ${BETA_fn}${roi_no}.nii -n -k ${GM_fn} -M` >> ${outdir}/avBETA_GM_cortex_0${roi_no}.txt
    done
done