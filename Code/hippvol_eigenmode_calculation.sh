#!/bin/env bash

start_time=$(date +%s)

export work_dir="/YourPath/hippocampus_cortex_gradient_youth" #Add your absolute path of your working directory
export hippunfold_dir="${work_dir}/Data/hippunfold_HCD"
export hippovol_dir="${work_dir}/Data/hippovol_native"

export hemispheres='L R'
export MNI2mm="${work_dir}/Surf_temp/vol_templates/MNI152_T1_2mm.nii.gz"

export num_modes=31
export normalization_type='none'
export normalization_factor=1

cat ./Subjlist-all-HCD | while read line ; do
Subj_id="sub-${line}"

for hemi in ${hemispheres}; do
    
# generate the hippocampal mask with resolution of 2mm 
if [ ! -f "${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hippo-2mm.nii.gz" ]; then
    fslmaths "${hippunfold_dir}/${Subj_id}/anat/${Subj_id}_hemi-${hemi}_space-T1w_desc-subfields_atlas-bigbrain_dseg.nii.gz" -bin \
             "${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hippo.nii.gz"
    flirt -in "${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hippo.nii.gz" -ref ${MNI2mm} \
          -out "${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hippo-2mm.nii.gz" -applyisoxfm 2
    fslmaths "${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hippo-2mm.nii.gz" -thr 0.25 -bin \
             "${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hippo-2mm.nii.gz"
fi

# calculte the vol eigenmode of hippo
nifti_input_filename="${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hippo-2mm.nii.gz"
nifti_output_filename="${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hippo-2mm_emode_${num_modes}.nii.gz"
output_eval_filename="${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hippo-2mm_eval_${num_modes}.txt"
output_emode_filename="${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hippo-2mm_emode_${num_modes}.txt"

if [ ! -f ${output_emode_filename} ]; then
    echo Processing ${Subj_id}-${hemi}
     python ${work_dir}/Dependencies/Python/volume_eigenmodes.py ${nifti_input_filename} ${nifti_output_filename} \
                                 ${output_eval_filename} ${output_emode_filename} \
                                 -N ${num_modes} -norm ${normalization_type} -normfactor ${normalization_factor}
fi

# mapping vol eigenmode to the midthickness surface of hipp & DG
hipp_emode_func="${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-hipp-2mm_emode.func.gii"
DG_emode_func="${hippovol_dir}/${Subj_id}/${Subj_id}_hemi-${hemi}_space-T1w_desc-dentate-2mm_emode.func.gii"

if [ ! -f ${hipp_emode_func} ]; then
    wb_command -volume-to-surface-mapping ${nifti_output_filename} ${hippunfold_dir}/sub-${1}/surf/sub-${1}_hemi-${hemi}_space-T1w_den-2mm_label-hipp_midthickness.surf.gii ${hipp_emode_func} -enclosing
fi


if [ ! -f ${DG_emode_func} ]; then
    wb_command -volume-to-surface-mapping ${nifti_output_filename} ${hippunfold_dir}/sub-${1}/surf/sub-${1}_hemi-${hemi}_space-T1w_den-2mm_label-dentate_midthickness.surf.gii ${DG_emode_func} -enclosing
fi

done
done

end_time=$(date +%s)
echo "It took $(($end_time - $start_time)) seconds for the calculation of eigenmodes of hippocampal volume"