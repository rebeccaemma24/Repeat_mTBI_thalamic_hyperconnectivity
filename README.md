# Repeat_mTBI_thalamic_hyperconnectivity
Corresponding code to manuscript 'Repeat traumatic brain injury exacerbates acute thalamic hyperconnectivity in humans', Woodrow et al., 2024


- conn_project_setup : to obtain subject-level beta files for each thalamic nucleus (post-preprocessing)
- ComBat_voxel + get_voxels_long : to run voxel-level multicentre harmonisation
- regression_batch_job : to run voxelwise level group comparisons
- get_avBETA : to obtain global average thalamocortical functional connectivity estimates
- GlobalThalamocorticalAnalysis_RepeatmTBI : analysis of global thalamocortical functional connectivity between groups 
