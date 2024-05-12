%-----------------------------------------------------------------------
% Job saved on 26-Apr-2023 13:52:47 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7487)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

% Edit input directories
%addpath(genpath('/spm12'));
root_dir = '/rootdir/';

% Setup indexes for subject inclusion
HC_subjects = cellstr(num2str([1:76].','%03d'));

matched_TBI = readtable('/Matched_mTBI.csv');
w = str2double(matched_TBI.Weight);
idx1 = find(w == 1);
mTBI_subjects = cellstr(num2str((idx1 + 76),'%03d')); % list of all mTBI subs numbers

num_rep = readtable('/num_repeats_covars.csv');
groups = str2double(num_rep.combined_concussion);
idx2 = find(groups == 2);
repeat_subjects = cellstr(num2str((idx2 + 184),'%03d')); % list of all subs with 2+ previous injuries n=21

n_files=16; %number of analyses to run 

%-----------------------------------------------------------------------
% Setup the regressors
regressors = readtable('/combat_covars.csv');
HC_age = table2array(regressors([1:76],5));
HC_sex = table2array(regressors([1:76],6));

mTBI_age = table2array(regressors([idx1 + 76],5));
mTBI_sex = table2array(regressors([idx1 + 76],6));

repeat_age = table2array(regressors([idx2 + 184],5));
repeat_sex = table2array(regressors([idx2 + 184],6));

age = [HC_age; mTBI_age; repeat_age];
sex = str2double([HC_sex; mTBI_sex; repeat_sex]);

% Regressor of interest: group
HCr = repmat(1,76,1);
mTBIr = repmat(2,21,1);
repeatr = repmat(3,21,1);
group = [HCr ; mTBIr ; repeatr];

%-----------------------------------------------------------------------
% Setup Loop to repeat this for all the components
for i = 1:n_files

    file_no = num2str(i,'%02d');
    file_mask = ['/OneSample/Analysis0' file_no '_mask.nii'];
    output_basename = ['Analysis_0' file_no];
    output_dir = ['/MatchedData/Regression2/Analysis_0' file_no '/'];
    mkdir(fullfile(output_dir))
    
    % Subjects Setup:
    % NB format: rootdir/Analysis00*_harmonized/BETA_Subject***_Source00*_harmonized.nii
    % Get the character length from spm_select, then select the required imgs

    file_folder = ['Analysis_0' file_no '_harmonized'];
    f = spm_select('FPList', fullfile(root_dir, file_folder), '.*_Subject001_Source0.._harmonized\.nii$');

    HC_files=blanks(length(f));
    for j=1:length(HC_subjects)
        subject = HC_subjects{j};
        fn_func = dir([root_dir 'Analysis_0' file_no '_harmonized/BETA_Subject' subject '_Source0' file_no '_harmonized.nii']);
        sub_file = cell2mat(fullfile({fn_func.folder}, {fn_func.name}));
        HC_files(j,:) = sub_file;
    end

    mTBI_files=blanks(length(f));
    for j=1:length(mTBI_subjects)
        subject = mTBI_subjects{j};
        fn_func = dir([root_dir 'Analysis_0' file_no '_harmonized/BETA_Subject' subject '_Source0' file_no '_harmonized.nii']);
        sub_file = cell2mat(fullfile({fn_func.folder}, {fn_func.name}));
        mTBI_files(j,:) = sub_file;
    end

    repeat_files=blanks(length(f));
    for j=1:length(repeat_subjects)
        subject = repeat_subjects{j};
        fn_func = dir([root_dir 'Analysis_0' file_no '_harmonized/BETA_Subject' subject '_Source0' file_no '_harmonized.nii']);
        sub_file = cell2mat(fullfile({fn_func.folder}, {fn_func.name}));
        repeat_files(j,:) = sub_file;
    end

    ALL_files = [HC_files; mTBI_files; repeat_files];

%-----------------------------------------------------------------------
% warm up the spm_jobman
spm('defaults','fmri'); 
spm_jobman('initcfg');
matlabbatch = [];

% Design specification
matlabbatch{1}.spm.stats.factorial_design.dir = { output_dir };
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = cellstr(ALL_files);
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.c = [group];
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.cname = 'Group';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov.iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [age];
matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'Age';
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).c = [sex];
matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'Sex';
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = { file_mask };
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
% Model Estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
% Contrast Manager
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Group+';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Group-';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 0 0 -1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

% Run
spm_jobman('run',matlabbatch);

% Clean up for next run
clear HC_files
clear mTBI_files
clear repeat_files
clear ALL_files
clear f
clear matlabbatch

end