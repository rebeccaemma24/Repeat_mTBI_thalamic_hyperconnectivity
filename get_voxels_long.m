% Script to get each BETA file in a column of voxels instead of 4d matrix
% One file for each of the beta value files, with each subject as a column

%---------------------------------------

outdir='/outdir/';
cwd = outdir;
rootdir='/rootdir/';
subjects = cellstr(num2str([1:43].','%03d'));
n_files = 16; %16
all_subs_long = zeros(1082035,length(subjects),n_files);

% setup loop to compute fnc matrices
for i=1:length(subjects)
    subject = subjects{i};

    for j=1:n_files
        n=num2str(j, '%02d');
    
        % get full file path for func (must already be unzipped)
        fn_nii = dir([rootdir 'BETA_Subject' subject '_Condition001_Source0' n '.nii']);
        fn_nii = fullfile({fn_nii.folder}, {fn_nii.name});
        nii_file=cell2mat(fn_nii);

        nii=niftiread(nii_file);
        % Change 0s to NaNs
        nii_signal=nii;
        nii_signal(nii_signal == 0) = NaN;
        % reshape to long
        nii_r = reshape(nii_signal, [], 1);
        all_subs_long(:,i,j) = nii_r;
    end    
end

% Now have the large mat file, split into its n_files to output
% Convert to cell array
C = num2cell(all_subs_long, [1 2]);
for k = 1:n_files
    n=num2str(k, '%02d');
    filename = ['BETA_long_Source0' n 'repeat.mat'];
    long_data = squeeze(C{k})';
    save(filename, 'long_data')
end
