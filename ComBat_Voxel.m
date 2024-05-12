%% Voxelwise Harmonization
%
% Packages: MATLAB2017b and above, ComBat, Tools for NIfTI and ANALYZE image
% Requirements: Data in long format (all voxels in one column per subject),
% batch regressors, example of output datatype for header/setup info

%% Load in the necessary data
addpath(genpath('/ComBatHarmonization-master'));
addpath(genpath('/NIfTI_20140122'));
%--------------------------------------------
% Only things to change:
rootdir = '/rootdir/';
n_files = 16;
regressors = readtable([rootdir 'combat_covars.csv']);
matlab_nii = load_untouch_nii('/data/BETA_Subject001_Condition001_Source001.nii'); % Existing nii file to setup the output options

%% Get the regressors into the necessary format
[Sites, ~, SiteCode] = unique(regressors(:,4));
SiteBatch = SiteCode.';
mod = regressors{:,[3,5,6]};

%% Load in Data
for i=1:n_files
    n=num2str(i, '%02d');
    
    input_fn = ['BETA_long_Source0' n '_ALL.mat'];
    load(input_fn) % named BETA_long data in working folder
    outdir = [rootdir 'Analysis_0' n '_harmonized/'];
    mkdir(fullfile(outdir))
% Deal with full rows of NaNs in data (to be added back in after for continuity)
    long_data = long_data.';
    NanRows = all(isnan(long_data),2); %logical array of where all values are NaN
    long_signal = long_data(~all(isnan(long_data),2),:); % remove rows of all NaN

    % If any NaNs remain, change these to 0 
    % NB Matlab's version of combat can't handle missing values
    long_signal(isnan(long_signal))=0;

% Run ComBat
    % data_harmonized = combat(dat, batch, mod, 1);
    long_signal_harmonized = combat(long_signal, SiteBatch, mod, 1);

% Add back in the NaN rows 

    pos     = find(NanRows == 1).'; % find the row positions of where Nans should be
    [r,c]   = size(long_signal_harmonized); % the number of rows and columns in old
    add     = numel(pos); % how much longer the new matrix will be
    long_harm = NaN(r + add,c); % Preallocate the size of the new matrix
    idx     = setdiff(1:r+add,pos); % all positions of new matrix except pos

    long_harm(idx,:) = long_signal_harmonized;

% Output this mat file as harmonized in current directory
    mat_file = ['BETA_long_Source0' n '_harmonized.mat'];
    save(mat_file, 'long_harm');

% Converting back to nifti: relies on existing file's setup

    % Loop through to write to individual files
    subjects = cellstr(num2str([1:227].','%03d'));

    for k = 1:length(subjects)
        % Change from long format back into the right dimensions
        long_harm_sub = reshape(long_harm(:,k), 97, 115, 97);
        % Change the img data
        matlab_nii.img = long_harm_sub;
        % Save file with output name
        subject = subjects{k};
        filename = [outdir 'BETA_Subject' subject '_Source0' n '_harmonized.nii'];
        save_untouch_nii(matlab_nii, filename);
    end

end