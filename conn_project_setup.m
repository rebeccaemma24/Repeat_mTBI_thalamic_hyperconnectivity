% Script to start up project definition and get ready for seed-based
% analysis. This imports new non-overlapping thalamus masks, and also runs
% ROI timeseries extraction on unsmoothed data. 


% add necessary paths
addpath(genpath('/spm12'));
addpath(genpath('/conn'));
%% 1. File setups before conn

% Go to working directory (i.e. the root path)
cd /rootpath
NSUBJECTS=43;
cwd = pwd; % Where you want this analysis saved

% Get all the possible files in the bids structure using the full path
% TBI
data = conn_bidsdir('/data/derivatives/fmriprep/');
allfiles = data.data.file;

% filter to just get the data with specific end patterns
idx = contains(allfiles,'_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii');
STRUCTURAL_FILE=allfiles(idx==1);
idx = contains(allfiles,'_space-MNI152NLin2009cAsym_label-CSF_probseg.nii');
CSFMASK_FILE=allfiles(idx==1);
idx = contains(allfiles,'_space-MNI152NLin2009cAsym_label-WM_probseg.nii');
WHITEMASK_FILE=allfiles(idx==1);
idx = contains(allfiles,'_space-MNI152NLin2009cAsym_label-GM_probseg.nii');
GREYMASK_FILE=allfiles(idx==1);
idx = contains(allfiles,'_space-MNI152NLin2009cAsym_desc-preproc_bold_smooth_clean_164.nii');
FUNCTIONAL_FILE=allfiles(idx==1);
idx = contains(allfiles,'_space-MNI152NLin2009cAsym_desc-preproc_bold_clean_164.nii');
FUNCTIONAL_UNSMOOTH=allfiles(idx==1);
clear idx

% Do a double check here that all the files are being found correctly
if rem(length(FUNCTIONAL_FILE),NSUBJECTS),error('mismatch number of functional files %n', length(FUNCTIONAL_FILE));end
if rem(length(FUNCTIONAL_UNSMOOTH),NSUBJECTS),error('mismatch number of functional unsmoothed files %n', length(FUNCTIONAL_UNSMOOTH));end
if rem(length(STRUCTURAL_FILE),NSUBJECTS),error('mismatch number of anatomical files %n', length(FUNCTIONAL_FILE));end
if rem(length(CSFMASK_FILE),NSUBJECTS),error('mismatch number of CSF files %n', length(FUNCTIONAL_FILE));end
if rem(length(WHITEMASK_FILE),NSUBJECTS),error('mismatch number of White Matter files %n', length(FUNCTIONAL_FILE));end
if rem(length(GREYMASK_FILE),NSUBJECTS),error('mismatch number of Grey Matter files %n', length(FUNCTIONAL_FILE));end

% Reshape data for CONN purposes
nsessions=length(FUNCTIONAL_FILE)/NSUBJECTS;
FUNCTIONAL_FILE=reshape(FUNCTIONAL_FILE,[nsessions,NSUBJECTS]);
FUNCTIONAL_UNSMOOTH=reshape(FUNCTIONAL_UNSMOOTH,[nsessions,NSUBJECTS]);
STRUCTURAL_FILE=reshape(STRUCTURAL_FILE,[nsessions, NSUBJECTS]);
GREYMASK_FILE=reshape(GREYMASK_FILE,[nsessions,NSUBJECTS]);
WHITEMASK_FILE=reshape(WHITEMASK_FILE,[nsessions,NSUBJECTS]);
CSFMASK_FILE=reshape(CSFMASK_FILE,[nsessions,NSUBJECTS]);

disp([num2str(size(FUNCTIONAL_FILE,1)),' sessions']); %Double-check
disp([num2str(size(FUNCTIONAL_FILE,2)),' subjects']);
TR=[];

% Groups setup
% vector of size [nsubjects,1] (with values from 1 to ngroup) defining subject groups
tbi=repmat(1,NSUBJECTS,1)';
GROUPS=tbi;

%% 2. CONN setup 
% New conn_*.mat experiment name
clear batch
batch.filename=fullfile(cwd,'conn_repeat_concussion.mat');
save('conn_repeat_concussion.mat')

% Setup interaction with cluster
batch.parallel.N=4; %NB this must match what is in the sbatch file under -n

% Setup
batch.Setup.isnew=1;
batch.Setup.nsubjects=NSUBJECTS;
batch.Setup.RT=TR;
%Point to func and struct files
batch.Setup.functionals=FUNCTIONAL_FILE;
batch.Setup.structurals=STRUCTURAL_FILE;
% Add secondary dataset of unsmoothed files
batch.Setup.secondarydatasets{1}.functionals_type = 4; % other, explicit
for nsub=1:NSUBJECTS
    batch.Setup.secondarydatasets{1}.functionals_explicit{nsub}{1}=FUNCTIONAL_UNSMOOTH{nsub};
end
% Setup Mask files for GM, WM, CSF
batch.Setup.masks.Grey = GREYMASK_FILE;
batch.Setup.masks.White = WHITEMASK_FILE;
batch.Setup.masks.CSF = CSFMASK_FILE;

% Setup conditions
batch.Setup.conditions.names={'rest'};

% Setup Subject Groups
batch.Setup.subjects.groups=GROUPS;
batch.Setup.subjects.group_names={'Group1'};

% Setup ROIs to add
% Thalamus_Najdenovska (L, R, and 7 subdivisions per hemisphere, max probability non-overlapping)
batch.Setup.rois.names={'Left_Thalamus', 'Right_Thalamus', 'LH_Pulvinar', 'LH_Anterior', 'LH_Medio_Dorsal', 'LH-Ventral_Latero_Dorsal', 'LH_Central', 'LH-Ventral_Anterior', 'LH-Ventral_Latero_Ventral', 'RH_Pulvinar', 'RH_Anterior', 'RH_Medio_Dorsal', 'RH-Ventral_Latero_Dorsal', 'RH_Central', 'RH-Ventral_Anterior', 'RH-Ventral_Latero_Ventral'};
batch.Setup.rois.files{1}={'/filepath/Thalamus_Lmask.nii'};
batch.Setup.rois.files{2}={'/filepath/Thalamus_Rmask.nii'};
batch.Setup.rois.files{3}={'/filepath/MaxProbMask01.nii'};
batch.Setup.rois.files{4}={'/filepath/MaxProbMask02.nii'};
batch.Setup.rois.files{5}={'/filepath/MaxProbMask03.nii'};
batch.Setup.rois.files{6}={'/filepath/MaxProbMask04.nii'};
batch.Setup.rois.files{7}={'/filepath/MaxProbMask05.nii'};
batch.Setup.rois.files{8}={'/filepath/MaxProbMask06.nii'};
batch.Setup.rois.files{9}={'/filepath/MaxProbMask07.nii'};
batch.Setup.rois.files{10}={'/filepath/MaxProbMask08.nii'};
batch.Setup.rois.files{11}={'/filepath/MaxProbMask09.nii'};
batch.Setup.rois.files{12}={'/filepath/MaxProbMask10.nii'};
batch.Setup.rois.files{13}={'/filepath/MaxProbMask11.nii'};
batch.Setup.rois.files{14}={'/filepath/MaxProbMask12.nii'};
batch.Setup.rois.files{15}={'/filepath/MaxProbMask13.nii'};
batch.Setup.rois.files{16}={'/filepath/MaxProbMask14.nii'};

batch.Setup.rois.dataset=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]; % Use 1st secondary dataset listed i.e. unsmoothed

% Setup analysis types to run
batch.Setup.analyses=[1,2]; % ROI-to-ROI and seed-to-voxel
% Setup units to allow for previous prepro and denoising
batch.Setup.analysisunits=2; % 2= raw units as using previously denoised data
batch.Setup.voxelresolution=3; % same as functionals

% Setup additional output file types 
batch.Setup.outputfiles(1)=0; %confound beta maps
batch.Setup.outputfiles(2)=1; %confound-corrected timeseries'
batch.Setup.outputfiles(3)=1; %seed-to-voxel r maps
batch.Setup.outputfiles(4)=1; %seed-to-voxel p maps
batch.Setup.outputfiles(5)=1; %seed-to-voxel FDR-p maps

batch.Setup.overwrite='Yes';
batch.Setup.done=1;

conn_batch(batch);

%% 3. Denoising
clear batch;
batch.filename=fullfile(cwd,'conn_repeat_concussion.mat');
batch.Denoising.confounds={}; % no removal of confounding effects
batch.Denoising.confounds.names={''}; % to indicate no confounds at all
batch.Denoising.filter=[0 inf]; % no filtering
batch.Denoising.detrending=0; % no detrending
batch.Denoising.despiking=0; % no despiking

batch.Denoising.overwrite='Yes';
batch.Denoising.done=1;

conn_batch(batch);

%% 4. First Level Analysis 
clear batch;
batch.filename=fullfile(cwd,'conn_repeat_concussion.mat');
batch.Analysis.measure=1; % analysis measure; 1 = correlation (bivariate)
batch.Analysis.sources={'Left_Thalamus', 'Right_Thalamus', 'LH_Pulvinar', 'LH_Anterior', 'LH_Medio_Dorsal', 'LH-Ventral_Latero_Dorsal', 'LH_Central', 'LH-Ventral_Anterior', 'LH-Ventral_Latero_Ventral', 'RH_Pulvinar', 'RH_Anterior', 'RH_Medio_Dorsal', 'RH-Ventral_Latero_Dorsal', 'RH_Central', 'RH-Ventral_Anterior', 'RH-Ventral_Latero_Ventral'}; % char array of source names (seeds)
batch.Analysis.type=3; % seed-to-voxel and ROI-to-ROI
batch.Analysis.weight=1; % no hrf weighting 

batch.Analysis.overwrite='Yes';
batch.Analysis.done=1;

conn_batch(batch);