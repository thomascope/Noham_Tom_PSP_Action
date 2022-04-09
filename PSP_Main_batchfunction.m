% This is a MEG/EEG preprocessing pipeline developed by Thomas Cope, based on a
% pipeline developed by Ed Sohoglu
%
% It assumes that you have already maxfiltered your data and created a
% trialDef file for each run containing the trial definitions (place these
% appropriately in each folder). I include an example python script to do
% this (you will need to change the paths etc).
%
% Before running, make sure that you have specified the parameters properly
% in this file, in the subjects_and_parameters script, and set the search terms in PSP_Preprocessing_mainfunction
% You will also need to move the bundled es_montage_all, tec_montage_all
% and MEGArtifactTemplateTopographies files to your pathstem folder.
% 
% Prerequisites:
% 
% A version of EEGLAB new enough to have the FileIO toolbox (the stock CBU version does not have this; I use 13_3_2b http://sccn.ucsd.edu/eeglab/downloadtoolbox.html )
% NB: You need to set the path to this in the PSP_Preprocessing_mainfunction case 'ICA_artifacts'
% The ADJUST toolbox for EEGLAB: (http://www.unicog.org/pm/pmwiki.php/MEG/RemovingArtifactsWithADJUST )
% 
%
% E-mail any questions to tec31@cam.ac.uk



%% Set up global variables
addpath(pwd)
pathstem = '/imaging/rowe/Noham_Tom_PSP_Action/preprocessed/';
maxfilteredpathstem = '/imaging/rowe/users/nw03/PSP_BP/data/scans/allscans/';
conditions = {'Reaction_time' 'Voluntary_button' 'Reaction_time_cue'};

contrast_labels = {'Sum RT';'Sum Voluntary';'Sum RT cue';'RT-voluntary';'voluntary-RT'};
contrast_weights = [1 0 0;0 1 0;0 0 1; 1 -1 0; -1 1 0];    

%% Specify participant parameters
PSP_demographics
group = [ones(1,length(patients))'; 2*ones(1,length(controls))']; % group 1 = patients, group 2 = controls
subjects = [patients; controls];
subjects = subjects';
group = group';

%% Specify preprocessing parameters

p.mod = {'MEGMAG' 'MEGPLANAR' 'EEG'}; % imaging modality (used by 'convert','convert+epoch','image','smooth','mask','firstlevel' steps) NB: EEG MUST ALWAYS BE LISTED LAST!!
p.ref_chans = {'EOG061','EOG062'}; %Reference channels in your montage in order horizEOG vertEOG ECG (NB: IF YOURS ARE DIFFERENT TO {'EOG061','EOG062','ECG063'} YOU WILL NEED TO MODIFY THE MONTAGE)

% cell array of experimental conditions (used by 'definetrials','convert+epoch' and 'sort' steps)
p.conditions = conditions;

%p.montage_fname = 'es_montage_MEG.mat'; % channel montage (used by 'convert','convert+epoch','artifact_ft','rereference' steps)
%p.montage_fname = 'es_montage_MEGPLANAR.mat'; % channel montage (used by 'convert','convert+epoch','artifact_ft','rereference' steps)
%p.montage_fname = 'es_montage_EEG.mat'; % channel montage (used by 'convert','convert+epoch','artifact_ft','rereference' steps)
% You need to set up some montages with and without the reference
% channels, and put them in the pathstem folder.
p.montage_fname = 'es_montage_all.mat'; % channel montage (used by 'convert','convert+epoch','artifact_ft','rereference' steps)
p.montage_rerefname = 'tec_montage_all.mat'; % channel montage (used by 'convert','convert+epoch','artifact_ft','rereference' steps)

p.fs = 1000; % original sample rate
p.fs_new = 250; % sample rate after downsampling in SPM (currently assumes that maxfilter HASN't downsampled data)

% for trial definitions
p.preEpoch = -2000; % pre stimulus time (ms)
p.postEpoch = 2000; % post stimulus time (ms)
p.triggers = [20 30 22]; % trigger values (correspond to p.conditions specified above)
p.minduration = 500; % if using definetrials_jp, minimum duration of a trial (ms)
p.maxduration = 100000; % if using definetrials_jp, maximum duration of a trial (ms)
%p.stimuli_list_fname = 'stimuli_list.txt';

% for robust averaging
p.robust = 1;

% for baseline correction
p.preBase = -2100; % pre baseline time (ms)
p.postBase = -2000; % post baseline time (ms)

% for combining planar gradiometer data
p.correctPlanar = 0; % whether to baseline correct planar gradiometer data after RMSing (using baseline period specified in preBase and postBase)

% for filtering 
p.filter = 'low'; % type of filter (lowpass or highpass)- never bandpass!
p.freq = 40; % filter cutoff (Hz)
%p.filter = 'high'; % type of filter (lowpass or highpass)- never bandpass!
%p.freq = 0.5; % filter cutoff (Hz)
%p.filter = 'stop';
%p.freq = [48 52];

% for computing contrasts of grand averaged MEEG data
p.contrast_labels = contrast_labels;
p.contrast_weights = contrast_weights;

% for image smoothing
p.xSmooth = 10; % smooth for x dimension (mm)
p.ySmooth = 10; % smooth for y dimension (mm)
p.zSmooth = 25; % smooth for z (time) dimension (ms)

% for making image mask
p.preImageMask = -2500; % pre image time (ms)
p.postImageMask = 1500; % post image time (ms)

% time windows over which to average for 'first-level' contrasts (ms)
p.windows = [-1000 0];

% set groups to input
p.group = group;
%% Event-related preprocessing steps

% note: should high-pass filter or baseline-correct before lowpass fitering
% to avoid ringing artefacts

% open up a parallel computing pool of appropriate size
% You should pilot one subject and see how much memory is required. This
% currently asks for 8Gb per run

workersrequested = length(subjects);

currentdr = pwd;
cd('/imaging/rowe/Noham_Tom_PSP_Action')
Poolinfo = cbupool(workersrequested,'--mem-per-cpu=8G --time=167:00:00');
parpool(Poolinfo,Poolinfo.NumWorkers);
cd(currentdr)

%Preprocessing pipeline below. This should be fully modular and intuitive.
%You will need to appropriately change the search terms at the start of the
%mainfunction to look for your dataname. If you change the order, these
%might need editing

%The mainfunction is called with the arguments
%PSP_Preprocessing_mainfunction(nextstep, previousstep, p, pathstem, maxfilteredpathstem, subjects{cnt}, cnt[, dates, blocksin, blocksout, rawpathstem, badeeg])
% previousstep can be specified using a name in the switch-case section of
% the mainfunction, or by entering a text searchstring (examples of each
% are below).

% Any functions present in the mainfunction and not listed below have not
% been optimised for SPM12 so will need some debugging before they will
% work.

data_exist = zeros(size(subjects));
parfor cnt = 1:size(subjects,2)
    try
        PSP_Preprocessing_mainfunction('definetrials_noham','convert',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
        data_exist(cnt) = 1;
    catch
        disp(['No data found for ' subjects{cnt}])
    end
end

% Exclude the subject for whom no data found
subjects(data_exist==0) = [];

%Now denoise the data
ICA_complete = zeros(size(subjects));
parfor cnt = 1:size(subjects,2)
   try
       PSP_Preprocessing_mainfunction('ICA_artifacts','convert',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt,dates,blocksin,blocksout,rawpathstem, badeeg, badchannels, runtodo)
       ICA_complete(cnt) = 1
       disp(['ICA complete for ' subjects{cnt}]);
   catch
       ICA_complete(cnt) = 0;
       disp(['ICA failed for ' subjects{cnt}]);
   end
end

% The first time you run through you should put a breakpoint on pause to check
% that the ICA worked correctly by inspecting todocomplete. The main
% reason for failure is if there are too many missing EEG electrodes for
% ADJUST to work. If this happens, then you can set the code to
% automatically continue without touching the EEG using a try - catch
% statement, or get it to omit EEG for specific subjects (example in the mainfunction) or try to fix it.

% pause

% This next bit is not yet fully parallelised (doesn't seem necessary
% because it doesn't take so long), so re-open a smaller matlab pool with
% one worker per subject to reduce occupation of the cluster

memoryrequired = num2str(8*size(subjects,1));

try
    matlabpool 'close'
    workerpool = cbupool(size(subjects,2));
    workerpool.ResourceTemplate=['-l nodes=^N^,mem=' memoryrequired 'GB,walltime=48:00:00'];
    matlabpool(workerpool)
catch
    fprintf([ '\n\nUnable to open up a cluster worker pool - opening a local cluster instead' ]);
    matlabpool(12)
end

parfor cnt = 1:size(subjects,2)
    %PSP_Preprocessing_mainfunction('definetrials','convert',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
   PSP_Preprocessing_mainfunction('epoch','ICA_artifacts',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt,dates,blocksin,blocksout,rawpathstem, badeeg);
end
parfor cnt = 1:size(subjects,2)
    PSP_Preprocessing_mainfunction('downsample','epoch',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end
% parfor cnt = 1:size(subjects,2)
%     PSP_Preprocessing_mainfunction('rereference','downsample',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
% end
parfor cnt = 1:size(subjects,2)
    PSP_Preprocessing_mainfunction('baseline','downsample',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end
parfor cnt = 1:size(subjects,2)
    PSP_Preprocessing_mainfunction('filter','baseline',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end
p.filestring_length = 4; %Specify that filenames are the first 4 unique letters and anything after this denotes a repeat (e.g. tray vs tray_1)
p.blocksout = blocksout;
parfor cnt = 1:size(subjects,2)
    PSP_Preprocessing_mainfunction('merge_recoded','filter',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end
p.conditions = postmerge_conditions;
parfor cnt = 1:size(subjects,2)
    PSP_Preprocessing_mainfunction('sort','merge',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end
parfor cnt = 1:size(subjects,2)
    PSP_Preprocessing_mainfunction('average','merge',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end
parfor cnt = 1:size(subjects,2)
    PSP_Preprocessing_mainfunction('filter','average',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end
parfor cnt = 1:size(subjects,2) 
    PSP_Preprocessing_mainfunction('combineplanar','fmcfbdeMr*.mat',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end

PSP_Preprocessing_mainfunction('grand_average','pfmcfbdeMr*.mat',p,pathstem, maxfilteredpathstem, subjects);
% This saves the grand unweighted average file for each group in the folder of the
% first member of that group. For convenience, you might want to move them
% to separate folders.


hasallconditions = zeros(1,size(subjects,2));
parfor cnt = 1:size(subjects,2)    
    try %Some participants didn't do all conditions, so can't be weighted with pre-specified contrasts.
   PSP_Preprocessing_mainfunction('weight','pfmcfbdeMr*.mat',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
   hasallconditions(cnt) = 1;
    catch
    end
end

PSP_Preprocessing_mainfunction('grand_average','wpfmcfbdeMr*.mat',p,pathstem, maxfilteredpathstem, subjects(logical(hasallconditions)));

parfor cnt = 1:size(subjects,2)
    PSP_Preprocessing_mainfunction('combineplanar_spm','fmcfbdeMr*.mat',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end
p.mod = {'MEGMAG' 'MEGCOMB'};
parfor cnt = 1:size(subjects,2)
    PSP_Preprocessing_mainfunction('image','PfmcfbdeMr*.mat',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end
parfor cnt = 1:size(subjects,2)
    % The input for smoothing should be the same as the input used to make
    % the image files.
    PSP_Preprocessing_mainfunction('smooth','PfmcfbdeMr*.mat',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end
for cnt = 1
    % The input for smoothing should be the same as the input used to make
    % the image files. Only need to do this for a single subject
    PSP_Preprocessing_mainfunction('mask','PfmcfbdeMr*.mat',p,pathstem, maxfilteredpathstem, subjects{cnt},cnt);
end  
% This saves the grand weighted average file for each group in the folder of the
% first member of that group. For convenience, you might want to move them
% to separate folders.

% now move all of the smoothed nifti images into folders marked
% controls/patients, either manually or with copyniftitofolder.py (you will
% need to change the paths and search characteristics appropriately).

matlabpool 'close';