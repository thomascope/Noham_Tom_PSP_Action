%% Set up global variables

%clear all

% make sure EEG modality of SPM software is selected
%spm('EEG');
%spm

% add required paths
addpath(pwd);
addpath('/group/language/data/ediz.sohoglu/matlab/utilities/');
% addpath('/opt/neuromag/meg_pd_1.2/');

% define paths
%define where the rawdata is kept (not actually important for this script,
%because maxfiltering done separately, e.g.
%rawpathstem = '/imaging/rowe/users/nw03/PSP_BP/data/scans/allscans';
maxfilteredpathstem = '/imaging/rowe/users/nw03/PSP_BP/data/scans/allscans/';

% define conditions
conditions = {'Reaction_time' 'Voluntary_button'};
% postmerge_conditions = {'Standard_onset fbdeMrun_play_raw_ssst' 'Standard_onset fbdeMrun_tray_raw_ssst' 'Standard_onset fbdeMrun_qway_raw_ssst' 'Deviant_D_onset fbdeMrun_play_raw_ssst'  'Deviant_D_onset fbdeMrun_tray_raw_ssst' 'Deviant_D_onset fbdeMrun_qway_raw_ssst' 'Deviant_T_onset fbdeMrun_play_raw_ssst'  'Deviant_T_onset fbdeMrun_tray_raw_ssst' 'Deviant_T_onset fbdeMrun_qway_raw_ssst' 'Standard_offset fbdeMrun_play_raw_ssst'  'Standard_offset fbdeMrun_tray_raw_ssst' 'Standard_offset fbdeMrun_qway_raw_ssst' 'Deviant_D_offset fbdeMrun_play_raw_ssst'  'Deviant_D_offset fbdeMrun_tray_raw_ssst' 'Deviant_D_offset fbdeMrun_qway_raw_ssst' 'Deviant_T_offset fbdeMrun_play_raw_ssst'  'Deviant_T_offset fbdeMrun_tray_raw_ssst' 'Deviant_T_offset fbdeMrun_qway_raw_ssst'};
% postmerge_conditions_tf = {'Standard_onset effbdMrun_play_raw_ssst' 'Standard_onset effbdMrun_play_raw_ssst' 'Standard_onset effbdMrun_qway_raw_ssst' 'Deviant_D_onset effbdMrun_play_raw_ssst'  'Deviant_D_onset effbdMrun_play_raw_ssst' 'Deviant_D_onset effbdMrun_qway_raw_ssst' 'Deviant_T_onset effbdMrun_play_raw_ssst'  'Deviant_T_onset effbdMrun_play_raw_ssst' 'Deviant_T_onset effbdMrun_qway_raw_ssst' 'Standard_offset effbdMrun_play_raw_ssst'  'Standard_offset effbdMrun_play_raw_ssst' 'Standard_offset effbdMrun_qway_raw_ssst' 'Deviant_D_offset effbdMrun_play_raw_ssst'  'Deviant_D_offset effbdMrun_play_raw_ssst' 'Deviant_D_offset effbdMrun_qway_raw_ssst' 'Deviant_T_offset effbdMrun_play_raw_ssst'  'Deviant_T_offset effbdMrun_play_raw_ssst' 'Deviant_T_offset effbdMrun_qway_raw_ssst'};

contrast_labels = {'Sum all standard onsets';'Sum all D deviant onsets';'Sum all T deviant onsets';'Sum all onsets';'Sum all offsets'; 'Onset_D-S'; 'Onset_Dd-S'; 'Onset_Dt-S'; 'Onset_Dd-Dt'; 'Offset_D-S'; 'Offset_Dd-S'; 'Offset_Dt-S'; 'Offset_Dd-Dt'; 'Play-Tray'; 'Play-Qway'; 'Tray-Qway'};
contrast_weights = [kron([1, 0, 0, 0, 0, 0],[1,1,1]);kron([0, 1, 0, 0, 0, 0],[1,1,1]);kron([0, 0, 1, 0, 0, 0],[1,1,1]);kron([1, 1, 1, 0, 0, 0],[1,1,1]);kron([ 0, 0, 0, 1, 1, 1],[1,1,1]);kron([ -1, 0.5, 0.5, 0, 0, 0],[1,1,1]);kron([ -1, 1, 0, 0, 0, 0],[1,1,1]);kron([ -1, 0, 1, 0, 0, 0],[1,1,1]);kron([ 0, 1, -1, 0, 0, 0],[1,1,1]);kron([ 0, 0, 0, -1, 0.5, 0.5],[1,1,1]);kron([ 0, 0, 0, -1, 1, 0],[1,1,1]);kron([ 0, 0, 0, -1, 0, 1],[1,1,1]);kron([ 0, 0, 0, 0, 1, -1],[1,1,1]);kron([1, 0, 0, 0, 0, 0],[1,-1,0]);kron([1, 0, 0, 0, 0, 0],[1,0,-1]);kron([1, 0, 0, 0, 0, 0],[0,1,-1])];    

% define subjects and blocks (group(cnt) = 1 for controls, group(cnt) = 2 for patients)
cnt = 0;
% 
cnt = cnt + 1;
subjects{cnt} = 'meg12_0458';
dates{cnt} = '080512';
blocksin{cnt} = {'0252_play', '0252_qway', '0252_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0259';
dates{cnt} = '080516';
blocksin{cnt} = {'0259_play', '0259_qway', '0259_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0260';
dates{cnt} = '080516';
blocksin{cnt} = {'0260_play', '0260_qway', '0260_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0264';
dates{cnt} = '080519';
blocksin{cnt} = {'0264_play', '0264_qway', '0264_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0265';
dates{cnt} = '080520';
blocksin{cnt} = {'0265_play', '0265_qway', '0265_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0266';
dates{cnt} = '080520';
blocksin{cnt} = {'0266_play', '0266_qway', '0266_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0292';
dates{cnt} = '080709';
blocksin{cnt} = {'0292_play', '0292_qway', '0292_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0293';
dates{cnt} = '080709';
blocksin{cnt} = {'0293_play', '0293_qway', '0293_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0295';
dates{cnt} = '080710';
blocksin{cnt} = {'0295_play', '0295_qway', '0295_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0307';
dates{cnt} = '080718';
blocksin{cnt} = {'0307_play', '0307_qway', '0307_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0311';
dates{cnt} = '080723';
blocksin{cnt} = {'0311_play', '0311_qway', '0311_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0330';
dates{cnt} = '080804';
blocksin{cnt} = {'0330_play', '0330_qway', '0330_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0338';
dates{cnt} = '080808';
blocksin{cnt} = {'0338_play', '0338_qway', '0338_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg08_0462';
dates{cnt} = '081121';
blocksin{cnt} = {'0462_play', '0462_qway', '0462_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 1;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg09_0183';
dates{cnt} = '090623';
blocksin{cnt} = {'0183_play', '0183_qway', '0183_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg09_0185';
dates{cnt} = '090625';
blocksin{cnt} = {'0185_play', '0185_qway', '0185_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg09_0186';
dates{cnt} = '090625';
blocksin{cnt} = {'0186_play', '0186_qway', '0186_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg09_0228';
dates{cnt} = '090811';
blocksin{cnt} = {'0228_play_raw_1', '0228_play_raw_2', '0228_play_raw_3', '0228_qway_raw_1', '0228_qway_raw_2', '0228_qway_raw_3', '0055_tray_rawreal_1', '0055_tray_raw_2', '0228_tray_raw_3'};
blocksout{cnt} = {'play_1', 'play_2', 'play_3', 'qway_1', 'qway_2', 'qway_3', 'tray_1', 'tray_2', 'tray_3'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ; badchannels{cnt, 4} = {} ;badchannels{cnt, 5} = {} ; badchannels{cnt, 6} = {}; badchannels{cnt, 7} = {}; badchannels{cnt, 8} = {}; badchannels{cnt, 9} = {};% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg10_0026';
dates{cnt} = '100304';
blocksin{cnt} = {'0526_play', '0526_qway', '0526_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg10_0314';
dates{cnt} = '101102';
blocksin{cnt} = {'0314_play_run1', '0314_play_run2', '0314_play_run3', '0314_tray_run1', '0314_tray_run2', '0314_tray_run3'};
blocksout{cnt} = {'play_1', 'play_2', 'play_3', 'tray_1', 'tray_2', 'tray_3'};
badchannels{cnt, 1} = {'0111'}; badchannels{cnt, 2} = {'0111'}; badchannels{cnt, 3} = {'0111'} ; badchannels{cnt, 4} = {'0111'}; badchannels{cnt, 5} = {'0111'}; badchannels{cnt, 6} = {'0111'} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 13; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg12_0208';
dates{cnt} = '120516';
blocksin{cnt} = {'0208_play', '0208_tray'};
blocksout{cnt} = {'play', 'tray'};
badchannels{cnt, 1} = {'1611', '2323'}; badchannels{cnt, 2} = {'1611', '2323'}; badchannels{cnt, 3} = {'1611', '2323'} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 32; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg12_0261';
dates{cnt} = '120614';
blocksin{cnt} = {'0261_play', '0261_qway', '0261_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {'1732'}; badchannels{cnt, 2} = {'1732'}; badchannels{cnt, 3} = {'1732'} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 32; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg12_0262';
dates{cnt} = '120614';
blocksin{cnt} = {'0262_play', '0262_qway', '0262_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 32; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg12_0359';
dates{cnt} = '120718';
blocksin{cnt} = {'0359_play', '0359_qway', '0359_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {}; badchannels{cnt, 2} = {}; badchannels{cnt, 3} = {} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 32; % delay time (ms) between trigger and stimulus

cnt = cnt + 1;
subjects{cnt} = 'meg12_0404';
dates{cnt} = '120815';
blocksin{cnt} = {'0404_play', '0404_qway', '0404_tray'};
blocksout{cnt} = {'play', 'qway', 'tray'};
badchannels{cnt, 1} = {'2621'}; badchannels{cnt, 2} = {'2621'}; badchannels{cnt, 3} = {'2621'} ;% define bad MEG (not EEG) channels here (if there are any)
badeeg{cnt} = {};
group(cnt) = 2;
p.delay{cnt} = 32; % delay time (ms) between trigger and stimulus