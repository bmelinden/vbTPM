% Sample runinput file for vbTPM, ML 2013-09-05


%%%%%% part I: analysis parameters
% data preprocessing
fSample    =30;     % Hz, sampling frequency
downSample =3 ;     % downsampling factor

driftcorrection=true;     % do drift-correction with Butterworth filter
fCut=0.05;                % filter cut-off frequency for [Hz]
include_calibration=true; % include calibration traces in the analysis
include_looping=true;     % include looping traces in the analysis

% analysis parameters %
%%%%%%%%%%%%%%%%%%%%%%%
restarts_cal=10; % number of independent analyses of each calibration trajectory
Ninit_cal=10;    % number of initial states 

restarts_trj=40; % number of independent analyses of each active trajectory
Ninit_trj=20;    % number of initial states

plotAnything=false; % If true, shows various debug plots. Turn off for batch analysis.
one_at_a_time=true; % Analyze just one trajectory per call to the analysis function.
                    % This prevents memory problems in matlab, but requires
                    % repeated calls (e.g. a batch script) to the analysis
                    % function.
                    
% file names and paths
source_path     =   'lacdata/';     % path to xy data (calibration and looping)
target_path     =   'HMMresults2/'; % path to put results in 
PIDprefix       =   'RI2'; % prefix of results files

% debug flags
haveMaxLength=false; % If true, use only the first maxLength data points in 
%maxLength=1000;     % each trajectory. Good for fast testruns.

% data files and field names
calibration_xyfield='x';    % name of xy-position field in calibration and 
looping_xyfield ='x';    	% looping data files, as in 
                            % load(filname,nn_xyfield). Empty strings
                            % (nn_xyfield='') means that no field is read,
                            % as in load(filname). 
        
% Data files: the idea here is that each bead get a 'dry run' (no activity)
% for calibration, and a couple of active runs (perhaps interrupted by
% refocusing etc), and the analysis can keep track of each bead separately. 
calibration_filename={...
'bead1_cal_raw.mat', ...
'bead2_cal_raw.mat', ...
'bead3_cal_raw.mat', ...
'bead4_cal_raw.mat', ...
'bead5_cal_raw.mat'};

% production runs: one more level of nested cells, to allow for multiple
% production trajectories for each calibration run.
looping_filename={...
{'bead1_trj_raw.mat'}, ...
{'bead2_trj_raw.mat'}, ...
{'bead3_trj_raw.mat'}, ...
{'bead4_trj_raw.mat'}, ...
{'bead5_trj_raw.mat'}};

% a more convenient way might be to create these cell vectors as part of
% the analysis scripts, and simply load them here, e.g., 
% filenames = load('path_to_data/filenames.mat');
% calibration_filename = filenames.calibration_filename;
% looping_filename = filenames.looping_filename;

% construction of prior parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For an explanation of the parameters and their meaning, see help texts of
% VB7_priorParent. The suggested values are default values, and you might
% need to change them depending on your system or choice of units etc.

% genuine states
fB=1;
B0=5e-5;
K0=0.6;
Kstd=0.5;

fPi=5;       % total prior strength of initial state vector
tD=1;tA=5;
KBscaling=2; % 1) total emission strength constant 
             % 2) strength of single state constant (recommended).                          

% emission model, spurious states (spurious states are not used in the
% standard analysis, only after manual state classification).
fBc=1;
Bc0=5e-5;
Kc0=0.6;
KcStd=0.5;
tDc0=10;     % average dwell time of spurious state 0 (the non-spurious state). 

% parameters for the construction of an initial guess             %
% will need tuning to your system to make the analysis efficient. %
% See the help text of VB7_initialGuess_KBregion                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau0=5;  % typical dwell times for genuine,
tauS=1;  % and spurious, states [s]
wKBT=10000; % weight (pseudo-counts) of initial guess
tSigma=2; % Gaussian window width [s], for RMSKBgaussfilter

% range of K,B values in which to look for initial guesses. Experiment with
% RMSKBgaussfilter to choose appropriate intervals that include all
% expected parameter values. I recommend a wider interval for KBSrange than
% for KBGrange, since spurious states are more scattered in the parameter
% plane. Note that both intervals are used during analysis.
KBGrange=[0.3 0.9 0.2e-4 5e-4]; % for the genuine states
KBSrange=[0 1 0 1e-3];          % more spurious states. 
KBbins=500; % number of bins in the KBGrange histogram used for the initial guess
