% skeleton of a runinput file for script. 
% The plan is that this script should supply all specific information for
% analyzing a particular set of data. It is then combined with a generic
% analysis script that runs the analysis based on the information in this
% file. 
% The runinput file is constructed by the GUI or a scripting user for each
% data set, and then processed by the analysis program.
% Hopefully, I will be able to supply several analysis scripts (e.g., with
% ot without using the parallel computing toolbox) that all work with the
% generic runinput file.
% M.L. 2011-02-07

%%%%%% part I: parameters
% data preprocessing
fSample    =30;     % Hz
downSample =10;     % downsampling factor

driftcorrection=false;      % true = do drift-correction with Butterworth filter
fCut=0.05;                  % filter cut-off frequency for [Hz]
include_calibration=true; % include calibration traces in the analysis
include_looping=true;     % include looping traces in the analysis

% construction of prior parameters
% emission model, genuine states
fPi=[]; % empty means use default values. 
fB=[];
B0=5e-5;
K0 =0.6;Kstd =0.5;
tD=2;       % average dwell time in prior distribution [s]
tA=2;       % strength of transition rate prior        [s]

KBscaling=2; % 1) total emission strength constant 
             % 2) strength of single state constant (recommended)
% emission model, spurious states (not activated)
fBc=[];
Bc0=5e-5;
Kc0=0.6;KcStd=0.5;
tDc0=10;    % average dwell time of spurious state 0 (the non-spurious state). 

% construction initial guess 
tau0=5;  % typical dwell times for genuine 
tauS=1;  % and spurious states
wKBT=[];
KBbins=[];
tSigma=2;       % Gaussian window width [s]
KBGrange=[0.3 0.9 0.2e-4 1e-4];
KBSrange=[0 1 0 5e-4];

% analysis
restarts_cal=10;  % number of analysis rounds per calibration trace
Ninit_cal=10;     % initial number of states for calibration traces

restarts_trj=30;  % number of analysis rounds per looping trace
Ninit_trj=20;     % initial number of states for looping traces

plotAnything=false;
one_at_a_time=true; % analyze just one trajectory per call to the analysis function 
                    % prevents memory problems in matlab, but requires
                    % repeated calls (e.g. a batch script) to the analysis
                    % function
% file names and paths
source_path     =   '../VB7_PUCdata/'; % where to find the data files
target_path     =   '../VB7_PUCdata/'; % where to put the result files
PIDprefix       =   'PUC306_100pMSJLacI_30'; % what to name the result files

% debug flags
haveMaxLength=false;
%maxLength=1000;

%%%%%% part II: data file names
% the data needs to be in the field x, so that the command
% load(filename,'x') loads the data to the workspace.
calibration_filename={...
'PUC306_100pMSJLacI_bead30_cal_corr.mat',...
'PUC306_100pMSJLacI_bead31_cal_corr.mat',...
'PUC306_100pMSJLacI_bead32_cal_corr.mat',...
'PUC306_100pMSJLacI_bead33_cal_corr.mat',...
'PUC306_100pMSJLacI_bead34_cal_corr.mat',...
'PUC306_100pMSJLacI_bead35_cal_corr.mat',...
'PUC306_100pMSJLacI_bead36_cal_corr.mat',...
'PUC306_100pMSJLacI_bead37_cal_corr.mat',...
'PUC306_100pMSJLacI_bead38_cal_corr.mat',...
'PUC306_100pMSJLacI_bead39_cal_corr.mat',...
};

looping_filename={...
{'PUC306_100pMSJLacI_bead30_trjA_corr.mat'},...
{'PUC306_100pMSJLacI_bead31_trjA_corr.mat'},...
{'PUC306_100pMSJLacI_bead32_trjA_corr.mat', 'PUC306_100pMSJLacI_bead32_trjB_corr.mat'},...
{'PUC306_100pMSJLacI_bead33_trjA_corr.mat', 'PUC306_100pMSJLacI_bead33_trjB_corr.mat'},...
{'PUC306_100pMSJLacI_bead34_trjA_corr.mat', 'PUC306_100pMSJLacI_bead34_trjB_corr.mat'},...
{'PUC306_100pMSJLacI_bead35_trjA_corr.mat'},...
{'PUC306_100pMSJLacI_bead36_trjA_corr.mat', 'PUC306_100pMSJLacI_bead36_trjB_corr.mat'},...
{'PUC306_100pMSJLacI_bead37_trjA_corr.mat', 'PUC306_100pMSJLacI_bead37_trjB_corr.mat'},...
{'PUC306_100pMSJLacI_bead38_trjA_corr.mat', 'PUC306_100pMSJLacI_bead38_trjB_corr.mat'},...
{'PUC306_100pMSJLacI_bead39_trjA_corr.mat', 'PUC306_100pMSJLacI_bead39_trjB_corr.mat'},...
};


