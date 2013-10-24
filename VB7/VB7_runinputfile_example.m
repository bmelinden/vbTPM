%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_VBEMiter_nomex.m, VBEM iteration without mex files, in the vbTPM package
% =========================================================================
% 
% Copyright (C) 2013 Martin Lindén
% 
% E-mail: bmelinden@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% % Public License for more details.
% 
% Additional permission under GNU GPL version 3 section 7
% 
% If you modify this Program, or any covered work, by linking or combining it
% with Matlab or any Matlab toolbox, the licensors of this Program grant you 
% additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code
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
% M.L. 2011-08-15   : added field names for xy data in calibration and
%                     looping data files.
% M.L. 2011-11-07   : reordered lines to reflect how often they need to be
%                     edited, and tried to add helpful comments.
% M.L. 2012-03-28   : changed some parameter values in line with what we
%                     have actually been using

%%%%%%%% part I: parameters that need updating for every data set %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How to find the data, and where to store the results. 
% Paths are relative to the location of the runinput file
source_path     =   'path_to_data/'; % (local) path to xy data (calibration and looping)
target_path     =   'path_to_results/'; % path to put results in (does not have to exist) 
PIDprefix       =   'myanalysis1'; % prefix for the result filenames

% Data files: the idea here is that each bead gets a 'dry run' (no activity)
% for calibration, and a couple of active runs (perhaps interrupted by
% refocusing etc), and the analysis can keep track of each bead separately. 
% The cell structure pattern must be followed here.
calibration_filename={...
'bead1_calibration.mat',...
'bead2_calibration.mat',...
'bead3_calibration.mat'};

looping_filename={...
{'bead1_trj1.mat','bead1_trj2.mat','bead1_trj3.mat'},...
{'bead2_trj1.mat'},...
{'bead3_trj1.mat','bead3_trj2.mat'}};

% Instead of actually typing the file names here, one could have them saved
% to a file and replace the above statements with load statements, e.g., 
% calibration_filename=load(mycalibrationfilenamefile);
% looping_filename=load(myloopingfilenamefile);

%%%% part II: parameters that probably stays the same througout the whole %
%%%% experiment, or should not be touched at all. %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data properties and preprocesing options %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How the position data is stored, namely the name of xy-position field in
% calibration and looping data files, as in   
% foo=load(filname,nn_xyfield);
% xy=foo.(xy_field);
% Empty strings (nn_xyfield='') means that no field is read, as in 
% xy=load(filname);
calibration_xyfield='x';    
looping_xyfield ='x';    	
% The data should be stored as a T by 2 matrix.

% Properties of the data, and preprocessing options
fSample    =30;     % sampling frequency of the data [Hz]
downSample =3 ;     % downsampling factor (1,2,3,..., 1 means no downsampling will be done).
% Downsampling can save computation time. Downsampling might hurt the fit
% quality, but also seems to alleviate some problems with analyzing very
% long (~100 000 data points) trajectories.
% 
% If you set downSample>1, make sure that the life time of the shortest
% events you need to detect are a lot longer than downSample/fSample.

driftcorrection=false;    % true/false. true = do drift-correction with Butterworth filter
fCut=0.05;                % filter cut-off frequency [Hz]
include_calibration=true; % include calibration traces in the analysis
include_looping=true;     % include looping traces in the analysis

% parameters for the construction of an initial guess %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range of K,B values in which to look for initial guesses. Choose a wider
% interval for KBSrange than for KBGrange. These are the most important
% initial guess parameters, and the most common source of convergence problems. 
KBGrange=[0.3 0.9 0.2e-4 5e-4]; % [Kmim Kmax Bmin Bmax] for the genuine states 
KBSrange=[0 1 0 1e-3];          % [Kmim Kmax Bmin Bmax] for spurious states. 

tSigma=2;   % Gaussian window width [s] to use with RMSKBgaussfilter.
KBbins=500; % number of bins in the KBGrange histogram used for the initial guess
% Experiment with RMSKBgaussfilter and actual data, and make sure that the
% specified intervals include most of the values in each trajectory. 

% Typical dwell times. Guessing on the short side are probably better. 
tau0=5;  % typical dwell times for genuine [s],
tauS=1;  % and spurious, states [s]

wKBT=10000; % weight (pseudo-counts) of initial guess. Do not make too small.

% Prior parameters, used by VB7_priorParent. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are not so important, rough order of magnitude estiamtes
% are good enough. 
K0=0.6;B0=5e-5; % mean values for the prior distributrions of the K and B 
% parameters. Getting the right order of magnitude depends on your sampling
% frequency, length units, tether length etc. The following relations are
% useful for estimates (or use RMSKBgaussfilter on some real data to see
% what typical values come out):
% Bead std, RMS = 1/sqrt(B0*(1-K0^2))
% Bead corner frequency, fc= - fSample*log(K0)
tD=1;tA=5;  % typical dwell time (tD [s]) and total prior strength for transitions (tA [s])

% the rest of the prior parameters are probably good as they are.
fB=[];    % prior strength (pseudocounts) for B parameters. Leave empty to use default value. 
Kstd=0.5; % standard deviation of the prior for the K parameter. Do not choose to small
fPi=[];     % prior strength of initial state. Leave empty to use default value.
% this parameter deals with how the total prior strength depends on the
% number of hidden states.
KBscaling=2;% 1) total emission strength constant 
            % 2) strength of single state constant (REALLY recommended).
% prior parameters for spurious states (Legacy parameters. Ignore.)
fBc=[];Bc0=5e-5;Kc0=0.6;KcStd=0.3;tDc0=10;

% parameters for how to run the analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restarts_trj=30; % number of independent analysis restarts per looping trajectory
Ninit_trj=50;    % maximum number of hidden states
restarts_cal=10; % number of idependent analyses of each calibration trajectory
Ninit_cal=10;    % maximum number of hidden states in calibration trajectories
% The algorithm starts with Ninit_xxx states and systematically tries to
% remove unnecessary ones. Hence, a large number increases the chance to
% get a good fit and might allow fewer restarts. But the computational cost
% might also go up. A rough sanity check is that the number of states
% actually found should be significantly lower than Ninit_xxx, perhaps by a
% factor two.

one_at_a_time=true; % analyze just one trajectory per call to VB7_batch_run, 
% and then quit Matlab (hence, a script that restarts matlab many times is
% needed). If true, the analysis will instead go on until all  data is
% analyzed. This is a workaround for bad memory management in Matlab. 

plotAnything=false; % Make a number of (not very well documented) plots 
% that might be useful for debugging parameter choices.
% debug flags
haveMaxLength=false; % If true, use only the first maxLength data points in 
maxLength=1000;      % each trajectory. A good way to make fast testruns in 
                     % order to check file names etc.


        
