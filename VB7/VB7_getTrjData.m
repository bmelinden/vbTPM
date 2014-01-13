% [trjdat,caldat]=VB7_getTrjData(runinputfile,k,b)
%
% retrieve looping and calibration data trajectory k,b for the analysis
% specified in runinput file. The objects trjdat,caldat have fields
%
% name      : name of data file
% data      : preprocessed data object, as used in the analysis
% x         : position trajectory
% RMS       : RMS trace
%

% change-log
% ML 2012-06-01 : handle runinput filenames ending with .m

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_getTrjData.m, access analysis results in the vbTPM package
% =========================================================================
% 
% Copyright (C) 2014 Martin Lindén
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
function [trjdat,caldat]=VB7_getTrjData(runinputfile,k,b)

% get settings for this runinput file
if(strcmp(runinputfile(end-1:end),'.m'))
    runinputfile=runinputfile(1:end-2);
end
eval(runinputfile);
% add file separator to source path if absent
if(~strcmp(source_path(end),filesep))
    source_path=[source_path filesep]
end

%% calibration data
% determine data file names and load data
caldat=struct;
if(include_calibration)
    caldat.name=[source_path calibration_filename{k}];
    caldat=load(caldat.name);
    % rename the coordinate field to 'x'
    if(~strcmp('x',calibration_xyfield))
        caldat.x=caldat.(calibration_xyfield);
        caldat=rmfield(caldat,calibration_xyfield);
    end
    % driftcorrect by lowpass filter, if specified in the runinputfile
    if(driftcorrection)
        caldat.x=BWdriftcorrect(caldat.x,fCut,fSample);
        caldat.fCut=fCut;
    end
    % preprocessed data
    caldat.data=VB7_preprocess(caldat.x,downSample,fSample);
    % RMS trace
    caldat.RMS=RMSKBgaussfilter(caldat.x,tSigma,fSample);
    caldat.tSigma=tSigma;
    % useful information
    caldat.fSample=fSample;
    caldat.downSample=downSample;
end
trjdat=struct;
if(include_looping)
    trjdat.name=[source_path looping_filename{k}{b}];
    trjdat=load(trjdat.name);
    if(~strcmp('x',looping_xyfield))
        trjdat.x=trj.(looping_xyfield);
        trjdat=rmfield(trjdat,looping_xyfield);
    end
    % driftcorrect by lowpass filter, if specified in the runinputfile
    if(driftcorrection)
        trjdat.x=BWdriftcorrect(trjdat.x,fCut,fSample);
        trjdat.fCut=fCut;
    end

    % preprocessed data
    trjdat.data=VB7_preprocess(trjdat.x,downSample,fSample);
    % RMS trace
    trjdat.RMS=RMSKBgaussfilter(trjdat.x,tSigma,fSample);
    trjdat.tSigma=tSigma;
    % useful information
    trjdat.fSample=fSample;
    trjdat.downSample=downSample;
end


