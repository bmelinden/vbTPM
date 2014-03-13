% [trjdat,caldat]=VB7_getTrjData(runinputfile,k,b)
%
% retrieve looping and calibration data trajectory k,b for the analysis
% specified in runinput file, i.e., from looping_filename{k}{b} and
% calibration_filename{k}. An options structure generated from a runinpit
% file can also be used as input.

% The objects trjdat,caldat have fields
%
% name      : name of runinput file
% data      : preprocessed data object, as used by the analysis routines of
%             vbTPM.
% x         : position trajectory
% RMS       : RMS trace, filtered as specified in the runinput file
% various additional parameters

% change-log
% ML 2012-06-01 : handle runinput filenames ending with .m
% ML 2014-03-13 : handle options structure input as well

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

if(isstruct(runinputfile))
    opt=runinputfile;
else
    % get settings for this runinput file
    if(strcmp(runinputfile(end-1:end),'.m'))
        runinputfile=runinputfile(1:end-2);
    end
    opt=VB7_getOptions(runinputfile);
end
%eval(runinputfile);
%% calibration data
% determine data file names and load data
caldat=struct;
if(opt.include_calibration)
    caldat.name=[opt.source_path opt.calibration_filename{k}];
    caldat=load(caldat.name);
    % rename the coordinate field to 'x'
    if(~strcmp('x',opt.calibration_xyfield))
        caldat.x=caldat.(opt.calibration_xyfield);
        caldat=rmfield(caldat,opt.calibration_xyfield);
    end
    % driftcorrect by lowpass filter, if specified in the runinputfile
    if(opt.driftcorrection)
        caldat.x=BWdriftcorrect(caldat.x,opt.fCut,opt.fSample);
        caldat.fCut=opt.fCut;
    end
    % preprocessed data
    caldat.data=VB7_preprocess(caldat.x,opt.downSample,opt.fSample);
    % RMS trace
    caldat.RMS=RMSKBgaussfilter(caldat.x,opt.tSigma,opt.fSample);
    caldat.tSigma=opt.tSigma;
    % useful information
    caldat.fSample=opt.fSample;
    caldat.downSample=opt.downSample;
end
trjdat=struct;
if(opt.include_looping)
    trjdat.name=[opt.source_path opt.looping_filename{k}{b}];
    trjdat=load(trjdat.name);
    if(~strcmp('x',opt.looping_xyfield))
        trjdat.x=trj.(opt.looping_xyfield);
        trjdat=rmfield(trjdat,opt.looping_xyfield);
    end
    % driftcorrect by lowpass filter, if specified in the runinputfile
    if(opt.driftcorrection)
        trjdat.x=BWdriftcorrect(trjdat.x,opt.fCut,opt.fSample);
        trjdat.fCut=opt.fCut;
    end

    % preprocessed data
    trjdat.data=VB7_preprocess(trjdat.x,opt.downSample,opt.fSample);
    % RMS trace
    trjdat.RMS=RMSKBgaussfilter(trjdat.x,opt.tSigma,opt.fSample);
    trjdat.tSigma=opt.tSigma;
    % useful information
    trjdat.fSample=opt.fSample;
    trjdat.downSample=opt.downSample;
end
