% VB7_runAnalysis(ri)
% start a new analysis job in pseudoparallel mode, based on he runinput
% file ri. This version only runs a single pass through the file list, and
% thus needs to be called multiple times.
% At least on Linux, it is good memory management to quit matlab in
% between, or to call the pack command, since Matlab tends to hoard memory
% while loading and clearing long data sets.

%% change-log
% M.L. 2011-02-09   : added sanity check on source and target paths, and
%                     runinput file name
% M.L. 2011-02-11   : single run option one_at_a_time in runinput file
% M.L. 2011-02-25   : make sure the prior option KBscaling exists
% M.L. 2011-03-01   : improved naming and saving scheme: savename based on
%                     data file, only save after actually analyzing
%                     something, get rid of the progress file
% M.L. 2011-08-15   : added the possibility to set the field name of the
%                     xy-position array, using variables
%                     calibration_xyfield and looping_xyfield in the
%                     runinput file. 
% M.L. 2011-10-25   : added driftcorrection flag and cut-off frequency to
%                     the list of saved variables, so that complete
%                     reproduction is possible.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_batch_run.m, analysis manager for the vbTPM package
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
function VB7_batch_run(ri)
VB7_printGPL('VB7_batch_run.m')

%% run runinput file
if(strcmp(ri(end-1:end),'.m'))
    ri=ri(1:end-2);
end
eval(ri);
if(~exist('KBscaling','var'));KBscaling=[];end
%% sanity check of source and target paths
if(isempty(source_path))
    source_path='./';
elseif(~strcmp(source_path(end),filesep))
    source_path=[source_path filesep];
end
if(isempty(target_path))
    target_path='./';
elseif(~strcmp(target_path(end),filesep))
    target_path=[target_path filesep];
end
mkdir(target_path)
%% set some defaults for backward compatibility
if(~exist('calibration_xyfield','var'))
    calibration_xyfield='x';
    warning('VB7_batch_run : using default name (x) for calibration_xyfield. Runinput file outdated?')
end
if(~exist('looping_xyfield','var'))
    looping_xyfield='x';
    warning('VB7_batch_run : using default name (x) for looping_xyfield. Runinput file outdated?')
end

%% file names
calSaveName=@(k)(sprintf('%s%s_bead%03d_cal.mat',target_path,PIDprefix,k));
calTmpName =@(k)(sprintf('%s%s_bead%03d_calTMP.mat',target_path,PIDprefix,k));
calLogName =@(k)(sprintf('%s%s_bead%03d_calLOG.txt',target_path,PIDprefix,k));
trjSaveName=@(k,j)(sprintf('%s%s_bead%03dr%03d_trj.mat',target_path,PIDprefix,k,j));
trjTmpName =@(k,j)(sprintf('%s%s_bead%03dr%03d_trjTMP.mat',target_path,PIDprefix,k,j));
trjLogName =@(k,j)(sprintf('%s%s_bead%03dr%03d_trjLOG.txt',target_path,PIDprefix,k,j));
resultFile=sprintf('%s%s_result.mat',target_path,PIDprefix);
% variables to save
calVariableList={'Wcal','NFcal','NFical','calLogFile','calDataFile','driftcorrection','fCut','calibration_xyfield'};
trjVariableList={'Wtrj','NFtrj','NFitrj','trjLogFile','trjDataFile','driftcorrection','fCut','looping_xyfield'};
%% run analysis
diary off
analyzed_something=false;
for k=1:max(length(looping_filename),length(calibration_filename))
    % construct prior parameter object
    Wparent=VB7_priorParent(fSample,downSample,...
        fPi,tD,tA,tDc0,fB,B0,fBc,Bc0,K0,Kstd,Kc0,KcStd,KBscaling);
    if(include_calibration)
        calSaveFile=calSaveName(k);
        calLogFile = calLogName(k);
        calTmpFile = calTmpName(k);
        calDataFile=[source_path calibration_filename{k}];
        % see if resultfile is already dealt with
        if(exist(calSaveFile,'file'))
            disp([calSaveFile ' already finished or checked out.'])
        else
            save(calSaveFile,'calLogFile') % reserve file
            diary(calLogFile)
            disp(['checking out ' calSaveFile ])
            disp('-----------------------------------------')
            
            % load data
            x=read_xyfield([source_path calibration_filename{k}],calibration_xyfield);
            % debug
            if(haveMaxLength)
                MT=min(maxLength,size(x,1));
                x=x(1:MT,:);
            end
            if(driftcorrection)
                x=BWdriftcorrect(x,fCut,fSample);
            end
            data=VB7_preprocess(x,downSample,fSample);
            
            % construct initial guess method
            if(plotAnything)
                figure(1)
                clf
                fig=subplot(2,1,1);
                [initFun,init_parameters]=VB7_initialGuess_KBregion(x,...
                    tau0,tauS,wKBT,KBGrange,KBbins,KBSrange,tSigma,...
                    fSample,downSample,fig);
                title(calibration_filename{k})
                pause(0.1)
            else
                [initFun,init_parameters]=VB7_initialGuess_KBregion(x,...
                    tau0,tauS,wKBT,KBGrange,KBbins,KBSrange,tSigma,...
                    fSample,downSample);
            end
            clear x;
            [Wcal,NFcal,NFical]= VB7_analyzeTrace(data,Wparent,...
                initFun,Ninit_cal,restarts_cal,calTmpFile);
            clear data;
            save(calSaveFile,calVariableList{1:end})
            delete(calTmpFile)
            disp(['saved analysis result to ' calSaveFile])
            disp('-----------------------------------------')
            diary off
            analyzed_something=true;
        end
    end
    if(include_looping)
        for j=1:length(looping_filename{k})
            trjSaveFile=trjSaveName(k,j);
            trjLogFile = trjLogName(k,j);
            trjTmpFile = trjTmpName(k,j);
            trjDataFile=[source_path looping_filename{k}{j}];
            % see if resultfile is already dealt with
            if(exist(trjSaveFile,'file'))
                disp([trjSaveFile ' already finished or checked out.'])
            else
                save(trjSaveFile,'trjLogFile') % reserve file
                diary(trjLogFile)
                disp(['checking out ' trjSaveFile ])
                disp('-----------------------------------------')                                
                % load data
                x=read_xyfield(trjDataFile,looping_xyfield);
                if(size(x,2) ~=2)
                    error('VB7_batch_run: x must be a 2 by T matrix')
                end
                % debug
                if(haveMaxLength)
                    MT=min(maxLength,size(x,1));
                    x=x(1:MT,:);
                end
                if(driftcorrection)
                    x=BWdriftcorrect(x,fCut,fSample);
                end
                data=VB7_preprocess(x,downSample,fSample);
                
                % construct initial guess method
                if(plotAnything)
                    fig=subplot(2,1,2);
                    [initFun,init_parameters]=VB7_initialGuess_KBregion(...
                        x,tau0,tauS,wKBT,KBGrange,KBbins,KBSrange,tSigma,...
                        fSample,downSample,fig);
                    title(looping_filename{k}{j})
                    pause(0.1)
                else
                    [initFun,init_parameters]=VB7_initialGuess_KBregion(...
                        x,tau0,tauS,wKBT,KBGrange,KBbins,KBSrange,tSigma,...
                        fSample,downSample);
                end
                clear x;
                [Wtrj,NFtrj,NFitrj]= ...
                    VB7_analyzeTrace(data,Wparent,initFun,Ninit_trj,restarts_trj,trjTmpFile);
                clear data;
                save(trjSaveFile,trjVariableList{1:end})
                delete(trjTmpFile)
                disp(['saved analysis result to ' trjSaveFile])
                disp('-----------------------------------------')
                diary off
                analyzed_something=true;
            end
        end
    end
    if(exist('one_at_a_time','var') && one_at_a_time && analyzed_something)
        break
    end
end
end
function x=read_xyfield(file,fieldname)
% x=VB7_read_xydata(file,fieldname)
% Reads the field fieldname from file file, and return in the variable x.
% Default: fieldname='x'. If fieldname is empty, the content of the whole
% file is returned. 

%% actual code
if(~isempty(fieldname))
    foo=load(file,fieldname);
    x=foo.(fieldname);
    clear foo
else
    x=load(file);
end
% check that x has correct dimensions
if(size(x,2) ~=2)
    error('VB7_batch_run:read_xyfield: xy-positions must be a T by 2 matrix.')
end

end
