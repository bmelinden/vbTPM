% VB7_runAnalysis(ri)
% start a new analysis job based on he runinput file ri. Runinput
% parameters are explained in the runinput file examples in example1/ 
%
% At least on Linux, it is good memory management to quit matlab in between
% running consecutive trajectories (using one_at_a_time=true; in the
% runinput file), since Matlab tends to hoard memory while loading and
% clearing long data sets. 
%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_batch_run.m, analysis manager for the vbTPM package
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
function VB7_batch_run(ri)
VB7_printGPL('VB7_batch_run.m')
%% run runinput file
if(strcmp(ri(end-1:end),'.m'))
    ri=ri(1:end-2);
end
opt=VB7_getOptions(ri);
%eval(ri); % replace this by the options struct construct instead

if(~isfield(opt,'KBscaling'));opt.KBscaling=[];end
%% sanity check of source and target paths
if(isempty(opt.source_path))
    opt.source_path='./';
elseif(~strcmp(opt.source_path(end),filesep))
    opt.source_path=[opt.source_path filesep];
end
if(isempty(opt.target_path))
    opt.target_path='./';
elseif(~strcmp(opt.target_path(end),filesep))
    opt.target_path=[opt.target_path filesep];
end
if(~exist(opt.target_path,'dir'))
    mkdir(opt.target_path)
end
%% set some defaults for backward compatibility
if(~isfield(opt,'calibration_xyfield'))
    opt.calibration_xyfield='x';
    warning('VB7_batch_run : using default name (x) for calibration_xyfield. Runinput file outdated?')
end
if(~isfield(opt,'looping_xyfield'))
    opt.looping_xyfield='x';
    warning('VB7_batch_run : using default name (x) for looping_xyfield. Runinput file outdated?')
end
%% file names
calSaveName=@(k)(sprintf('%s%s_bead%03d_cal.mat',opt.target_path,opt.PIDprefix,k));
calTmpName =@(k)(sprintf('%s%s_bead%03d_calTMP.mat',opt.target_path,opt.PIDprefix,k));
calLogName =@(k)(sprintf('%s%s_bead%03d_calLOG.txt',opt.target_path,opt.PIDprefix,k));
trjSaveName=@(k,j)(sprintf('%s%s_bead%03dr%03d_trj.mat',opt.target_path,opt.PIDprefix,k,j));
trjTmpName =@(k,j)(sprintf('%s%s_bead%03dr%03d_trjTMP.mat',opt.target_path,opt.PIDprefix,k,j));
trjLogName =@(k,j)(sprintf('%s%s_bead%03dr%03d_trjLOG.txt',opt.target_path,opt.PIDprefix,k,j));
resultFile=sprintf('%s%s_result.mat',opt.target_path,opt.PIDprefix);
% variables to save
% saving these variables are probably redundant as they are present in the
% opt struct, but until I have checked that they are not used somewhere
% downstream, they stay:
driftcorrection=opt.driftcorrection;fCut=opt.fCut;
calibration_xyfield=opt.calibration_xyfield;looping_xyfield=opt.looping_xyfield;
calVariableList={'Wcal','NFcal','NFical','calLogFile','calDataFile','driftcorrection','fCut','calibration_xyfield','opt'};
trjVariableList={'Wtrj','NFtrj','NFitrj','trjLogFile','trjDataFile','driftcorrection','fCut','looping_xyfield','opt'};
%% run analysis
diary off
analyzed_something=false;
for k=1:max(length(opt.looping_filename),length(opt.calibration_filename))
    % construct prior parameter object
    Wparent=VB7_priorParent(opt.fSample,opt.downSample,...
        opt.fPi,opt.tD,opt.tA,opt.tDc0,opt.fB,opt.B0,opt.fBc,opt.Bc0,opt.K0,opt.Kstd,opt.Kc0,opt.KcStd,opt.KBscaling);
    if(opt.include_calibration)
        calSaveFile=calSaveName(k);
        calLogFile = calLogName(k);
        calTmpFile = calTmpName(k);
        calDataFile=[opt.source_path opt.calibration_filename{k}];
        % see if resultfile is already dealt with
        if(exist(calSaveFile,'file'))
            disp([calSaveFile ' already finished or checked out.'])
        else
            save(calSaveFile,'calLogFile') % reserve file
            diary(calLogFile)
            disp(['checking out ' calSaveFile ])
            disp('-----------------------------------------')
            
            % load data
            x=read_xyfield([opt.source_path opt.calibration_filename{k}],opt.calibration_xyfield);
            % debug
            if(opt.haveMaxLength)
                MT=min(opt.maxLength,size(x,1));
                x=x(1:MT,:);
            end
            if(driftcorrection)
                x=BWdriftcorrect(x,fCut,opt.fSample);
            end
            data=VB7_preprocess(x,opt.downSample,opt.fSample);
            
            % construct initial guess method
            if(opt.plotAnything)
                figure(1)
                clf
                fig=subplot(2,1,1);
                [initFun,init_parameters]=VB7_initialGuess_KBregion(x,...
                    opt.tau0,opt.tauS,opt.wKBT,opt.KBGrange,opt.KBbins,opt.KBSrange,opt.tSigma,...
                    opt.fSample,opt.downSample,fig);
                title(opt.calibration_filename{k})
                pause(0.1)
            else
                [initFun,init_parameters]=VB7_initialGuess_KBregion(x,...
                    opt.tau0,opt.tauS,opt.wKBT,opt.KBGrange,opt.KBbins,opt.KBSrange,opt.tSigma,...
                    opt.fSample,opt.downSample);
            end
            clear x;
            [Wcal,NFcal,NFical]= VB7_analyzeTrace(data,Wparent,...
                initFun,opt.Ninit_cal,opt.restarts_cal,calTmpFile);
            clear data;
            save(calSaveFile,calVariableList{1:end})
            delete(calTmpFile)
            disp(['saved analysis result to ' calSaveFile])
            disp('-----------------------------------------')
            diary off
            analyzed_something=true;
        end
    end
    if(opt.include_looping)
        for j=1:length(opt.looping_filename{k})
            trjSaveFile=trjSaveName(k,j);
            trjLogFile = trjLogName(k,j);
            trjTmpFile = trjTmpName(k,j);
            trjDataFile=[opt.source_path opt.looping_filename{k}{j}];
            % see if resultfile is already dealt with
            if(exist(trjSaveFile,'file'))
                disp([trjSaveFile ' already finished or checked out.'])
            else
                save(trjSaveFile,'trjLogFile') % reserve file
                diary(trjLogFile)
                disp(['checking out ' trjSaveFile ])
                disp('-----------------------------------------')                                
                % load data
                x=read_xyfield(trjDataFile,opt.looping_xyfield);
                if(size(x,2) ~=2)
                    error('VB7_batch_run: x must be a 2 by T matrix')
                end
                % debug
                if(opt.haveMaxLength)
                    MT=min(opt.maxLength,size(x,1));
                    x=x(1:MT,:);
                end
                if(driftcorrection)
                    x=BWdriftcorrect(x,fCut,opt.fSample);
                end
                data=VB7_preprocess(x,opt.downSample,opt.fSample);
                
                % construct initial guess method
                if(opt.plotAnything)
                    fig=subplot(2,1,2);
                    [initFun,init_parameters]=VB7_initialGuess_KBregion(...
                        x,opt.tau0,opt.tauS,opt.wKBT,opt.KBGrange,opt.KBbins,opt.KBSrange,opt.tSigma,...
                        opt.fSample,opt.downSample,fig);
                    title(opt.looping_filename{k}{j})
                    pause(0.1)
                else
                    [initFun,init_parameters]=VB7_initialGuess_KBregion(...
                        x,opt.tau0,opt.tauS,opt.wKBT,opt.KBGrange,opt.KBbins,opt.KBSrange,opt.tSigma,...
                        opt.fSample,opt.downSample);
                end
                clear x;
                [Wtrj,NFtrj,NFitrj]= ...
                    VB7_analyzeTrace(data,Wparent,initFun,opt.Ninit_trj,opt.restarts_trj,trjTmpFile);
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
    if(isfield(opt,'one_at_a_time') && opt.one_at_a_time && analyzed_something)
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
