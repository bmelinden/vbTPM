% [resultFile,cal,trj]=VB7_batch_manage(runinputfile,a1,a2,...)
%
% manage the results files from the analysis jobs described by
% runinputfile, by performing actions a1,a2,... 
% 
% actions:
% 'reconverge'      : reconverge all models before saving results, and also
%                     rewrite the reconverged models (good for small
%                     parameter changes or bugfixes).
% 'collect'         : read analysis results and collect them in cell
%                     vectors that are saved to the resultFile
% 'cleanup'         : look through the list of result files, and delete those
%                     that does not contain finished results, plus assorted
%                     log files and temporary files. Useful if something
%                     crashes and the analysis needs to be restarted.
% 'autoclean'       : do not ask before deleting files in the cleanup
%                     action (use w caution...).
% 'save',bool       : save results to file or not (default true)
% 'file',string     : file name to write results to. [] for 
%                     default name, constructed from runinput file.
% 'path','string'   : path to write results to. [] default is target_path,
%                     from the runinput file. 
% 'saveall'         : save everything, not only the calibration and looping
%                     structures, to the resultFile. For debugging.
% 'displaylevel',d  : controls the amount of information output
%                     1 : nothing (only errors and warnings)
%                     2 : 1+summary (default)   
%                     3 : 2+unfinished beads 
%                     4 : everything

%% change-log
% M.L. 2011-03-01   : started
% M.L. 2011-03-11   : added save options
% M.L. 2011-06-08   : corrected (?) the counter for unanalyzed files.
% M.L. 2011-07-25   : corrected the auto-clean function to actually remove
%                     unfinished result files
% M.L. 2011-10-28   : control the amount of output
% M.L. 2011-11-12   : add fields to the trj and cal structures.
% M.L. 2013-04-25   : option to reconverge models finished before collecting

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_batch_manage.m, manage calculations and analysis results,
% in the vbTPM package
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
function [resultFile,cal,trj]=VB7_batch_manage(varargin)
VB7_printGPL('VB7_batch_manage.m')
%% call runinput file
resultFile=[];
if(nargin>0)
    runinputfile=varargin{1};
    opt=VB7_getOptions(runinputfile);
    %eval(runinputfile);
end
if(nargin<=1)
    disp('no action specified')
    return
end
if(~exist('KBscaling','var'));KBscaling=[];end
%% file names - same convetion as used by VB7_batch_run
calSaveName=@(k)(sprintf('%s%s_bead%03d_cal.mat',opt.target_path,opt.PIDprefix,k));
calTmpName =@(k)(sprintf('%s%s_bead%03d_calTMP.mat',opt.target_path,opt.PIDprefix,k));
calLogName =@(k)(sprintf('%s%s_bead%03d_calLOG.txt',opt.target_path,opt.PIDprefix,k));
trjSaveName=@(k,j)(sprintf('%s%s_bead%03dr%03d_trj.mat',opt.target_path,opt.PIDprefix,k,j));
trjTmpName =@(k,j)(sprintf('%s%s_bead%03dr%03d_trjTMP.mat',opt.target_path,opt.PIDprefix,k,j));
trjLogName =@(k,j)(sprintf('%s%s_bead%03dr%03d_trjLOG.txt',opt.target_path,opt.PIDprefix,k,j));
% variables to load and collect
% saving these variables are probably redundant as they are present in the
% opt struct, but are needed downstream for compatibility: removing them
% consistently from the analysis might make my old analysis files
% unreadable. Hence, this is left to update in a future version
driftcorrection=opt.driftcorrection;fCut=opt.fCut;
calibration_xyfield=opt.calibration_xyfield;looping_xyfield=opt.looping_xyfield;
calVariableList={'Wcal','NFcal','NFical','calLogFile','calDataFile','driftcorrection','fCut','calibration_xyfield','opt'};
trjVariableList={'Wtrj','NFtrj','NFitrj','trjLogFile','trjDataFile','driftcorrection','fCut','looping_xyfield',    'opt'};
cal=cell(1,1);
trj=cell(1,1);
finished=[];
%% loop through misc arguments and decide on actions to take
collect=false;
cleanup=false;
autoclean=false;
minimalSave=true;
saveresults=true;
reconverge=false;
% default output destination
save_path=opt.target_path;
save_name=[opt.PIDprefix '_result.mat'];

displayFinished=false;
displayUnfinished=false;
displaySummary=true;


na=2;
while(na<=nargin)
    parameter=varargin{na};
    if(strcmp(parameter,'collect'))
        collect=true;
    elseif(strcmp(parameter,'cleanup'))
        cleanup=true;
    elseif(strcmp(parameter,'autoclean'))
        cleanup=true;
        autoclean=true;
    elseif(strcmp(parameter,'saveall'))
        minimalSave=false;
    elseif(strcmp(parameter,'save'))
        na=na+1;
        saveresults=varargin{na};        
        saveresults && true; % error if saveresults is a boolean
    elseif(strcmp(parameter,'file'))
        na=na+1;
        save_name=varargin{na};
    elseif(strcmp(parameter,'path'))
        na=na+1;
        save_path=varargin{na};
        if(save_path(end)~=filesep)
            save_path=[save_path filesep];
        end
        elseif(strcmp(parameter,'displaylevel'))
        na=na+1;
        displaylevel=varargin{na};
        if(~isnumeric(displaylevel)||round(displaylevel)~=displaylevel)
            error(['VB7_batch_manage : ' varargin{na} ' not a valid displaylevel'])
        end
            switch displaylevel
                case 1
                    displayFinished=false;
                    displayUnfinished=false;
                    displaySummary=false;
                case 2
                    displayFinished=false;
                    displayUnfinished=false;
                    displaySummary=true;
                case 3
                    displayFinished=false;
                    displayUnfinished=true;
                    displaySummary=true;
                case 4
                    displayFinished=true;
                    displayUnfinished=true;
                    displaySummary=true;
                otherwise
                    error(['VB7_batch_manage : ' varargin{na} ' not a valid displaylevel'])
            end
    elseif(strcmp(parameter,'reconverge'))
        reconverge=true;
    end
    na=na+1;

end

if(saveresults && isempty(resultFile))
    resultFile=sprintf('%s%s',save_path,save_name);
end
%% run analysis
finished=zeros(size(opt.looping_filename));
for k=1:max(length(opt.looping_filename),length(opt.calibration_filename))
    if(opt.include_calibration)
        calSaveFile=calSaveName(k);
        calLogFile = calLogName(k);
        calTmpFile = calTmpName(k);
        calDataFile=[opt.source_path opt.calibration_filename{k}];
        
        calFileExist=exist(calSaveFile,'file');
        calFileIsFinished=false;
        if(calFileExist)
            caltmp=load(calSaveFile); % load ananlysis of single calibration trajectory
            calfields=fieldnames(caltmp);
            if(length(calfields)>1) % then analysis is finished
                calFileIsFinished=true;
            end
        end        
        if(calFileExist && cleanup && ~calFileIsFinished)
            dodelete=false;
            if(autoclean)
                dodelete=true;
            else
                f=input(['delete files ' calSaveFile ' etc? (y/N) : '],'s');
                if(strcmp(f,'Y') || strcmp(f,'y'))
                    dodelete=true;
                end
            end
            if(dodelete)
                delete(calSaveFile);
                delete(calLogFile);
                delete(calTmpFile);
            end
        end
        if(calFileExist && reconverge && calFileIsFinished)
            [~,caldat]=VB7_getTrjData(runinputfile,k,1); % possible error if trj k bead 1 does not exist.
            caltmp.Wcal=VB7iterator(caltmp.Wcal,caldat.data);
            caltmp.calDataFile=calDataFile;
            caltmp.calLogFile=calLogFile;
            save(calSaveFile,'-struct','caltmp');
            disp(['wrote reconverged calibration model to ' calSaveFile])
        end
        if(calFileExist && collect && calFileIsFinished)
            cal{1,k}=load(calSaveFile);
            % check that all fields are present, and try to add them if
            % not
            for mm=1:length(calVariableList)
                if(~isfield(cal{k},calVariableList{mm}) && exist(calVariableList{mm},'var'))
                    cal{k}.(calVariableList{mm})=eval(calVariableList{mm});
                end
            end
            
        end
        if(calFileIsFinished)
            if(displayFinished)
                disp([calDataFile ' : OK'])
            end
        else
            if(displayUnfinished)
                disp([calDataFile ' : not finished'])
            end
        end
    end
    if(opt.include_looping)
        trjFileIsFinished=zeros(size(opt.looping_filename{k}));
        for j=1:length(opt.looping_filename{k})
            trjSaveFile=trjSaveName(k,j);
            trjLogFile = trjLogName(k,j);
            trjTmpFile = trjTmpName(k,j);
            trjDataFile=[opt.source_path opt.looping_filename{k}{j}];
            
            trjFileExist=exist(trjSaveFile,'file');            
            if(trjFileExist)
                trjtmp=load(trjSaveFile);
                trjfields=fieldnames(trjtmp);
                if(length(trjfields)>1) % then this is not yet analyzed
                    trjFileIsFinished(j)=true;                    
                end
            end
            if(trjFileExist && cleanup && ~trjFileIsFinished(j))
                dodelete=false;
                if(autoclean)
                    dodelete=true;
                else
                    f=input(['delete files ' trjSaveFile ' etc? (y/N) : '],'s');
                    if(strcmp(f,'Y') || strcmp(f,'y'))
                        dodelete=true;
                    end
                end
                if(dodelete)
                    delete(trjSaveFile);
                    delete(trjLogFile);
                    delete(trjTmpFile);
                end
            end
            if(trjFileExist && reconverge && trjFileIsFinished(j))
                [trjdat,~]=VB7_getTrjData(runinputfile,k,j); % possible error if trj k bead 1 does not exist.
                trjtmp.Wtrj=VB7iterator(trjtmp.Wtrj,trjdat.data);
                trjtmp.trjDataFile=trjDataFile;
                trjtmp.trjLogFile=trjLogFile;                
                save(trjSaveFile,'-struct','trjtmp');
            disp(['wrote reconverged looping model to ' trjSaveFile])
            end
            if(trjFileExist && collect && trjFileIsFinished(j))
                trj{k}{j}=load(trjSaveFile);
                % check that all fields are present, and try to add them if
                % not
                for mm=1:length(trjVariableList)
                   if(~isfield(trj{k}{j},trjVariableList{mm}) && exist(trjVariableList{mm},'var'))
                        trj{k}{j}.(trjVariableList{mm})=eval(trjVariableList{mm});
                   end
                end
                trj{k}{j}.calk=k; % calibration trajectory is cal{k}
            else
                trj{k}{j}=[];
            end
            if(trjFileIsFinished(j))
                if(displayFinished)
                    disp([trjDataFile ' : OK'])
                end
            else
                if(displayUnfinished)
                    disp([trjDataFile ' : not finished'])
                end
            end
        end
        finished(k)=(sum(trjFileIsFinished==true)==length(trjFileIsFinished)); % check if all trj of this bead were analyzed
    end
end
if(displaySummary)
    disp(['Found ' int2str(sum(finished==true)) ' finished beads (of ' int2str(length(finished)) ')'])
end
if(saveresults)
    if(displaySummary)
        disp('saving results . . . ')
    end
    if(~minimalSave)
        % remove temporary variables created by this script, save
        % everything else
        clear caltmp calFileExist calLogName calSaveName calTmpName calFileIsFinished calDataFile calTmpFile calLogFile calSaveFile calVariableList  calfields calDataFile
        clear trjtmp trjFileExist trjLogName trjSaveName trjTmpName trjFileIsFinished trjDataFile trjTmpFile trjLogFile trjSaveFile trjVariableList  trjfields trjDataFile        
        save(resultFile)
    else
        save(resultFile,'cal','trj');
    end
end

