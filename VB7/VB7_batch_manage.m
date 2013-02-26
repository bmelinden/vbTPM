function [resultFile,cal,trj]=VB7_batch_manage(varargin)
% [resultFile,cal,trj]=VB7_batch_manage(runinputfile,a1,a2,...)
%
% manage the results files from the analysis jobs described by
% runinputfile, by performing actions a1,a2,... 
% 
% actions:
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
% M.L. 2011-11-12   : add fields to the trj anc cal structures.
%% call runinput file
resultFile=[];
if(nargin>0)
    runinputfile=varargin{1};
    if(strcmp(runinputfile(end-1:end),'.m'))
        runinputfile=runinputfile(1:end-2);
    end
    eval(runinputfile);
end
if(nargin<=1)
    disp('no action specified')
    return
end
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
%% file names
calSaveName=@(k)(sprintf('%s%s_bead%03d_cal.mat',target_path,PIDprefix,k));
calTmpName =@(k)(sprintf('%s%s_bead%03d_calTMP.mat',target_path,PIDprefix,k));
calLogName =@(k)(sprintf('%s%s_bead%03d_calLOG.txt',target_path,PIDprefix,k));
trjSaveName=@(k,j)(sprintf('%s%s_bead%03dr%03d_trj.mat',target_path,PIDprefix,k,j));
trjTmpName =@(k,j)(sprintf('%s%s_bead%03dr%03d_trjTMP.mat',target_path,PIDprefix,k,j));
trjLogName =@(k,j)(sprintf('%s%s_bead%03dr%03d_trjLOG.txt',target_path,PIDprefix,k,j));
% variables to load and collect
calVariableList={'Wcal','NFcal','NFical','calLogFile','calDataFile','driftcorrection','fCut','calibration_xyfield'};
trjVariableList={'Wtrj','NFtrj','NFitrj','trjLogFile','trjDataFile','driftcorrection','fCut','looping_xyfield'};
cal=cell(1,1);
trj=cell(1,1);
finished=[];
%% loop through misc arguments and decide on actions to take
collect=false;
cleanup=false;
autoclean=false;
minimalSave=true;
saveresults=true;
% default output destination
save_path=target_path;
save_name=[PIDprefix '_result.mat'];

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
        saveresults && true;
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
    
    end
    na=na+1;

end

if(saveresults && isempty(resultFile))
    resultFile=sprintf('%s%s',save_path,save_name);
end
%% run analysis
finished=zeros(size(looping_filename));
for k=1:max(length(looping_filename),length(calibration_filename))
    if(include_calibration)
        calSaveFile=calSaveName(k);
        calLogFile = calLogName(k);
        calTmpFile = calTmpName(k);
        calDataFile=[source_path calibration_filename{k}];
        
        calFileExist=exist(calSaveFile,'file');
        calFileIsFinished=false;
        if(calFileExist)
            caltmp=load(calSaveFile);
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
    if(include_looping)
        trjFileIsFinished=zeros(size(looping_filename{k}));
        for j=1:length(looping_filename{k})
            trjSaveFile=trjSaveName(k,j);
            trjLogFile = trjLogName(k,j);
            trjTmpFile = trjTmpName(k,j);
            trjDataFile=[source_path looping_filename{k}{j}];
            
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

