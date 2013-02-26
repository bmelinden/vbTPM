function TotalResultFile=VB7_pseudoparallel_manage(varargin)
% TotalResultFile=VB7_cleanup_pseudoparallel(runinputfile,opt)
%
% manage output from VB7_pseudoparallel_run
%
%
% runinputfile  : runinput file to manage
% opt           : list of actions to take (in order). Possible choices:
%
% reset checkout :  Mark checked-out but not finished beads as not checked
%                   out, making them available for analysis on subsequent
%                   calls to VB7_pseudoparallel_run. (Warning: even beads
%                   that are still running are affected, and hence might be
%                   analyzed several times.)
% resultfiles   : check if all finished jobs have existing results files,
%                 and update progress record for those that do not.
% collect       : collect all finished analyses in a single result file 
%                 ~_result.mat

% removes entries in the progress file that has only been checked out, but
% not checkid in, in case a job crashes or is interupted before completion.

%% change-log
% M.L. 2011-02-07
% M.L. 2011-02-16   : added output filename to total results

%% read parameters and locate progress file
TotalResultFile=[];
if(nargin>0)
    runinputfile=varargin{1};
    if(strcmp(runinputfile(end-1:end),'.m'))
        runinputfile=runinputfile(1:end-2);
    end
    eval(runinputfile);
    progress_filname=[target_path PIDprefix '_progress.mat'];
    if(~exist(progress_filname,'file'))
        error([progress_filname ' does not exist.'])
    end
end
if(nargin<=1)
    disp('no action specified')
    return
end
%% loop through misc arguments and decide on actions to take
resetCheckout=false;
checkResultsFiles=false;
collectResults=false;
for na=2:nargin
    if(strcmp(varargin{na},'reset checkout'))
        resetCheckout=true;
    end
    if(strcmp(varargin{na},'resultfiles'))
        checkResultsFiles=true;
    end
    if(strcmp(varargin{na},'collect'))
        collectResults=true;
        checkResultsFiles=true;
    end
end

if(checkResultsFiles) % check that all PID result files are present
    PIDfilesPresent=true;
    PF=load(progress_filname);
    for k=1:length(PF.finished)
        if(PF.finished(k)>0)
            PIDresults  =sprintf('%s%s_results%03d.mat',target_path,PIDprefix,PF.finished(k));
            if(~exist(PIDresults,'file'));
                disp(['Could not find results file ' PIDresults])
                PIDfilesPresent=false;
                PF.finished(k)=0;
            end
        end
    end
    if(~PIDfilesPresent)
        save(progress_filname,'-struct','PF');
    end
    clear PF;
end
if(resetCheckout) % reset started record for jobs that are not yet finished
    PF=load(progress_filname);
    ind=find(PF.finished==0);
    PF.started(ind)=0;
    save(progress_filname,'-struct','PF');
    if(isempty(ind))
        disp(['no unfinished jobs in ' runinputfile ])
    else
        disp(['unfinished jobs in ' runinputfile ' : ' int2str(ind)])
    end
    clear PF;
end
if(collectResults)  %% collect all data into one result file
    PF=load(progress_filname);
    ind=find(PF.finished>0);
    if(isempty(ind))
        disp(['Found no results for  '  runinputfile])
        return
    end
    % read first result file
    PID=PF.finished(ind(1));
    PIDresults1  =sprintf('%s%s_results%03d.mat',target_path,PIDprefix,PID);
    R=load(PIDresults1);
    disp(['read result file ' PIDresults1])
    R.finished=ind;
    % merge all subsequent results
    for k=2:length(ind)
        PID=PF.finished(ind(k));
        PIDresults  =sprintf('%s%s_results%03d.mat',target_path,PIDprefix,PID);
        D=load(PIDresults);
        disp(['read result file ' PIDresults])
        pids=find(PF.finished==PID);
        if(D.include_calibration)
            for kk=pids
                R.Wcal{kk}=D.Wcal{kk};
                R.NFcal{kk}=D.NFcal{kk};
                R.NFical{kk}=D.NFical{kk};
            end
        end
        if(D.include_looping)
            for kk=pids
                R.Wtrj{kk}{1}   =   D.Wtrj{kk}{1};
                R.NFtrj{kk}{1}  =  D.NFtrj{kk}{1};
                R.NFitrj{kk}{1} = D.NFitrj{kk}{1};
                for j=2:length(D.Wtrj{kk})
                    if(~isempty(D.Wtrj{kk}{j}))
                        R.Wtrj{kk}{j}   =   D.Wtrj{kk}{j};
                        R.NFtrj{kk}{j}  =  D.NFtrj{kk}{j};
                        R.NFitrj{kk}{j} = D.NFitrj{kk}{j};
                    end
                end
            end
        end
    end
    % write to file
    TotalResultFile=sprintf('%s%s_result.mat',target_path,PIDprefix);
    disp(['Saving analysis results from ' int2str(sum(PF.finished>0)) ...
        ' beads to ' TotalResultFile])
    save(TotalResultFile,'-struct','R');
    clear PF;
end
end



