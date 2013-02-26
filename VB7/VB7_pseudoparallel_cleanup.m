function VB7_pseudoparallel_cleanup(runinputfile)
% VB7_cleanup_pseudoparallel(runinputfile)
% removes entries in the progress file that has only been checked out, but
% not checkid in, in case a job crashes or is interupted before completion.

%% change-log
% M.L. 2011-02-07

%% code

eval(runinputfile);
progress_filname=[target_path PIDprefix '_progress.mat'];
if(~exist(progress_filname,'file'))
    error([progress_filname ' does not exist.'])
end
for attempt=1:20
    try
        PF=load(progress_filname);
        ind=find(PF.finished==0);
        PF.started(ind)=0;
        save(progress_filname,'-struct','PF');
        if(isempty(ind))
            disp(['no unfinished jobs in ' runinputfile ])
        else
            disp(['unfinished jobs in ' runinputfile ' : ' int2str(ind)])
        end
        break
    catch me
        disp(['Failed to load or save ' progress_filname])
        disp('... will try again in a little while.')
        me.getReport
        pause(5+5*rand)
    end
end
if(attempt==20)
    error(['Failed to read ' progress_filname ...
        ' too many times.'])
end

%% check that all PID files are present
allfound=true;
for k=1:length(PF.finished)
    if(PF.finished(k)>0)
        PIDresults  =sprintf('%s%s_results%03d.mat',target_path,PIDprefix,PF.finished(k));
        if(~exist(PIDresults,'file'));
            disp(['Could not find results file ' PIDresults])
            allfound=false;
            PF.finished(k)=0;
        end
    end
end
if(~allfound)
    save(progress_filname,'-struct','PF');
    VB7_pseudoparallel_cleanup(runinputfile);
    return
end
%% collect all data into one result file
ff=dir([target_path PIDprefix '_results*.mat']);
PID0result=[target_path ff(1).name];
R=load(PID0result);
PID0=R.PID;
if(isempty(find(PF.finished==0,1))) % then collect it in a total results file
    PID0result=sprintf('%s%s_result.mat',target_path,PIDprefix);
end

for k=1:length(PF.finished)    
    PID=PF.finished(k);
    if(PID==PID0 || PF.finished(k)==0)
        continue
    end
    PIDresults  =sprintf('%s%s_results%03d.mat',target_path,PIDprefix,PID);
    D=load(PIDresults);
    
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
    PF.started(pids)=PID0;
    PF.finished(pids)=PID0;
    save(PID0result,'-struct','R');
    save(progress_filname,'-struct','PF');
    disp([PIDresults ' ---> ' PID0result])
    disp(['delete(' PIDresults ')']);
    delete(PIDresults)
end



