function VB7_pseudoparallel_run(ri)
% VB7_pseudoparallel_run(ri)
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

%% run runinput file
if(strcmp(ri(end-1:end),'.m'))
    ri=ri(1:end-2);
end
eval(ri);
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
%% part III: analysis
%% construct results data structure
if(include_looping)
    Wtrj=cell(size(looping_filename));
    NFtrj=cell(size(looping_filename));
    NFitrj=cell(size(looping_filename));
end
if(include_calibration)
    Wcal=cell(size(calibration_filename));
    NFcal=cell(size(calibration_filename));
    NFical=cell(size(calibration_filename));    
end
%% create progess file if not present
progress_filname=[target_path PIDprefix '_progress.mat'];
if(~exist(progress_filname,'file'))
    started=zeros(1,length(looping_filename));
    finished=zeros(1,length(looping_filename));
    save(progress_filname,'started','finished')
    clear started finished
end
%% get process ID number and create PID file names
ff=dir([target_path PIDprefix '_results*.mat']);
PID=1;
for k=1:length(ff)
   PID=max(PID,1+str2double(ff(k).name(end-6:end-4)));
end
disp(['process ID ' int2str(PID)])
PIDresults  =sprintf('%s%s_results%03d.mat',target_path,PIDprefix,PID);
PIDtmp      =sprintf('%s%s_tmp%03d.mat',target_path,PIDprefix,PID);
PIDlog      =sprintf('%s%s_log%03d.txt',target_path,PIDprefix,PID);
save(PIDresults) 
%% run analysis 
diary off
diary(PIDlog)
for k=1:length(looping_filename)
    % reserve a file number for analysis (for multiple CPUs running the same
    % data set). Try to deal with conflicts by waiting some time in case of
    % load/save errors
    found_assignment=false;
    clear PF
    for attempt=1:20
        try
            PF=load(progress_filname);
            if(PF.started(k)==0 && PF.finished(k)==0)
                PF.started(k)=PID;
                save(progress_filname,'-struct','PF');
                found_assignment=true;
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
    if(~found_assignment)
        pause(1+rand);
        if(PF.started>0)
            disp([ri ': Bead no. ' int2str(k) ' already checked out.'])
        elseif(PF.finished(k)>0)
            disp([ri ': Bead no. ' int2str(k) ' already finished.'])
        end            
        continue % go on to look at the next one
    else
        clear PF;
        %% start analyzing this bead
        if(~exist('KBscaling','var'));KBscaling=[];end
        disp([ri ' : Checked out bead no. ' int2str(k) ', in ' progress_filname])
        % construct prior parameter object
        Wparent=VB7_priorParent(fSample,downSample,...
            fPi,tD,tA,tDc0,fB,B0,fBc,Bc0,K0,Kstd,Kc0,KcStd,KBscaling);
        if(include_calibration)
            % load data ML: this should call a wrapper function
            % x=xloader(filename), to allow loading different formats
            load([source_path calibration_filename{k}],'x');            
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
            [Wcal{k},NFcal{k},NFical{k}]= VB7_analyzeTrace(data,Wparent,initFun,Ninit_cal,restarts_cal,PIDtmp);
            clear data;
            save(PIDresults)
        end
        if(include_looping)
            for j=1:length(looping_filename{k})
                % load data ML: this should call a wrapper function
                % x=xloader(filename), to allow loading different formats
                load([source_path looping_filename{k}{j}],'x');
            
                % debug
                if(haveMaxLength)
                    MT=min(maxLength,size(x,1));
                    x=x(1:MT,:);
                end
                if(driftcorrection)
                    x=BWdriftcorrect(x,fCut,fSample);
                end               
                data=VB7_preprocess(x,downSample,fSample,[source_path looping_filename{k}{j}]);
                
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
                [Wtrj{k}{j},NFtrj{k}{j},NFitrj{k}{j}]= ...
                    VB7_analyzeTrace(data,Wparent,initFun,Ninit_trj,restarts_trj,PIDtmp);
                clear data;
                save(PIDresults)
            end
        end        
        %% check in result
        for attempt=1:20
            try
                PF=load(progress_filname);
                PF.finished(k)=PID;
                save(progress_filname,'-struct','PF');
                clear PF;
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
        save(PIDresults);   % save results
        delete(PIDtmp);     % delete temporary file
        disp('--------------------------------------')
        if(exist('one_at_a_time','var') && one_at_a_time==true)
            break
        end
    end
end
diary off
end
