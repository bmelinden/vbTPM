function [Wbest,NF,NFiMax,rndState0]=VB7_analyzeTrace(data,Wparent,initFun,N0,restarts,savename,rndState0)
% [Wbest,NF,NFiMax,rndState0]=VB7_analyzeTrace(data,Wparent,initFun,N0,restarts,savename,rndState0)
% 
% analyzes a single trace using the VB7 method.
% data      : preprocessed single trajectory in VB7 format.
% Wparent   : VB7 opject with priorParameterOptions field 
%             (see VB7_priorParent)
% initFun   : function handle to a  constructor of initial guesses (see
%             e.g., VB7_initialGuess_KBregion).
% N0        : initial number of genuine states for search algoruthm
%             (defalut 20)
% restarts  : number of independent starts from scratch (default 20).
% savename  : Optional. If given, the local workspace is saved to this
%             file after each round of analysis (for debugging purposes).
% rndState0 : Optional. Sets the state of the pseudo random number
%             generator before analysis (for debugging).
% output:
% Wbest     : best converged VB7 structure
% NF        : [N Nc F round]. size, lower bound, and round when tested, for
%             all tested models.
% NFiMax    : [N Nc F round NF-line], collecting the NF values for the best
%             model of each size, and on which NF line it is.
% rndState0 : Optional. State of RandStream.getDefaultStream before
%             starting analysis. (for debugging). 
%% change-log
% M.L. 2010-10-20   : translated from VBE1 -> VB3
% M.L. 2010-11-08   : improved error handling
% m.l. 2011-02-04   : translated from VB3 -> VB7, removed CC objects for
%                     memory reasons.
% M.L. 2011-02-08   : improved time reporting
% M.L. 2011-02-14   : changed name of error message
% M.L. 2011-02-18   : change error handling so that just that round is
%                     aborted upon encountering an error
%% reset random numbers and save starting state
if(exist('rndState0','var'))
    rst = RandStream.create('mt19937ar','seed',sum(100*clock));
    RandStream.setGlobalStream(rst);
    global random_number_generator_reset_this_session;
    rst.State=rndState0;
    clear rst
end
if(~exist('random_number_generator_reset_this_session','var'))
    rst = RandStream.create('mt19937ar','seed',sum(100*clock));
    RandStream.setGlobalStream(rst);
    global random_number_generator_reset_this_session;
    rndState0=rst.State;
    clear rst
end
%% search parameters
% use default convergence parameters
relTolF=[];
tolPar=[];
maxIter=[]; 
outputLevel=0;
GtoS=false;
% wrapper for greedy search function
% [Wbest,N1,F1,C1]=VB7_greedySearch(W0,data,maxIter,relTolF,tolPar,outputLevel,GtoS)
searchFun=@(W,data)(VB7_greedySearch(W,data,maxIter,relTolF,tolPar,outputLevel,GtoS));
%% run VBEM iterations for various model sizes
NF=[];
W0=VB7_priorParent(Wparent,N0,1);  % initialize VB7 object with priors only
tic
tsave=toc;
tCompute0=toc;
Wbest.F=-Inf;
addWarning=false;
for cycle=1:restarts
    encounteredError=false;
    tCompute1=toc;
    Winit=initFun(W0); % initial guess
    disp(['initialized N = ' int2str(Winit.N) ', T = ' int2str(length(data.R2)) ...
        ', round ' int2str(cycle) ' of ' int2str(restarts)])
    % run actual iterations
    t0=toc;
    try
        [W,N1,F1]=searchFun(Winit,data);
    catch me
        me.getReport
        crashfile=sprintf('VBanalyzeTrace_error_%f.mat',now);
        save(crashfile)
        warning('VBanalyzeTrace:VBsearch_error',['VB7_greedySearch encountered an error. Saving state to ' crashfile])
        encounteredError=true;
        addWarning=true;
    end    
    if(addWarning)
        Wbest.warning='Analysis errors occurred during analysis of this trajectory';
    end        
    if(~encounteredError)
        if(W.F>Wbest.F) % then we have a new optimal
            Wbest=W;
        end
        
        % save all models from the greedy search
        for k=1:size(N1,1)
            NF(end+1,:)=[N1(k,:) F1(k) cycle];
        end
        NFiMax=[];
        for k=1:max(NF(:,1))
            ind=find(NF(:,1)==k);
            if(~isempty(ind))
                [~,iMax]=max(NF(ind,3));
                NcMax=NF(ind(iMax),2);
                iMax=ind(iMax);
                NFiMax(end+1,:)=[NF(iMax,:) iMax];
            end
        end
        tCompute=toc;
        [~,ib]=max(NFiMax(:,3));
        
        disp(['Round ' int2str(cycle) ' finished in ' num2str((tCompute-tCompute1)/60,3) ' min.' ...
            ' (average ' num2str((tCompute-tCompute0)/60/cycle,3) ' min./round). ' ...
            'Best fit so far : ' int2str(Wbest.N) ...
            ' states, from round ' int2str(NFiMax(ib,4)) '.'])
        if( (toc-tsave>120 || cycle==restarts) && exist('savename','var') && ~isempty(savename))
            ts0=toc;
            save(savename);%,'NF','Wparent','initFun','Wbest','NFiMax')
            ts1=toc;
            disp(['saved intermediate results: ' num2str(ts1-ts0,3) ' s.'])
            tsave=toc;
        end
    end    
end

