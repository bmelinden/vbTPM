function [Wbest,N1,F1,C1]=VB7_greedySearch_recursive(W0,data,maxIter,relTolF,tolPar,outputLevel,GtoS)
%% [Wbest,N1,F1,C1]=VB7_greedySearch(W0,data,maxIter,relTolF,tolPar,outputLevel,GtoS)
%
% Starting from model W0 (not necessarily converged), performs a greedy
% search for models with fever states and some dirt states. Working in
% order of increasing state occupancy of genuine states, the algorithm
% tries to 1) remove states, 2) split states, and 3) convert 'genuine'
% states to 'spurious' states (if GtoS=true, default false)  
% 
% The search is a greedy recursive search, i.e., as soon as a change
% increases the lower bound F of the fit, the algorithm accepts that change
% and starts the next round of refinements.
% It stops when no more states can be removed or converted.
% 
% INPUT:
% W0        : a VB7 model to start search from
% data      : the data set, after preprocessing with VB7_preprocess
% 
% maxIter, relTolF, tolPar: optional convergence parameters. See help
%                           VB7iterator for details and default values.
% OUTPUT: 
% Wbest     : the optimal model 
% N1,F1     : colum matrices with model dimensions and lower bounds on all
%             models encountered during the greedy search.
% C1        : cell vector of convergence structures from VB7iterator of all
%             models encountered during the search.

%% change-log
% M.L. 2010-07-02   : added prior generator function as to input list
% M.L. 2010-09-30   : spawned off VBE1_greedyReduce, that does a greedy
%                     search of state reduction, added help text
% M.L. 2010-10-17   : made sure that recursion stops at N=1 (bugfix), and
%                     that only one N=1 model is tried (saves time).
% M.L. 2010-10-21   : translated to VB3 (and future iteration schemes)
% M.L. 2010-11-03   : adapted to new version of VBiterator, saved old
%                     version in other file. 
% M.L. 2010-11-26   : tweaked the state reduction scheme: now it only tried
%                     to add extra transition pseudocount if the removed
%                     state was occupied at all; this both speeds up
%                     computations and improves the search.
% M.L. 2011-01-06   : translated to VB5
% M.L. 2011-01-07   : added and tested a state-splitting move (looping
%                     states only).
% M.L. 2011-02-03    : translated to VB7, rewrote help text
% M.L. 2011-02-07    : added upper bound to the number of states

%% default parameters
% here we are mainly interested in convergence of F.
if(~exist('maxIter','var')  || isempty(maxIter));      maxIter=1000;  end
if(~exist('relTolF','var')  || isempty(relTolF));      relTolF=1e-9;  end
if(~exist('tolPar','var')   || isempty( tolPar));       tolPar=1e-3;  end
if(~exist('outputLevel','var')||isempty(outputLevel));outputLevel=1;  end
if(~exist('GtoS','var')   || isempty( GtoS));                GtoS=false;end
%% prepare computations, and converge initial guess.
C1=cell(1,1);
W0=VB7_VBEMiter(W0,data);
%W0=VB5_pSort(W0,W0.est.sAverage,W0.est.cAverage);
[W0,C1{1},F1]=VB7iterator(W0,data,maxIter,relTolF,tolPar,outputLevel);
N1=[W0.N W0.Nc];
W0=VB7_VBEMiter(W0,data);

if(outputLevel>0)
    disp(['Starting greedy state reduction from N = ' int2str([W0.N W0.Nc])])
end
% recursive search
[N1,F1,C1]=greedySearch(C1,N1,F1,W0,data,maxIter,relTolF,tolPar,outputLevel,2*W0.N);
% sort results
C1=C1(end:-1:1);
N1=N1(end:-1:1,:);
F1=F1(end:-1:1,:);

[~,ind]=max(F1);
Wbest=C1{ind}.W1;

[~,b]=max(F1);
if(outputLevel>0)
    disp(['Ending greedy search with N = [' int2str(N1(b,:)) '] states.'])
end
%% recursion function
function [N1,F1,C1]=greedySearch(C1,N1,F1,W0,data,maxIter,relTolF,tolPar,outputLevel,Nmax)
if(W0.N<=1 || W0.N >= Nmax) % stop
    return
else
    [~,h]=sort(W0.est.sAverage); % order in increasing occupance
    
    foundImprovement=false;
    
    FbestThisLevel=-Inf;
    NestThisLevel=[W0.N W0.Nc];
    CbestThisLevel=W0;
    for k=1:length(h)
        if(outputLevel>0)
            disp(['Greedy attempt to modify state ' int2str(h(k)) ' ( ' int2str(k) ' of ' int2str(length(h)) ').' ])
        end
        %% try to remove a looping state
        Wt=VB7_removeState(W0,h(k),[]);        
        Wt=VB7_VBEMiter(Wt,data);
        [Wt,Ct]=VB7iterator(Wt,data,maxIter,relTolF,tolPar,outputLevel);                            
        %% test for improvement
        if(Wt.F>FbestThisLevel);
            FbestThisLevel=Wt.F;
            NestThisLevel=[Wt.N Wt.Nc];
            CbestThisLevel=Ct;
        end
        if(Wt.F>= W0.F) % after a successful move, start over
            if(outputLevel>0)
                disp('removing a state helped')
                disp('-------------------------------------------------------------------------')
            end
            C1{end+1}=Ct;
            N1(end+1,:)=[Wt.N Wt.Nc];
            F1(end+1,1)=Wt.F;
            [N1,F1,C1]=greedySearch(C1,N1,F1,Wt,data,maxIter,relTolF,tolPar,outputLevel,Nmax);
            foundImprovement=true;
            break % do not make further attempts at this model size
        end
        %% if some significant occupation weight was removed, try adding some extra
        % transitions to all states, to compensate for potential lost
        % connectivity. This makes no sense: if no significant weight was
        % removed, and that did not improve F, then something is wrong.
        if(W0.est.sAverage(h(k))>0.1/W0.T)
            Wtt=Wt;
            Wtt.E.wA=Wtt.E.wA+2;
            [Wtt,Ctt]=VB7iterator(Wtt,data,maxIter,relTolF,tolPar,outputLevel);
            if(Wtt.F>Wt.F)
                Wt=Wtt;
                Ct=Ctt;
            end
            %% test for improvement
            if(Wt.F>FbestThisLevel);
                FbestThisLevel=Wt.F;
                NestThisLevel=[Wt.N Wt.Nc];
                CbestThisLevel=Ct;
            end
            if(Wt.F>= W0.F) % after a successful move, start over
                if(outputLevel>0)
                    disp('removing a state and adding extra counts helped')
                    disp('-------------------------------------------------------------------------')
                end
                C1{end+1}=Ct;
                N1(end+1,:)=[Wt.N Wt.Nc];
                F1(end+1,1)=Wt.F;
                [N1,F1,C1]=greedySearch(C1,N1,F1,Wt,data,maxIter,relTolF,tolPar,outputLevel,Nmax);
                foundImprovement=true;
                break % do not make further attempts at this model size
            end
        end
        %% try to split state in two (TBA)        
        Wt=VB7_split(W0,h(k),2);
        Wt=VB7_VBEMiter(Wt,data);
        [Wt,Ct]=VB7iterator(Wt,data,maxIter,relTolF,tolPar,outputLevel);                    
        %% test for improvement
        if(Wt.F>FbestThisLevel);
            FbestThisLevel=Wt.F;
            NestThisLevel=[Wt.N Wt.Nc];
            CbestThisLevel=Ct;
        end        
        if(Wt.F>= W0.F) % after a successful move, start over
            if(outputLevel>0)
                disp('splitting a state helped!')
                disp('-------------------------------------------------------------------------')
            end
            C1{end+1}=Ct;
            N1(end+1,:)=[Wt.N Wt.Nc];
            F1(end+1,1)=Wt.F;
            [N1,F1,C1]=greedySearch(C1,N1,F1,Wt,data,maxIter,relTolF,tolPar,outputLevel,Nmax);
            foundImprovement=true;
            break % do not make further attempts at this model size
        end
        
        if(GtoS)% try to convert a state to a dirt state
            Wt=VB7_GSconversion(W0,h(k),[]);
            Wt=VB7_VBEMiter(Wt,data);
            [Wt,Ct]=VB7iterator(Wt,data,maxIter,relTolF,tolPar,outputLevel);
            %% test for improvement
            if(Wt.F>FbestThisLevel);
                FbestThisLevel=Wt.F;
                NestThisLevel=[Wt.N Wt.Nc];
                CbestThisLevel=Ct;
            end
            if(Wt.F>= W0.F) % after a successful move, start over
                if(outputLevel>0)
                    disp('dirt state conversion helped')
                    disp('-------------------------------------------------------------------------')
                end
                C1{end+1}=Ct;
                N1(end+1,:)=[Wt.N Wt.Nc];
                F1(end+1,1)=Wt.F;
                [N1,F1,C1]=greedySearch(C1,N1,F1,Wt,data,maxIter,relTolF,tolPar,outputLevel,Nmax);
                foundImprovement=true;
                break % do not make further attempts at this model size
            end
        end
    end
    % no reduction or conversion from W0 was an improvement. Still, record
    % the least bad attempt.
    if(~foundImprovement)
        if(outputLevel>0)
            disp(['Greedy search could not improve N = ' int2str([W0.N W0.Nc]) ' model.'])
        end
        C1{end+1}=CbestThisLevel;
        N1(end+1,:)=NestThisLevel;
        F1(end+1,1)=FbestThisLevel;
    end
end
end
end

