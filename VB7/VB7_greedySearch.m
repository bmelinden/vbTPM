% [Wbest,N1,F1,C1]=VB7_greedySearch(W0,data,maxIter,relTolF,tolPar,outputLevel,GtoS)
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
% outputlevel: controls the amount of text output. See VB7iterator for
%              details. Default 0 (no output).
% GtoS      : attempt to optimize by converting states to spurious states
%             explicitely. Does not work well with current models (default:
%             false).
% 
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
% M.L. 2011-02-03   : translated to VB7, rewrote help text
% M.L. 2011-02-07   : added upper bound to the number of states
% M.L. 2011-02-11   : changed recursive function to while loop, and made
%                     storing C1 cell array optional, to try to reduce
%                     memory usage. Replaced recursive version of
%                     algorithm.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_greedySearch.m, model search in the vbTPM package
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
function [Wbest,N1,F1,C1]=VB7_greedySearch(W0,data,maxIter,relTolF,tolPar,outputLevel,GtoS)

% default parameters
% here we are mainly interested in convergence of F.
if(~exist('maxIter','var')  || isempty(maxIter));      maxIter=1000;  end
if(~exist('relTolF','var')  || isempty(relTolF));      relTolF=1e-9;  end
if(~exist('tolPar','var')   || isempty( tolPar));       tolPar=1e-3;  end
if(~exist('outputLevel','var')||isempty(outputLevel));outputLevel=1;  end
if(~exist('GtoS','var')   || isempty( GtoS));              GtoS=false;end

Nmax=50;
maxMoves=100;
saveC=(nargout==4);

%% prepare computations, and converge initial guess.
if(saveC);C1=cell(1,1);end
W0=VB7_VBEMiter(W0,data);
%W0=VB5_pSort(W0,W0.est.sAverage,W0.est.cAverage);
[W0,C1{1},F1]=VB7iterator(W0,data,maxIter,relTolF,tolPar,outputLevel);
N1=[W0.N W0.Nc];
W0=VB7_VBEMiter(W0,data);

if(outputLevel>0)
    disp(['Starting greedy state reduction from N = ' int2str([W0.N W0.Nc])])
end
% recursive search
%[N1,F1,C1]=greedySearch(C1,N1,F1,W0,data,maxIter,relTolF,tolPar,outputLevel,2*W0.N);

%% looping search
Wbest=W0;
oneMoreTry=true;
move=0;
while(oneMoreTry && move < maxMoves && Wbest.N <= Nmax)
    foundImprovement=false;
    move=move+1;
    W0=Wbest;
    
    [~,h]=sort(W0.est.sAverage); % order in increasing occupancy
    for k=1:length(h)
        if(W0.N>1)
            %% try to remove a looping state
            Wt=VB7_removeState(W0,h(k),[]);
            Wt=VB7_VBEMiter(Wt,data);
            [Wt,Ct]=VB7iterator(Wt,data,maxIter,relTolF,tolPar,outputLevel);
            if(Wt.F>= W0.F) % after a successful move, start over
                if(outputLevel>0)
                    disp('removing the state helped')
                    disp('-------------------------------------------------------------------------')
                end
                Wbest=Wt;
                if(saveC); C1{end+1}=Ct; end
                N1(end+1,:)=[Wt.N Wt.Nc];
                F1(end+1,1)=Wt.F;
                foundImprovement=true;
                break % do not make further attempts at this model size
            end
            %% if some significant occupation weight was removed, try adding some extra
            % transitions to all states, to compensate for potential lost
            % connectivity. ML: I think this should not be needed because the
            % most recent version of VB7_removeState already does this kind of
            % compensation. Still, it cannot hurt (except computational time)
            % to keep this move.
            if(W0.est.sAverage(h(k))>0.1/W0.T)
                Wt.E.wA=Wt.E.wA+2;
                [Wt,Ct]=VB7iterator(Wt,data,maxIter,relTolF,tolPar,outputLevel);
                % test for improvement
                if(Wt.F>W0.F) % after a successful move, start over
                    if(outputLevel>0)
                        disp('removing the state and adding extra counts helped')
                        disp('-------------------------------------------------------------------------')
                    end
                    Wbest=Wt;
                    if(saveC); C1{end+1}=Ct; end;
                    N1(end+1,:)=[Wt.N Wt.Nc];
                    F1(end+1,1)=Wt.F;
                    foundImprovement=true;
                    break % do not make further attempts at this model size
                end
            end
        end
        %% try to split state in two
        try
            Wt=VB7_split(W0,h(k),2);
        catch me
            me.getReport
            crashfile=sprintf('VB7_split_error_%f.mat',now);
            save(crashfile)
            error('VB7_greedySearch:VB7_split_error',['VB7_greedySearch encountered a VB7_split error. Saving state to ' crashfile])
        end
        Wt=VB7_VBEMiter(Wt,data);
        [Wt,Ct]=VB7iterator(Wt,data,maxIter,relTolF,tolPar,outputLevel);
        % test for improvement
        if(Wt.F>W0.F) % after a successful move, start over
            if(outputLevel>0)
                disp('splitting the state helped')
                disp('-------------------------------------------------------------------------')
            end
            Wbest=Wt;
            if(saveC); C1{end+1}=Ct; end;
            N1(end+1,:)=[Wt.N Wt.Nc];
            F1(end+1,1)=Wt.F;
            foundImprovement=true;
            break % do not make further attempts at this model size
        end
        if(GtoS)% try to convert genuine state to spurious state
            Wt=VB7_GSconversion(W0,h(k),[]);
            Wt=VB7_VBEMiter(Wt,data);
            [Wt,Ct]=VB7iterator(Wt,data,maxIter,relTolF,tolPar,outputLevel);
            % test for improvement
            if(Wt.F>W0.F) % after a successful move, start over
                if(outputLevel>0)
                    disp('conversion to spurious state helped')
                    disp('-------------------------------------------------------------------------')
                end
                Wbest=Wt;
                if(saveC); C1{end+1}=Ct; end;
                N1(end+1,:)=[Wt.N Wt.Nc];
                F1(end+1,1)=Wt.F;
                foundImprovement=true;
                break % do not make further attempts at this model size
            end
        end
    end
    
    % give up if W0 could not be improved upon
    if(~foundImprovement)
        oneMoreTry=false;
    end
end

%% sort results
if(saveC);C1=C1(end:-1:1);end
N1=N1(end:-1:1,:);
F1=F1(end:-1:1,:);

[~,b]=max(F1);
if(outputLevel>0)
    disp(['Ending greedy search with N = [' int2str([Wbest.N Wbest.Nc]) '] states.'])
end
end

