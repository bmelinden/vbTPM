function [W,C,F]=VB7iterator_oldsym(W,x,maxIter,relTolF,tolPar,outputLevel)
%%
% [W,C,F]=VB7iterator(W,x,maxIter,relTolF,tolPar,outputLevel)
% run VBE iterations using VB7_VBEMiter_oldsym(W,x), starting from VB7 structure W
% and data x, using convergence criteria explained below. C is an output
% structure with information about iteration parameters and convergence
% status, and F is the lower bound.
%
% maxIter   : maximum number of VBE iterations to run.
% relTolF   : convergence criterion for relative change in likelihood bound.
% tolPar    : convergence criterion for M-step parameters.
%             iterate until 
%             |dF/F| <= relTolF and 
%             max|dPar(i)/Par(i)| <= tolPar, 
%             i = X.p, X=M,Mc,
%             p=n,c,v,mu,wA,wR
% Default values are maxIter=1000, all relTolF=1e-6, tolPar=1e-3
% outputLevel: 0: no output, 1: display convergence, not progress (default), 
%              2: display convergence measures for each iteration
%
% note: this file is for debugging purposes only, and should be deleted
% eventually.

%% change log
% M.L. 2010-03-25 : first version
% M.L. 2010-07-02 : decreased minimum number of iterations to 2, to save
%                    time.
% M.L. 2010-10-21 : new file VBiterator, accepting a function handle for
%                    the actual iterations (thus forward compatible).
% M.L. 2010-10-27 : Handling errors from VB4_VBEMiter.
% M.L. 2010-11-03 : Adding more convergence criteria,
%                    and iterate until everything is converged. Saved old
%                    version to other file.
% M.L. 2010-11-08 : saving crash file in case of maximum number of
%                    iterations as well, and decreasing maxIter to 1e3.
% M.L. 2010-12-27 : translated to VB4, and got rid of the VBiterFun option.
% M.L. 2011-01-06 : translated to VB5
% M.L. 2011-02-02 : translated to VB7, turned off saving state when
%                   reaching maximum iterations
% M.L. 2011-03-01 : version that uses the VBEMite_oldsym iterations, used
%                   to debug dusing an update of the notation.

%% set convergence parameters
if(~exist('maxIter','var')  || isempty(maxIter));maxIter=1000;end
if(~exist('relTolF','var')  || isempty(relTolF));relTolF=1e-6;end
if(~exist('tolPar','var')   || isempty(tolPar)); tolPar =1e-3;end
C.maxIter =maxIter;
C.relTolF =relTolF;
C.tolPar  =tolPar;
C.iter    =0;
C.converged=false;
C.W0=W;
if(isfield(W,'minimalStorage')) % save memory by saving a small object...
    C.W0=W.minimalStorage(); 
end
%% run the VBE iterations until convergence
if(~exist('outputLevel','var') || isempty(outputLevel))
    outputLevel=1;
end
displayProgress=false;
displayExit=true;
if(outputLevel==0) % adjust output settings
    displayExit=false;
elseif(outputLevel==2)
    displayProgress=true;
elseif(outputLevel~=1)
    warning('VB7iterator: Undefined outputLevel. Using default value')
end
runMore=true;
C.exitStatus='';
Wm5=[];Wm4=[];Wm3=[];Wm2=[];Wm1=[];
while(runMore)
    % keep a short history in case something goes wrong...
    Wm5=Wm4;Wm4=Wm3;Wm3=Wm2;Wm2=Wm1;Wm1=W;
    Wold=W;
    try
        W=VB7_VBEMiter_oldsym(W,x);
    catch me
        me.getReport
        crashfile=sprintf('VB7iterator_VBiter_error_%f.mat',now);
        save(crashfile)
        runMore=false;
        C.exitStatus=[C.exitStatus 'VB7_VBEMiter_oldsym returned error'];
        error(['VB7_VBEMiter_oldsym returned error. Saved state to ' crashfile])
    end
    C.iter=C.iter+1;

    if(isfield(Wold,'F')) % then we can check convergence
        %% converge criterion in terms of relative changes in F and parameter values
        if(~isfinite(W.F)) % check for problem
            crashfile=sprintf('VB7iterator_VBiter_N aNInf_%f.mat',now);
            save(crashfile);
            runMore=false;
            error(['VB7iterator found W.F=NaN or Inf. Saving state to ' crashfile])
            C.exitStatus=[C.exitStatus 'W.F=NaN or Inf. '];
            pause(0.5)
        end
        C.dFrel=(W.F-Wold.F)/abs(W.F);
        
        fM=fields(W.M);
        C.dPrel=-Inf;
        C.limitPar='';
        for k=1:length(fM)
            Pnew=W.M.(fM{k})(1:end);
            Pold=Wold.M.(fM{k})(1:end);
            ind=find(Pnew~=0);
            dPrel=max(abs((Pnew(ind)-Pold(ind))./Pnew(ind)));
            if(dPrel>C.dPrel)
                C.limitPar=['M.' fM{k}];
                C.dPrel=dPrel;
            end
        end
        fM=fields(W.Mc);
        for k=1:length(fM)
            Pnew=W.Mc.(fM{k})(1:end);
            Pold=Wold.Mc.(fM{k})(1:end);
            ind=find(Pnew~=0);
            dPrel=max(abs((Pnew(ind)-Pold(ind))./Pnew(ind)));
            if(dPrel>C.dPrel)
                C.limitPar=['Mc.' fM{k}];
                C.dPrel=dPrel;
            end
        end
        
        Fconverged=(abs(C.dFrel)<C.relTolF);
        Pconverged=C.dPrel<C.tolPar;
        if(Fconverged && Pconverged)
            C.exitStatus=['Converged normally after ' int2str(C.iter) ' iterations, relTolF = ' num2str(relTolF) ', tolPar = ' num2str(tolPar) ', (' C.limitPar ' limiting).'];
            C.converged=true;
            runMore=false;
        end
        if(C.iter>=maxIter) % check for max iterations
            runMore=false;
            C.exitStatus=[C.exitStatus 'Maximum number of iterations reached. '];
%            crashfile=sprintf('VB7iterator_maxiter_%f.mat',now);
%            disp('Maximum number of iterations reached. saving state to')
%            disp(crashfile)
%            save(crashfile)
        end
        if(displayProgress) % display convergence progress
            if(mod(C.iter,10)==1)
                displayHeader();
            end
            displayConvergence();
        end
        if(C.iter<2) % run at least 2 iterations
            runMore=true;
            C.exitStatus='';
        end
    end
end
C.W1=W.minimalStorage(); % save final state parameters
C.F=W.F;
F=W.F;
C.exitStatus=C.exitStatus(1:end-1);
if(displayExit) % display exit message
    displayHeader()
    displayConvergence()
    fprintf('%s \n',C.exitStatus)
end
%% auxiliary functions
    function displayConvergence()
        fprintf('%02d %02d % 5d % 0.2e ',[W.N W.Nc C.iter C.dFrel]);
        fprintf('%0.2e  %s ',C.dPrel,C.limitPar)
        fprintf(' p(s,c) :')
        fprintf(' %0.2e ',[W.est.sAverage(W.est.sAverage>1/W.T) W.est.cAverage([false W.est.cAverage(2:end)>1/W.T]) ])
        %fprintf('% 0.2e ',[W.F-Wold.F W.est.Fterms-Wold.est.Fterms]);
        fprintf('\n')
    end
    function displayHeader()
        disp([' N NH     it   dF/|F|   d<p>/p   p '])
    end
end