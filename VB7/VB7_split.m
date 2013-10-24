%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_VBEMiter_nomex.m, VBEM iteration without mex files, in the vbTPM package
% =========================================================================
% 
% Copyright (C) 2013 Martin Lindén
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
function W1=VB7_split(W0,ss,TsplitFactor)
% W=VB7_split(W0,ss,TsplitFactor)
% splits genuine state ss into two distinct genuine states, with mean dwell
% times differing by a factor TsplitFactor (each new state gets the
% original ss a factor sqrt(TsplitFactor) bigger or smaller dwell time a
% factor sqrt(TsplitFactor).

%% change-log
% M.L. 2011-01-06   : started
% M.L. 2011-02-03   : translated from VB5 to VB7, new rescaling formulae
% M.L. 2011-02-14   : throw error on split of empty state. Debug thing,
%                     because that problem seems to depend in initiakl
%                     guess, not data set.
%% actual code
N=W0.N+1;
Nc=W0.Nc;
W1=VB7_priorParent(W0,N,Nc); % create extended model

%% expand emission parameters
sk=[1:W0.N ss];
W1.M.wPi=W0.M.wPi(sk);
W1.M.n  =  W0.M.n(sk);
W1.M.c  =  W0.M.c(sk);
W1.M.mu = W0.M.mu(sk);
W1.M.v  =  W0.M.v(sk);
% transitions to spurious states
W1.Mc=W0.Mc;
W1.Mc.wA=W1.Mc.wA(sk,:);
%% expand and rescale transition probabilities
if(W0.N>1)
    % rescale dwell times of the new states
    wA=W0.M.wA-W0.PM.wA; % observed transitions before conversion
    if(sum(wA(ss,:))==0 || sum(wA(:,ss))==0)
        error('VB7_split:unpopulated',['VB7_split: split of unpopulated state causes division by zero'])
        %warning('VB7_split: split of unpopulated state. Including pseudocounts in split')
        wA=wA+W0.PM.wA;
    end
    % compute new staying and interconversion probabilities
    p11=wA(ss,ss)/sum(wA(ss,:)); % old staying probability of state ss
    f=sqrt(TsplitFactor);
    pMin=0.01;
    pL=max(pMin,1-f*(1-p11)); % f-fold decrease in mean dwell time (with lower bound)
    pH=1-1/f^2*(1-pL);        % f^2-fold increase in mean dwell time compared to pL
    % old formulae
    % pH=p11*f/(1+p11*(f-1));   % new staying probabilities, with rescaled
    % pL=p11/f/(1+p11*(1/f-1)); % dwell times
    % pHL=(1-p11)/(W0.N-1);     % interconversion rate among new states
    
    % compute new conversion probaility vectors so that exit from each
    % state has 50% chance to interconvert between the newly split states
    sRow=wA(ss,:); % outbound transitions
    sRow(ss)=0;
    sRow=sRow/sum(sRow);           % normalize
    sRowH=[sRow 0]*(1-pH)/2;
    sRowL=[sRow 0]*(1-pL)/2;
    sRowH([ss end])=[   pH    (1-pH)/2];
    sRowL([ss end])=[(1-pL)/2    pL   ];
    
    W1.M.wA=wA(sk,sk);
    W1.M.wA(ss,:)=sRowH;
    W1.M.wA(end,:)=sRowL;
    W1.M.wA=W1.M.wA+W1.PM.wA;
else
    % with no other genuine states, there is no dwell time to start
    % from, so we have to make one up. I choose the time scale
    % to 10*(noise correlation time), or 10, whichever is bigger
    n=10*max(1,-1/log(W0.est.sKaverage)/W0.priorParameterOptions.downSample);
    W1.M.wA=10000*[1-1/n   1/n; 1/n    1-1/n];
    
    % we also need to differentiate the emission parameters
    % slightly to introduce some asymmetry as a basis for splitting
    K0=W0.M.mu;
    K1=K0*[0.9 1.1];
    if(K1(2)>0.98)
        K1=K1/K1(2)*0.98;
    end
    W1.M.mu=K1;
    W1.M.v=0.5*[1 1]*W0.M.v;
    
    B0=W0.M.n./W0.M.c;
    a=-B0/(1-K0);      % from the observation that B ~ a*K+b for genuine states
    B1=B0+a*(K1-K0);
    W1.M.n=0.5*W0.M.n*[1 1];
    W1.M.c=W1.M.n./B1;
    
end
