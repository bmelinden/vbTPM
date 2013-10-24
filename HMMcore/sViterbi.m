% [s,tWarn]=sViterbi(Q,qst)
% most likely trajectory by the Viterbi algorithm
% input Q, qst are the structures W.Q, W.qst in a converged VBE1 object
% (probably works for other algorithms too), i.e.,
% Q is proportional to a transition matrix, and
% qst is proportional to a local emission likelihood
% (see VBE1_VBEiter4.m for code that computes Q, qst that this algorithm is
% intended for).
% tWarn is a list of time points where normalization roblems occur, because
% of numerical underflow.
%
% This is the matlab version of VBviterbi.c, and is kept for debugging and
% benchmarking reasons (the c version is MUCH faster).

% M.L. 2011-12-20 : created a C version, and checked that they agree
%                   numerically.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sViterbi.m, viterbi algorithm, part of the vbTPM package
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
function [s,tWarn]=sViterbi(Q,qst) 

[T,N]=size(qst);
pt=zeros(T,N);
maxPrevious=zeros(T,N);
pt(1,:)=qst(1,:)/sum(qst(1,:)); % initial probability
pp=zeros(1,N);
tWarn=[];
for tV=2:T
    for jV=1:N      % current target state
        for kV=1:N  % previous state
            pp(kV)=pt(tV-1,kV)*Q(kV,jV)*qst(tV,jV); % probability of most likely path that ends with kV -> jV
        end
        if(sum(pp==0)>0)            
            warning(['non-normalizeable: pp(' int2str(tV) ') = [' num2str(pp,3) '].'])
            tWarn(end+1)=tV;
        end
        % probability of previous state before ending up at jV.
        [pt(tV,jV),maxPrevious(tV,jV)]=max(pp); 
    end
    pt(tV,:)=pt(tV,:)/sum(pt(tV,:)); % normlize, since we only need to keep 
                                     % track of relative probabilities.
end
s=zeros(T,1);
[~,s(T)]=max(pt(T,:));  
for tV=T-1:-1:1
    s(tV)=maxPrevious(tV+1,s(tV+1));
end
if(~isempty(tWarn))
    tWarn=union(tWarn(1),tWarn);
end
end
