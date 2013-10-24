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
function W1=VB7_removeState(W0,rmG,rmS)
%% W1=VB7_removeState(W0,rmG,rmS)
%
% function to remove single states (genuine or spurious)
% 
% W0    : parent VB7 structure
% rmG,rmS: indices of clean/dirt states that are to be changed. For
% now, only one state can be changed at a time.

%% change-log
% M.L. 2010-12-27   : started on VB4 converter
% M.L. 2011-01-05   : started on VB5 converter
% M.L. 2011-02-03   : translated to VB7 converter, and then to VB7 remover
%% action
if(~exist('rmS'))
    rmS=[];
end

if(length(rmG)==1 && isempty(rmS) && W0.N>1 ) % remove a genuine state
    clear rmS;
    N=W0.N-1;
    Nc=W0.Nc;
    W1=VB7_priorParent(W0,N,Nc); % new prior parameters
    if(~isfield(W1,'PM'))
        warning('VB7_convertStates: incoming state did not come with correct prior parameters. Using defaults.')
        W1=VB7_priorParent(W1,N,Nc);
    end   
    % transfer spurious state
    W1.Mc.wPi= W0.Mc.wPi;
    W1.Mc.n  = W0.Mc.n;
    W1.Mc.c  = W0.Mc.c;
    W1.Mc.mu = W0.Mc.mu;
    W1.Mc.v  = W0.Mc.v;

    
    % remove clean state    
    sk=[1:rmG-1 rmG+1:W0.N];
    W1.M.wPi=W0.M.wPi(sk);
    W1.M.n  =  W0.M.n(sk);
    W1.M.c  =  W0.M.c(sk);
    W1.M.mu = W0.M.mu(sk);
    W1.M.v  =  W0.M.v(sk);
    
    % transfer observed transitions
    wA =W0.M.wA -W0.PM.wA;
    wAc=W0.Mc.wA-W0.PMc.wA;
    
    W1.Mc.wR=W0.Mc.wR;            % no change in S->S transitions
    W1.Mc.wA=wAc(sk,:)+W1.PMc.wA; % remove old genuine state         
    % try to compensate for transitions via the removed state
    toR=wA(sk,rmG);
    frR=wA(rmG,sk);    
    W1.M.wA=wA(sk,sk)+W1.PM.wA+toR*frR;
elseif(isempty(rmG)   && length(rmS)==1 && rmS~=1 ) % convert a dirt state to clean
    clear rmG;
    N=W0.N;
    Nc=W0.Nc-1;
    W1=VB7_priorParent(W0,N,Nc); % new prior parameters
    if(~isfield(W1,'PM'))
        warning('VB7_convertStates: incoming state did not come with correct prior parameters. Using defaults.')
        W1=VB7_priorParent(W1,N,Nc);
    end    
    % transfer genuine states
	W1.M.wPi= W0.M.wPi;
    W1.M.n  = W0.M.n;
    W1.M.c  = W0.M.c;
    W1.M.mu = W0.M.mu;
    W1.M.v  = W0.M.v;

    % remove spurious state
    sk=[1:rmS-1 rmS+1:W0.Nc]; % dirt states to keep
    W1.Mc.wPi=W0.Mc.wPi(sk);
    W1.Mc.n = W0.Mc.n(sk);
    W1.Mc.c = W0.Mc.c(sk);
    W1.Mc.mu=W0.Mc.mu(sk);
    W1.Mc.v = W0.Mc.v(sk);
    
    % move transition counts around
    W1.M.wA=W0.M.wA;        % no change in G->G transitions
    wAc=W0.Mc.wA-W0.PMc.wA; % 
    wR =W0.Mc.wR-W0.PMc.wR;   
    
    % try to compensate for transitions via the removed state
    toR=wR(sk,rmS);
    frR=wR(rmS,sk);    
    W1.Mc.wR=wR(sk,sk)+W1.PMc.wR+toR*frR;

    toR=wAc(:,rmS);
    frR=wR(rmS,sk);
    W1.Mc.wA=wAc(:,sk)+W1.PMc.wA+toR*frR;    
    
else
    disp(['VB7_convertStates : rmG  = [ ' int2str(rmG(1:end)) ' ]'])
    disp(['VB7_convertStates : rmS = [ ' int2str(rmS(1:end)) ' ]'])
    error('VB7_convertStates : only one state can be converted at a time, and dirt state 1 canot be converted.')
end
