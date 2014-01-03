% S=VBviterbi_log0(lnQ,lnqst) 
% most likely trajectory by the Viterbi algorithm, using log of transition
% matrix lnQ and emission likelihood lnqst. This is the matlab version of
% VBviterbi_log.c, created for debugging purposes, but very much slower.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VBlogViterbi.m, viterbi algorithm, part of the vbTPM package
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
function S=VBviterbi_log0(lnQ,lnqst) 

[T,N]=size(lnqst);
% algorithm version without large lnpt array
lnPP=zeros(1,N);
lnP0=zeros(1,N);
lnP1=zeros(1,N);
MaxPrev=zeros(T,N,'uint16'); % state variable
lnP1=lnqst(1,:)-mean(lnqst(1,:)); % initial probability, not normalized

for tV=2:T
    lnP0=lnP1;
    lnP1=zeros(1,N);
    for jV=1:N      % current target state
        for kV=1:N  % previous state
            lnPP(kV)=lnP0(kV)     +lnQ(kV,jV)+lnqst(tV,jV); % probability of most likely path that ends with kV -> jV            
        end
        % probability of previous state before ending up at jV.
        [lnP1(jV),         MaxPrev(tV,jV)]=max(lnPP);         
    end
    lnP1=lnP1-mean(lnP1); % rescale to avoid numerical under- or overflow.
end
S=zeros(T,1,'uint16');
[~,S(T)]=max(lnP1);  
for tV=T-1:-1:1
    S(tV)=MaxPrev(tV+1,S(tV+1));
end
end

