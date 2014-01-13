% [STATE,DWELL,Dstart,Dend,wA]=getDwellTRJ(s,spurious)
% converts a state trajectory s into a list of sojourn states STATE, and a
% list of dwell times DWELL.
% dwell m occurred from index Dstart(m) to index Dend(m)
% wA is a matrix of transition counts, comparable to the one used in the
% VBEM algorithm, wA(k,j) = number of k -> j transitions in s.
%
% spurious : (optional) list of spurious states in the data. Dwells with spurious
%            states are divided between the surrounding non-spurious
%            states. Default: empty.

% M.L. 2010-06-17   : added Dstart,Dend.
% M.L. 2011-01-21   : added wA
% M.L. 2012-03-20   : added ability to remove spurious states

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getDwellTRJ.m, part of the vbTPM package
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

function [STATE,DWELL,Dstart,Dend,wA]=getDwellTRJ(s,spurious)

% remove the ignore states
if(exist('ignore','var') && ~isempty(ignore))
    t1=find(~ismember(s,ignore),1,'first');
    s1=s(t1); % first non-spurious state
    s(1:t1-1)=s1; % first spurious dwell overwritten by next non-spurious state
    
    for t=t1+1:length(s)
        if(ismember(s(t),ignore)) % then 
            if(~ismember(s(t-1),ignore)) % then a spurious dwell just started
                t1=t;      % beginning of spurious dwell
                s1=s(t-1); % remember last genuine state
            end
            % if both s(t) and s(t-1) are spurious, then we should keep
            % looking for the 
        else
            if(ismember(s(t-1),ignore)) % then a spurious dwell just ended
               t2=t-1;
               s2=s(t);
               tmid=floor(0.5*(t1+t2));
               s(t1:tmid)=s1;   % divide up the dwell between surrounding 
               s(tmid+1:t2)=s2; % genuine states
            end
            % if both s(t) and s(t-1) are genuine states, then everything
            % is find
        end
    end
    clear s1 s2 t1 t2 tmid
end

sNow=s(1);
dwellSoFar=1;
STATE=[];
DWELL=[];
Dstart=1;
Dend=[];
N=max(s);
wA=zeros(N,N);
for k=2:length(s)
    wA(s(k-1),s(k))=wA(s(k-1),s(k))+1;
    if(sNow==s(k)) % one more step in the same state
        dwellSoFar=dwellSoFar+1;
    else % transition happened
         % store old dwell time
         STATE(end+1)=sNow;
         Dend(end+1)=k-1;
         Dstart(end+1)=k;
         % start new dwell time count
         sNow=s(k);
         dwellSoFar=1;        
    end
end
Dend(end+1)=length(s);
STATE(end+1)=s(end);
DWELL=1-Dstart+Dend;



end
