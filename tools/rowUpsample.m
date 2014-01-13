% upsamples the rows of the matrix v by a factor ds, and then truncates it
% to Tmax rows. The first row is not upsampled. This is meant to make the
% results of a downsampled VB3 analysis comparable to data that has not
% been downsampled. 
% M.L.2010-12-03

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rowUpsample.m, part of the vbTPM package
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
function V=rowUpsample(v,ds,Tmax)
if(~exist('Tmax','var')|| isempty(Tmax))
    Tmax=1+ds*(size(v,1)-1);
end

V=zeros(Tmax,size(v,2));

V(1,:)=v(1,:);
for k=2:size(v,1)
   r1=1+(k-2)*ds;
   for m=1:ds
       V(r1+m,:)=v(k,:);
   end
end
V=V(1:Tmax,:);

