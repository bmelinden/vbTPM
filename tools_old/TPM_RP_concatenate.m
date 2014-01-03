% x=TPM_RP_concatenate(lacX,lacY)
%
% function to transform raw data in the form 
%  lacX: {[M1 x120x B1 double]  ... [M4 x120x B4 double] }
%  lacY: {[M1 x120x B1 double]  ... [M4 x120x B4 double] }
% to 
%  x{i}{j} = [Mi*120 x 2 double], where i is the set index, and j is the bead
%
% This is the raw data format of xy-coordinate data saved by the RP group
% analysis software (as of now)
%
% M.L. 2010-06-04

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TPM_RP_concatenate.m, converge VBEM iterations in the vbTPM package
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
function x=TPM_RP_concatenate(lacX,lacY)

x=cell(size(lacX));
for i=1:length(lacX) % loop over bead sets
    if(~isempty(lacX{i}))
        x{i}=cell(1,size(lacX{i},3));
        for j=1:length(x{i}); % loop over beads in set k
            x{i}{j}=zeros(120*size(lacX{i},1),2);
            for k=1:size(lacX{i},1)
                x{i}{j}(120*(k-1)+(1:120),1:2)=[lacX{i}(k,:,j)' lacY{i}(k,:,j)'];
            end
        end
    end
end
