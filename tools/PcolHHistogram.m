%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PcolHHistogram.m, part of the vbTPM package
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
function [h,x,y,z]=PcolHHistogram(K,B,KBrange,NKB,fig)


if(~exist('NKB','var')|| isempty(NKB));NKB=[100 100];end
if(length(NKB)==1); NKB=NKB*[1 1]; end
if(~exist('KBrange','var')|| isempty(KBrange));KBrange=[-1 1 0 max(B) ];end
   


if(exist('fig','var'))
    try
        figure(fig);
    catch
        axes(fig);
    end
end

Krange=linspace(KBrange(1),KBrange(2),NKB(1));
Brange=linspace(KBrange(3),KBrange(4),NKB(2));

[a,b]=hist3([K B],{Krange, Brange});
x=b{1};
y=b{2};
z=a'/sum(sum(a));
h=pcolor(x,y,z);
axis(KBrange);
set(h,'edgecol','non')
