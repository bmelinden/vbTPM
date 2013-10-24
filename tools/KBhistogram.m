% [h,K,B,RMS,Krange,Brange,KBn]=KBhistogram(X,tSigma,fSample,KBrange,NKB,fig)
%
% h : plot handle
% K,B,RMS . gaussian filtered estimates
% KBrange =[Kmin Kmax Bmin Bmax] range for the histogram, which is plotted
% in figure fig (or axes(fig)), if the argument is given.
% NBK = [NK NB] is teh number of lattice points in each dimension of the
% histogram
% 
% M.L. 2010-10-17

%% change log
% M.L. 2011-02-02 : improved help text

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manualSegmentML.m, part of the vbTPM package
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
function [h,K,B,RMS,Krange,Brange,KBn]=KBhistogram(X,tSigma,fSample,KBrange,NKB,fig)

if(~exist('tSigma','var')|| isempty(tSigma));  tSigma=4;end
if(~exist('fSample','var')|| isempty(fSample));fSample=30;end
[RMS,K,B]=RMSKBgaussFilter(X,tSigma,fSample);

if(~exist('NKB','var')|| isempty(NKB));NKB=[100 100];end
if(length(NKB)==1); NKB=NKB*[1 1]; end
if(~exist('KBrange','var')|| isempty(KBrange));KBrange=[-1 1 0 max(B) ];end
   
Krange=linspace(KBrange(1),KBrange(2),NKB(1));
Brange=linspace(KBrange(3),KBrange(4),NKB(2));

[a,b]=hist3([K B],{Krange, Brange});
KBn=a';

if(exist('fig','var') && ~isempty(fig))
    try
        figure(fig);
    catch
        axes(fig);
    end
    h=pcolor(Krange,Brange,KBn/sum(sum(KBn)));
    axis(KBrange);
    ylabel('B');
    xlabel('K');
    set(h,'edgecol','non')
else
    h=[];
end

