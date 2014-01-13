% [X,T,dt,Xraw]=getBeadXY(source) 
% recover xy-trace from the data files in directory source
% the raw positions are filtered with a 0.1 Hz first order Butterworth
% filter to remove drift.
% M.L. 2010-01-09   : added a loop to cut out the first non-NaN interval
%                     of each time series.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getBeadXY.m, part of the vbTPM package
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

function [X,T,dt,Xraw]=getBeadXY(source) 

% get all POS file names
if(~exist('source','var'))
   source=pwd; 
end

%% read all data files and concatenate data
a=dir([source '/*POS.mat']);
NF=length(a);                    % number of files 
foo=load([source '/' a(1).name]);NB=foo.nroi; % number of beads
x=zeros(NF*120,NB);y=x;
disp(['Found ' int2str(NF) ' files.']);
for k=1:NF
    b=load([source '/' a(k).name]);
    % figure out which position this file is, i.e., find n in 
    % filename_n_POS.mat, where n is an integer of undetermined length 
    foo=regexp(a(k).name,'_','split');
    n=str2double(foo{end-1});
    
    x(1+(n-1)*120:n*120,1:NB)=b.pm_xy1(:,:,1);
    y(1+(n-1)*120:n*120,1:NB)=b.pm_xy1(:,:,2);
end

%% next design a high-pass Butterworth filter as described by Han et al
Fsample=30;  % sampling frequency (~30 Hz), rad/s
Fnyquist=Fsample/2;

Fcutoff=0.1;  % cut-off at 0.1 Hz [rad/s]
Wn=Fcutoff/Fnyquist;
[Z,P,K]=butter(1,Wn,'high');
[sos,g]=zp2sos(Z,P,K);
Hd=dfilt.df2tsos(sos,g);   % Create a dfilt object

dt=1/Fsample;
%% take out and filter the first interval of non-NaN of each trace
X=cell(1,NB);
Xraw=cell(1,NB);
T=zeros(1,NB);
offset=0;
for k=1:NB
    n2x=find(isnan(x(:,k)),1)-1;
    n2y=find(isnan(x(:,k)),1)-1;
    n2=min(n2x,n2y);
    if(length(n2)==1)
        X{k-offset}=[filter(Hd,x(1:n2,k)) filter(Hd,y(1:n2,k))];
        Xraw{k-offset}=[x(1:n2,k) y(1:n2,k)];
        T(k-offset)=n2;
    else
        offset=offset+1;
    end
end
X=X(1:NB-offset);
Xraw=Xraw(1:NB-offset);
T=T(1:NB-offset);
