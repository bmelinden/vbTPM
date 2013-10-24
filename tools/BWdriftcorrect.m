% x=BWdriftcorrect(x0,Fcut,Fsample)
%
% filter the columns of x0 with a hign pass 1st order Butterworth filter
% with cut-off frequency Fcut (default 0.05 Hz).
% Fsample = smapling frequency of x0 (default 30 Hz).


%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BWdriftcorrect.m, part of the vbTPM package
% =========================================================================
% 
% Copyright (C) 2013 Martin Linden
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

function x=BWdriftcorrect(x0,Fcut,Fsample)
% design a high-pass Butterworth filter as described by Han et al
if(~exist('Fsample','var')) Fsample=30;end  % sampling frequency (~30 Hz) [Hz]
if(~exist('Fcut','var') ) Fcut=0.05;end     % canonical value used by Phillips group

Wn=Fcut/(Fsample/2); % normalized cut-off frequency
[Z,P,K]=butter(1,Wn,'high'); % signal processing toolbox
[sos,g]=zp2sos(Z,P,K);       % signal processing toolbox?
Hd=dfilt.df2tsos(sos,g);     % Create a dfilt object, signal processing toolbox

% filter indata columns
x=zeros(size(x0));
for k=1:size(x0,2)
    x(:,k)=filter(Hd,x0(:,k));
end
end
