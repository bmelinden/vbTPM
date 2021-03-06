%% [RMSf,fRMS,t]=gaussFilter(x,tSigma,tSample)
% Gaussian lowpass filter with window std tSigma, or a signal with sampling
% time tSample (default 1/30 s). The radius rho(t) is defined as
% rho(t)^2=x(t)^2+y(t)^2. If <.> denotes the action of the Gaussian
% smoothing, then 
% RMSf=sqrt(<rho(t)^2>)  RMS of filtered rho^2-signal
% fRMS=<sqrt(rho(t)^2)>  filtered sqrt(rho^2)-signal
% The present analysis uses RMSf, which tends to have slightly larger
% amplitudes and smaller end effects. 
%
% x =[x(t) y(t)]    : two T by 2 matrix of positions
% tSigma            : std of Gaussian smoothing kernel, default 4 [s].
% fSample           : sampling frequency, default 30 [Hz].
%
% This algorithm uses straight multiplication in frequency space, computed
% with fft, with no bells and whistles to handle things like end effects
% etc. Hence, the signal falls off rapidly in the first and last tSigma seconds. 

%% change-log
% M.L. 2010-07-06   : added default value (4 s) for tSigma.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gaussFilter.m, part of the vbTPM package
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

function [RMSf,fRMS,t,r2]=gaussFilter(x,tSigma,fSample)
%% handle partameters
if(size(x,1)==2 && shape(x,2)>2) % then transpose x
    x=x';
end
N=size(x,1);
if(~exist('fSample','var'))
    fSample=30; % Hz
end
if(~exist('tSigma','var'))
    tSigma=4; % s
end
%% Gaussian filter parameters
fc=sqrt(log(2))/2/pi/tSigma; % 3dB cut-off frecuency [Hz]
fG=@(f)(exp(-log(2)/2*(f./fc).^2));
NFFT=2^nextpow2(N);                 % number of frequency points in fft
f=fSample/2*linspace(0,1,NFFT/2+1); % spectral frequcies [Hz]
ff=[f -f(end-1:-1:2)]';             % frequencies in fft
t=(N-1)/fSample*linspace(0,1,N);    % time points
%% apply Gaussian filter in Fourier space
r2=sum(x.^2,2); % rho^2(t)

fr2=fft(r2,NFFT);
foo=sqrt(ifft(fG(ff).*fr2));
RMSf=foo(1:N); % remove spurious points

fr= fft(sqrt(r2),NFFT);
foo=ifft(fG(ff).*fr);
fRMS=foo(1:N); % remove spurious points


