function [RMS,kG,bG,t]=RMSKBgaussfilter(x,tSigma,fSample)
% [RMS,kG,bG,t]=RMSKBgaussfilter(x,tSigma,fSample)
%
% a running average estimate of RMS,K,B, and their local std.
% kG~<K(t)>, bG~<B(t)>, averaged using a gaussian kernel of width tSigma s.
%
% fSample is the sampling frequency.
% default: tSigma=4 s, fSample=30 Hz.
%

%% change-log
% 2010-10-18 M.L.   : found factor 2 bug in B(t), and verified the
%                     algorithm numerically with synthetic data.
%                     Changed RMS length to coincide with the other
%                     properties
%% parameters
if(~exist('fSample','var')|| isempty(fSample)); fSample=30;end    
if(~exist('tSigma','var')|| isempty(tSigma)); tSigma=4;end
if(size(x,1)==1)
    RMS=[];kG=[];bG=[];t=[];dk=[];db=[];
    return
end
%% Gaussian filter parameters
N=size(x,1);
fc=sqrt(log(2))/2/pi/tSigma; % 3dB cut-off frecuency [Hz]
fG=@(f)(exp(-log(2)/2*(f./fc).^2));
NFFT=2^nextpow2(N);                 % number of frequency points in fft
f=fSample/2*linspace(0,1,NFFT/2+1); % spectral frequcies [Hz]
ff=[f -f(end-1:-1:2)]';             % frequencies in fft
t=(N-1)/fSample*linspace(0,1,N);    % time points
%% apply Gaussian filter in Fourier space
r2=sum(x.^2,2); % rho^2(t)
fr2=fft(r2,NFFT);
r2=ifft(fG(ff).*fr2);
RMS=sqrt(r2(1:N));
S1=r2([1 1:N-1]); % r2(t-1)
S2=r2(1:N);   % r2(t)

S12=sum(x([1 1:end-1],:).*x([2 2:end],:),2); % sum X(t)*x(t-1), double-count t=1,2
fS12=fft(S12,NFFT);
S12=ifft(fG(ff).*fS12);
S12=S12(1:N);
% n=ones(size(S1));

c=(S2-S12.^2./S1);
v=S1./c;

kG=S12./S1;
dk=sqrt(0.5./v);
bG=1./c;
db=bG;





