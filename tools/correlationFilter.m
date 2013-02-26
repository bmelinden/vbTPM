function [fcMS,t,RMS]=correlationFilter(x,tSigma,fSample)
% [RMSfc,t]=correlationFilter(x,tSigma,fSample)
% 
% Analoguous to the Gaussian filtered RMS trace for rho^2, one could try
% to plot -fSample*log( <x(t) . x(t-dt) >) ~ <f_c> + const, which would be
% a measure of local correlation time, if <...> means a suitable
% time-average. 
%

%% handle partameters
if(size(x,1)==2 && shape(x,2)>2) % then transpose x
    x=x';
end
N=size(x,1)-1;   % since we compute a correlation function
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
K=sum(x(1:end-1,:).*x(2:end,:),2);  % x(t)*x(t-1)
R2=sum(x(1:end-1,:).^2,2);          % x(t)^2

fr=fft(R2,NFFT);
foo=(ifft(fG(ff).*fr));
RMS=foo(1:N);

fr=fft(K,NFFT);
foo=(ifft(fG(ff).*fr));

% (<x(t)*x(t-1)> / <x(t)^2>) ~ exp(- fc/fSample )
fcMS=-fSample*log(foo(1:N)./RMS);
RMS=sqrt(RMS);


end
