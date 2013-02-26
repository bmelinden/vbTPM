function x=BWdriftcorrect(x0,Fcut,Fsample)
% x=BWdriftcorrect(x0,Fcut,Fsample)
%
% filter the columns of x0 with a hign pass 1st order Butterworth filter
% with cut-off frequency Fcut (default 0.05 Hz).
% Fsample = smapling frequency of x0 (default 30 Hz).

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
