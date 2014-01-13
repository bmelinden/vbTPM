% [w,mu,sig,y1,y0]=gmFit(x,y,h0,mu0,sig0)
% fit sum of Gaussians to histogram y=freq(x)
% M.L. 2010-06-04

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gmFit.m, part of the vbTPM package
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

function [w,mu,sig,y1,y0]=gmFit(x,y,h0,mu0,sig0)
M=length(mu0); % number of peaks

hms=lsqcurvefit(@H,[h0 mu0 sig0],x,y);
h=hms(1:M);
mu=hms(M+1:2*M);
sig=hms(2*M+1:3*M);
y0=H([h0 mu0 sig0],x);
y1=H(hms,x);

w=zeros(1,M);
for k=1:M
   yk= H([h(k) mu(k) sig(k)],x);
   w(k)=sum(yk);
end

figure(33)
clf, hold on
plot(x,y,'k')
plot(x,y0,'g')
plot(x,y1,'r')
end

function h=H(wms,X)
h=0*X;
K=length(wms)/3;
if(K~=round(K))
    error('argument to H must have length be multiple of three')
end
w=wms(1:K);
m=wms(K+1:2*K);
s=wms(2*K+1:3*K);
for kk=1:length(m)
    h=h+w(kk)*exp(-0.5*((X-m(kk))/s(kk)).^2);
end
end
