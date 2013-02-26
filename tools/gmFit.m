function [w,mu,sig,y1,y0]=gmFit(x,y,h0,mu0,sig0)
% [w,mu,sig,y1,y0]=gmFit(x,y,h0,mu0,sig0)
% fit sum of Gaussians to histogram y=freq(x)
% M.L. 2010-06-04

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
