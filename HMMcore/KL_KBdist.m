function KL=KL_KBdist(n,c,mu,v,n0,c0,mu0,v0)
% KL=KL_KBdist(n,c,mu,v,n0,c0,mu0,v0)
%
% compute elementwise Kullback-Leibler divergence 
%
% KL(n,c,mu,v|n0,c0,mu0,v0) = int q*log(q/q0) dKdB

KL=-(n+0.5)./c.*(c-c0-v0.*(mu-mu0).^2)...
   +0.5*log(v./v0)+(n0+0.5).*log(c./c0)...
   -gammaln(n+0.5)+gammaln(n0+0.5)...
   +(n-n0).*psi(n+0.5)+v0./v/2-0.5;

