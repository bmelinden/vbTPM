function KL=KL_dirichlet(w,u)
% KL=KL_dirichlet(w,u)
% 
% compute the Kullback-Leibler divergence of two Dirichlet distributions
% with parameters w and u,
% KL = KL(w||u) = int p(x|w) ln ( p(x|w) / p(x|u) ) dx
% 
% In addition, elements with u(j)=0 are ignored, equivalent to 
% u=u(find(u~=0));w=w(find(u~=0));
%
% (In the VBHMM setting, w is the variational distribution, and u is the
% prior.)
%
% M.L. 2012-05-02

% find non-zero indices
ind=find(u~=0);
u=u(ind);
w=w(ind);

w0=sum(w);
u0=sum(u);

KL=gammaln(w0)-gammaln(u0)...
    +sum(gammaln(w)-gammaln(u)...
    +(w-u).*(psi(w)-psi(w0)));





