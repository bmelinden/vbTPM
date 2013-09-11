function [n0,c0,v0,mu0]=VB7_EBupdate_KB(n,c,v,mu)
% [n0,c0,v0,mu0]=VB7_EBupdate_KB(n,c,v,mu)
% perform the hierarchical prior update steps on the KB-distribution in the
% Gaussian TPM model, with emission probability parameters n,c,v,mu.
% PM is an optional prior structure (with fields n,c,mu,v) used as initial
% guess in the numerical optimization procedure.


% ML 2012-07-04 : inserted more check-points, caught output bug, and inactivated the numerical
% refinement, which only seems to make things worse.
% ML 2012-07-04 : switched to a log-equation for n, to avoid x<0 errors in
% the numerical solution

opt1=optimset('Display','off','TolFun',1e-15,'TolX',1e-15,'MaxFunEvals',1e5,'MaxIter',1e5);
opt2=optimset('Display','on','TolFun',1e-15,'TolX',1e-15,'MaxFunEvals',1e5,'MaxIter',1e5);
M=size(mu,1);

% analytical solution
a=(n+0.5)./c;
lnb=psi(n+0.5)-log(n+0.5);

mu0=sum(mu.*a,1)./sum(a,1);


dmu2=(mu-ones(M,1)*mu0).^2;
VS=1/2./v+a.*dmu2;
v0=1/2./(sum(VS,1)/M);

a0=mean(a,1);
lnb0=mean(log(a)+lnb,1)-log(a0);
n0=min(n,[],1);
res=0*lnb0;
for m=1:length(lnb0)
    % try using a logarithmic unknown to enforce exp(x)>0
    leq=@(x)(psi(exp(x)+0.5)-log(exp(x)+0.5)-lnb0(m));  
    
    x0=log(0.1*n0(m));
    if(~isreal(x0))
        disp('im initial guess')
        x0=-10;
    end    
    [logn,~,nflag,noutput]=fsolve(leq,x0,opt1);
    n1(m)=exp(logn);
    res1(m)=leq(logn);

    % uning a linear unknown
    eq=@(x)(psi(x+0.5)-log(x+0.5)-lnb0(m));
    try
        [n0(m),~,nflag,noutput]=fsolve(eq,0.1*n0(m),opt1);
        res(m)=eq(n0(m));
        
        if(nflag<=0)
            disp(['VB7_EBupdate_KB fsolve output flag, linear eq: ' int2str(nflag)])
            noutput
            disp(['message : ' noutput.message])
            disp('--------------------')
        end
    catch me
        me
    	disp('could not solve KB equations linear in n. going w log equations only.')
        res(m)=inf;
        n0(m)=n1(m);
    end
    
   
    
end
dnrel=(n1-n0)./n1;
disp([' dn/n ' num2str(dnrel,3)])
disp([' nres lin:  ' num2str(res,3)])
disp([' nres log: ' num2str(res1,3)])

%if(max(abs(dnrel))>0.1)
%    keyboard
%end

% generally use the logarithmic equation
n0=n1;

c0=(n0+0.5)./a0;
