function [n0,c0,v0,mu0]=VB7_EBupdate_KB(n,c,v,mu)
% [n0,c0,v0,mu0]=VB7_EBupdate_KB(n,c,v,mu)
% perform the hierarchical prior update steps on the KB-distribution in the
% Gaussian TPM model, with emission probability parameters n,c,v,mu.

% M.L. 2012-05-02

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

% choose smallest residual:
% if this is necessary, then there is probably something wrong with the
% input models
%for m=1:length(lnb0)
%        if(abs(res1(m))<abs(res(m)))
%            n0(m)=n1(m);
%        end
%        
%    %if(n0(m)>0 && n1(m)>0)%
%
%    elseif(n1(m)>0 && n0(m)<0)
%        n0(m)=n1(m);
%    elseif(n0(m)<0 && n1(m) <0)
%       save foo
%       error('negative n found for both equations')
%    end
%end

c0=(n0+0.5)./a0;




if(0) % numerical refinement: this actually seems to make things worse...
    
    for m=1:size(n,2)
        %disp(['state : ' int2str(m)])
        f=@(x)(sum(KL_KBdist(n(:,m),c(:,m),mu(:,m),v(:,m),x(1),x(2),x(3),x(4))));
        x0=[n0(m) c0(m) mu0(m) v0(m)];
        x=fminsearch(f,x0,opt2);
        n1(m) =x(1);
        c1(m) =x(2);
        mu1(m)=x(3);
        v1(m) =x(4);
        
        df(m)=f(x)-f(x0);
        disp(['state : ' int2str(m) ', [dF dF/|F| ] = ' num2str([df(m) df(m)/abs(f(x))])])
        disp('----------')
    end
    
    disp(['dn : ' num2str(n1./n0-1)])
    disp(['dc : ' num2str(c1./c0-1)])
    disp(['dmu: ' num2str(mu1./mu0-1)])
    disp(['dv : ' num2str(v1./v0-1)])
end
