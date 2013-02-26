function W=convert_VBE1toVB2(WE,W20)
%%
% W=VBE1_to_VB2(WE,W20)
% 
% convert a VBE1 structure to the 'closest' VB2 structure, loosely defined
% as the one where the M-step trial distributions are as close as possible.
%
% WE    : VBE1 structure
% W20   : VB2 start structure (only the priors are used)

%% change-log
% M.L. 2010-07-19   : start

%% actual computation
W.PM=W20.PM;
W.N=W20.N;
W.a=W20.a;

W.M.wPi=WE.M1.wPi-WE.PM1.wPi+W.PM.wPi;
W.M.m=WE.M1.m-WE.PM1.m+W.PM.m;
W.M.b=WE.M1.b-WE.PM1.b+W.PM.b;
W.M.wA=WE.M1.wA;

W.M.n=WE.M2.n-WE.PM2.n+W.PM.n;
W.M.c=WE.M2.c-WE.PM2.c+W.PM.c;

W.mu=0*W20.mu0;
W.mu(1,:)=WE.Kaverage;
for j=1:W.N
    W.Sigma(:,:,j)=eye(length(W.a))*1e-10;
    W.Sigma(1,1,j)=WE.Kvariance(j)*W.M.c(j)/2/(W.M.n(j)-1-length(W.a)/2);
    W.M.V(:,:,j)=inv(W.Sigma(:,:,j));
    W.M.w(:,j)=W.M.V(:,:,j)*W.mu(:,j);
end

