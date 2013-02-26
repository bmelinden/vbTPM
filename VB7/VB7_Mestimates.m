function est=VB7_Mestimates(M,fSample,downSampling)
% Compute some estimates based on variational distribution
% parametrers (no hyper-parameters substracted).
%
% ML 2012-12-06

% occupation probabilities
est.Pocc=sum(M.wA,2)/sum(sum(M.wA,2));
% transition rates and dwell times: s(t)
est.A=rowNormalize(M.wA);
est.dA=zeros(size(M.wA));
a0=sum(M.wA,2);

for k=1:length(M.mu)
    est.dA(k,:)=sqrt(est.A(k,:).*(1-est.A(k,:))/(1+a0(k)));
end
est.tD=1./(1-diag(est.A))/(fSample/downSampling);

% emission parameters
est.Kaverage=M.mu;
est.Baverage=(M.n+1/2)./M.c;
est.Kstd=sqrt(M.c/2./M.v./(M.n-0.5)); % sqrt(Var(K))
est.Bstd=sqrt(M.n+1/2)./M.c;            % sqrt(Var(B))
est.KBcov=est.Kaverage.*est.Baverage; % <K*B>

est.RMS=sqrt(1./est.Baverage./(1-est.Kaverage.^2));
est.TC=-1./log(M.mu)/(fSample); % real time units; correction for downsampling not needed
est.NC=-1./log(M.mu);                              % sampling time units

