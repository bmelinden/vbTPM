function W=VB7_VBEMiter(W,X,varargin)
%% W=VB7_VBEMiter(W,X,...)
%
% Perform a full VBEM iteration, with the K,B-correated Gaussian noise
% model, on the VB7 structure W, with data structure X.
% X=VB3_preprocess(x,downSample);
% 
% If W.E and W.Ec exist, then a full VBE iteration (M-step + E-step) is
% performed. If W.E or W.Ec fields are missing, then only the E-step is
% performed. This is useful for starting up an iteration, since initial
% guesses are only needed for the emission parameters, or to restore an
% estimate if only the M-step fields were stored, for example using the
% function VB7_minimalStorage().
%
% possible options:  'estimate' adds a field est2, which contains memory
% intensive estiamtes, such as the viterbi path.
%
% Fields in the VB7 tructure W:
% W.M.wPi;         : initial state
% W.M.wA;         : transition counts for s(t) transition matrix
% W.Mc.wA;        : transition counts for onset of dirt events (c(t) >1)
% W.Mc.wR;          : transition counts for ending of dirt events
% W.M.n; W.M.c;  : B-distribution for s(t)
% W.M.v; W.M.mu; : K-distribution for s(t)
% W.Mc.n; W.Mc.c;  : B-distribution for c(t)
% W.Mc.v; W.Mc.mu; : K-distribution for c(t)
% Prior parameters are stored in W.PMs and W.PMc
% E-step parameters are under W.E
% W.F              : lower bound on the likelihood
% W.est            : various useful estimates that do not take up a lot of
%                    memory 
% W.est2           : more memory intensensive estimates
%
% this function uses mex files for the computer intensive inner loops:
% ../HMMcore/VBforwback.c and ../HMMcore/VBwAcount.c
% If you have no binaries for your system, try to compile them with the
% compilation script in ../HMMcore/HMMcore_compile.m, or send an e-mail to
% bmelinden@gmail.com

%% change log and notes
% M.L. 2012-01-12   : moved inner loop binaries (VBwaCount and VBforwback)
%                     to HMMcore directory, to make it easier to reuse them
%                     in other HMM projects.
% M.L. 2011-12-20   : stored log(qst) and log(Q) in W.est, to try and hunt
%                     down weird Viterbi paths for really long
%                     trajectories. Did not help, I think weird Viterbi
%                     paths happens when a combination of very long
%                     trajectories and weak transition probability priors
%                     lead to numerical underflow in the VBE iterations (it
%                     still runs though, and the Pmax paths look OK).
% M.L. 2011-02-02   : started from VB5_VBEMiter.m
% M.L. 2011-02-21   : ~20% optimization by skipping some steps in the
%                     transitions count loop when the model has no spurious
%                     states.
% M.L. 2011-02-25   : Further optimizations by rewriting large parts of the
%                     VBE step in c.
% M.L. 2011-03-01   : successfully renamed E.C -> E.U, E.S -> E.C E.M ->
%                     E.V, E.P -> E.M, and made corresponding changes in
%                     the Ec field names. 
% M.L. 2011-07-25   : removed the internal function minimalStorage from being saved to the VB7 structure, since
%                     it turns out to be very bulky to store. Created an
%                     external VB7 function instead.
% M.L. 2011-08-15   : reformulated VBE step slightly to avoid division by
%                     zero in case of points x=(0,0).

%% options
do_estimates=false;
if(nargin>2) % then parse options
    kmax=nargin-2;
    k=1;
    while(k<=kmax)
        option=varargin{k};
        if(strcmpi(option,'estimate'))
            do_estimates=true;
        elseif(strcmpi(option,'purge'))
            if(isfield(W,'est'))
                W=rmfield(W,'est');
            end
            if(isfield(W,'est2'))
                W=rmfield(W,'est2');
            end
        else
            error(['VB%_VBEMiter: option ' option ' not recognized.'])
        end
    k=k+1;
    end
end
%% sanity check on input
if(~isempty(find((W.Mc.wA(W.PMc.wA==0)>0), 1)))
   warning('VB7:zeropattern','VB7_VBEMiter.m : W.Mc.wA and  W.PMc.wA have different zero patterns. Check initial guess.')
end
W.T=size(X.R2,1);
NN=W.N*W.Nc;  % number of states in composite hidden process
%% deal with VB7 object with old notation
if( isfield(W,'E') && isfield(W.E,'P') ) % then this was created before the name change
    W.E.U=W.E.C;
    W.E.C=W.E.S;
    W.E.V=W.E.M;
    W.E.M=W.E.P;
    W.E=rmfield(W.E,{'S','P'});
    W.Ec.U=W.Ec.C;
    W.Ec.C=W.Ec.S;
    W.Ec.V=W.Ec.M;
    W.Ec.M=W.Ec.P;    
    W.Ec=rmfield(W.Ec,{'S','P'});
end
%% start of actual code
%% state convention z(ht)=( c(t) , s(t) ), in the order 
% (1,1), (1,2), ... , (1,N), (2,1), ... , (2,N), (3,1) , ... (Nh,N)
    function z=MZ(c,s)
        % translate from (c,s) to internal state index
        z=(c-1)*W.N+s;
    end
sMap=zeros(1,NN);
cMap=zeros(1,NN);
for j=1:W.N
    for k=1:W.Nc
        cMap(MZ(k,j))=k;
        sMap(MZ(k,j))=j;
    end
end
%% M-step
if(isfield(W,'E') && isfield(W,'Ec')) % then start with M step
    %% transitions and initial conditions
    W.M.wPi = W.PM.wPi + W.E.ds_1;
    W.Mc.wPi= W.PMc.wPi+W.Ec.ds_1;
    
    W.M.wA =  W.PM.wA + W.E.wA;
    W.Mc.wA= W.PMc.wA +W.Ec.wA;    
    W.Mc.wR= W.PMc.wR +W.Ec.wR;    
    %% emission model
    W.M.n  = W.PM.n  + W.E.M;
    W.M.v  = W.PM.v  + W.E.V;
    W.M.mu = (W.PM.v.*W.PM.mu + W.E.U)./(W.PM.v+W.E.V);
    W.M.c  = W.PM.c  + W.E.C + W.PM.v.*W.PM.mu.^2 - W.M.v.*W.M.mu.^2;

    W.Mc.n(2:end)  = W.PMc.n(2:end)  + W.Ec.M(2:end);
    W.Mc.v(2:end)  = W.PMc.v(2:end)  + W.Ec.V(2:end);
    W.Mc.mu(2:end) =(W.PMc.v(2:end).*W.PMc.mu(2:end) + W.Ec.U(2:end))./(W.PMc.v(2:end)+W.Ec.V(2:end));
    W.Mc.c(2:end)  = W.PMc.c(2:end)  + W.Ec.C(2:end) + W.PMc.v(2:end).*W.PMc.mu(2:end).^2 - W.Mc.v(2:end).*W.Mc.mu(2:end).^2;
        
    % check for problems
    isNanInf=(sum(~isfinite([W.M.wPi W.Mc.wPi W.M.wA(1:end) reshape(W.Mc.wA,1,NN) ...
        W.M.n  W.Mc.n W.M.c W.Mc.c W.M.v W.Mc.v W.M.mu W.Mc.mu]))>1);
    if(isNanInf)
        error('VB7_VBEMiter:Mfield_not_finite','Nan/Inf generated in VBM step')
    end
end
%% E-step starts here
%% trial emission priobability q(s,t)
%% variational pointwise contributions, same as VB4
% t=1
qstL1=zeros(1,NN);
for j=1:W.N
    for k=1:W.Nc
        qstL1(1,MZ(k,j))=psi( W.M.wPi(j))-psi(sum( W.M.wPi)) ...
                        +psi(W.Mc.wPi(k))-psi(sum(W.Mc.wPi));
    end
end
% t>1: we compute log(qst) first, but use the same variable to save memory.
qst   =zeros(W.T,NN);
qstL0s=psi(W.M.n+1/2 )-log(pi*W.M.c); % same for all t>1
qstL0c=psi(W.Mc.n+1/2)-log(pi*W.Mc.c);% same for all t>1
% k=1 :
for j=1:W.N
    qst(:,j)=qstL0s(j)-X.R2m1/2/W.M.v(j)...
        -(W.M.n(j)+1/2)/W.M.c(j)*(...
        X.R2+W.M.mu(j)^2*X.R2m1-2*W.M.mu(j)*X.X12);
end
for k=2:W.Nc
    qst_k=qstL0c(k)-X.R2m1/2/W.Mc.v(k)...
        -(W.Mc.n(k)+1/2)/W.Mc.c(k)*(...
        X.R2+W.Mc.mu(k)^2*X.R2m1-2*W.Mc.mu(k)*X.X12);
    for j=1:W.N
        qst(:,MZ(k,j))=qst_k;
    end
end
qst=qst*X.downSampling;
% the down-sampling factor should not affect t=1, which is not blocked (and
% does not depend on the data anyway). Hence, we reinstert our original
% result here.
qst(1,:)=qstL1; 
qstLMax=max(qst,[],2);
lnqst =zeros(W.T,NN);
for jk=1:NN % now we compute the actual pointwise emission contribution
    lnqst(:,jk)=qst(:,jk)-qstLMax;
    %qst(:,jk)=exp(qst(:,jk)-qstLMax);
end
qst=exp(lnqst);
%% variational transition probabilities (lnQcc changed in VB7)
% first assemble the individual transition matrices
lnQss=zeros(W.N,W.N);
for j=1:W.N
    for jm=1:W.N
        lnQss(jm,j)=psi(W.M.wA(jm,j))-psi(sum(W.M.wA(jm,:)));
    end
end
lnQsc=zeros(W.N,W.Nc);
for j=1:W.N
    for k=1:W.Nc
        lnQsc(j,k)=psi(W.Mc.wA(j,k))-psi(sum(W.Mc.wA(j,:)));
    end
end
lnQcc=zeros(W.Nc,W.Nc); % all transitions are allowed in VB7
for j=2:W.Nc
    for k=1:W.Nc
        lnQcc(j,k)=psi(W.Mc.wR(j,k))-psi(sum(W.Mc.wR(j,:)));
    end
end
% then build the composite transition matrix
lnQ=zeros(NN,NN);
for jm=1:W.N
    for j=1:W.N
        for km=1:W.Nc
            for k=1:W.Nc
                if(km==1) % then include sticking rate
                    lnQ(MZ(km,jm),MZ(k,j))=lnQss(jm,j)+lnQsc(jm,k);
                else      % then include unsticking rate
                    lnQ(MZ(km,jm),MZ(k,j))=lnQss(jm,j)+lnQcc(km,k);                    
                end                
            end
        end
    end
end
lnQmax=max(max(lnQ));
Q=exp(lnQ-lnQmax);
%% forward sweep  (depends only on Q and qst)
[Za,alpha,~,beta,qt]=VBforwback(Q,qst);
%% compute quantities for next M-steps, same as VB5
% ML: problems with numerical underflow in Q or qst can manifest itself as
% zero elements in qt. a cheap fix (not implemented yet), is to add a very
% small number ( sqrt(realmin) ?) to all elements.
%qt=alpha.*beta;  % qt(t,j) = <(z(t)> = P(z(t)=j)
%qtZ=sum(qt,2);
%isNanInf=(~isempty(find(qtZ==0,1))); % check that qt is normalizable 
%if(isNanInf)
%        error('VB7_VBEMiter:qt_not_normalizeable','non-overlapping alpha(t,:) and beta(t,:) generated in VBE step.')
%end
%for j=1:NN
%   qt(1:end,j)=qt(1:end,j)./qtZ(1:end,1); % add realmin to avoid qt=0? 
%end
%% transition counts
% helper matrices
Ms=zeros(NN,W.N); % integrates out (c,s) -> s 
for k=1:W.Nc 
   Ms(MZ(k,1:W.N),1:W.N)=eye(W.N);
end
Mc=zeros(NN,W.Nc); % integrates out (c,s) -> c 
for k=1:W.Nc 
    Mc(MZ(k,1:W.N),k)=1;
end
MsT=Ms'; % pre-compute matrix transposes to save time
McT=Mc'; 

W.E.ds_1 =qt(1,:)*Ms;          % <s_1>=<\delta_{j,s_1}>_{q(s)}
W.Ec.ds_1=qt(1,:)*Mc;

W.E.wA =zeros( W.N,W.N );
W.Ec.wA=zeros( W.N,W.Nc);
W.Ec.wR=zeros(W.Nc,W.Nc);
WA=VBwAcount(alpha,beta,qst,Q);
if(W.Nc>1)
    W.E.wA=MsT*WA*Ms;    % all s(t-1) -> s(t) transitions
    W.Ec.wR=McT*WA*Mc;  % all c(t-1) -> c(t) transitions
    W.Ec.wA=WA(1:W.N,:)*Mc;% only include transtion out of c(t-1)=1
    W.Ec.wR(1,:)=0; % remove transitions out of c(t-1)=1    
else % speedup if there are only genuine states
    W.E.wA=WA;
    W.Ec.wA=sum(W.E.wA,2);
end
%% for the emission models (same as VB4)
W.E.M=sum(qt(2:end,MZ(1,1:W.N)),1);   % sum_{t=2}^T p(s_t,c_t==1)
W.E.C=zeros(1,W.N); % W.E.C = sum_{t=2}^T P(s(t)=j,c(t)=1).*x(t)^2
W.E.V=zeros(1,W.N); % W.E.V = sum_{t=2}^T P(s(t)=j,c(t)=1).*x(t-1)^2
W.E.U=zeros(1,W.N); % W.E.U = sum_{t=2}^T P(s(t)==j)*(x(t)*x(t-1))
for j=1:W.N
    W.E.C(j)=sum(qt(2:end,j).*X.R2(2:end,1),1);
    W.E.V(j)=sum(qt(2:end,j).*X.R2m1(2:end,1),1);    
    W.E.U(j)=sum(qt(2:end,j).*X.X12(2:end,1));
end

W.Ec.M=zeros(1,W.Nc);
W.Ec.C=zeros(1,W.Nc); % W.E.C = sum_{t=2}^T P(c(t)=k).*x(t)^2
W.Ec.V=zeros(1,W.Nc); % W.E.V = sum_{t=2}^T P(c(t)=k).*x(t-1)^2
W.Ec.U=zeros(1,W.Nc);
for k=2:W.Nc
    W.Ec.M(k)=sum(sum(qt(2:end,MZ(k,1:W.N)),1));
    W.Ec.C(k)=sum(sum(qt(2:end,MZ(k,1:W.N)),2).*X.R2(2:end,1),1);
    W.Ec.V(k)=sum(sum(qt(2:end,MZ(k,1:W.N)),2).*X.R2m1(2:end,1),1);
    W.Ec.U(k)=sum(sum(qt(2:end,MZ(k,1:W.N)),2).*X.X12(2:end,1),1);
end

% account for down sampling
W.E.M=W.E.M*X.downSampling;
W.E.C=W.E.C*X.downSampling;
W.E.V=W.E.V*X.downSampling;
W.E.U=W.E.U*X.downSampling;
W.Ec.M=W.Ec.M*X.downSampling;
W.Ec.C=W.Ec.C*X.downSampling;
W.Ec.V=W.Ec.V*X.downSampling;
W.Ec.U=W.Ec.U*X.downSampling;
% check for problems
isNanInf=(sum(~isfinite([W.E.M W.E.C W.E.V W.E.U]))>1);
if(isNanInf)
    error('VB3_VBEMiter:Efield_not_finite','Nan/Inf generated in VBE step')
end
%% compute free energy
%% forward sweep normalization constant (same as VB4)
lnZQ=(W.T-1)*lnQmax;
lnZq=sum(qstLMax);
lnZz=sum(log(Za));
lnZ=lnZQ+lnZq+lnZz;
W.estF.Zterms=[lnZQ lnZq lnZz];
%% KL divergence of transition probabilities of s(t) (same as VB4)
KL_Ak=zeros(1,W.N); 
for k=1:W.N
    wAk0=sum(W.M.wA(k,:));
    uAk0=sum(W.PM.wA(k,:));
    nz=find(W.PM.wA(k,:)>0);
%    KL_Ak(k)=gammaln(wAk0)-gammaln(uAk0)-(wAk0-uAk0)*psi(wAk0)...
%        -sum(gammaln(W.M.wA(k,:))-gammaln(W.PM.wA(k,:))...
%        -(W.M.wA(k,:)-W.PM.wA(k,:)).*psi(W.M.wA(k,:)));
    KL_Ak(k)=gammaln(wAk0)-gammaln(uAk0)-(wAk0-uAk0)*psi(wAk0)...
        -sum(gammaln(W.M.wA(k,nz))-gammaln(W.PM.wA(k,nz))...
        -(W.M.wA(k,nz)-W.PM.wA(k,nz)).*psi(W.M.wA(k,nz)));
end
KL_A=sum(KL_Ak);
W.estF.Aterms=-KL_Ak;
clear wAk0 uAk0 nz;
%% KL divergence of transition probabilities of c(t)
% sticking rates, same as VB5
KL_Ak=zeros(1,W.N); 
for k=1:W.N
    wAk0=sum(W.Mc.wA(k,:));
    uAk0=sum(W.PMc.wA(k,:));
    KL_Ak(k)=gammaln(wAk0)-gammaln(uAk0)-(wAk0-uAk0)*psi(wAk0)...
        -sum(gammaln(W.Mc.wA(k,:))-gammaln(W.PMc.wA(k,:))...
        -(W.Mc.wA(k,:)-W.PMc.wA(k,:)).*psi(W.Mc.wA(k,:)));
end
KL_Ac=sum(KL_Ak(1:W.N));
clear wAk0 uAk0;
% unsticking rates, same as VB5, as long as W.PMc.wR has correct non-zero
% pattern
KL_Ak=zeros(1,W.Nc); 
for k=2:W.Nc
    cCols=W.PMc.wR(k,:)>0; % only count transitions where PMc.wR>0
    wRk0=sum(W.Mc.wR(k,cCols));
    uRk0=sum(W.PMc.wR(k,cCols));
    KL_Ak(k)=gammaln(wRk0)-gammaln(uRk0)-(wRk0-uRk0)*psi(wRk0)...
        -sum(gammaln(W.Mc.wR(k,cCols))-gammaln(W.PMc.wR(k,cCols))...
        -(W.Mc.wR(k,cCols)-W.PMc.wR(k,cCols)).*psi(W.Mc.wR(k,cCols)));
end
KL_R=sum(KL_Ak(2:W.Nc));
W.estF.Acterms=[-KL_Ak -KL_R];
%% KL divergence of initial state probability for s(t),c(t) (same as VB4)
u0Pi=sum(W.PM.wPi);
KL_pi = log(u0Pi)-psi(u0Pi)-1/u0Pi...
    +sum((W.M.wPi-W.PM.wPi).*psi(W.M.wPi)...
    -gammaln(W.M.wPi)+gammaln(W.PM.wPi));
u0cPi=sum(W.PMc.wPi);
KL_pih = log(u0cPi)-psi(u0cPi)-1/u0cPi...
    +sum((W.Mc.wPi-W.PMc.wPi).*psi(W.Mc.wPi)...
    -gammaln(W.Mc.wPi)+gammaln(W.PMc.wPi));
W.estF.piTerms=[-KL_pi -KL_pih];
%% KL divergence of emission parameters for s(t) (same as VB4)
mKL_BKv=zeros(1,W.N);
mKL_BKc=zeros(1,W.N);
mKL_BKn=zeros(1,W.N);
for j=1:W.N
    mKL_BKn(j)=gammaln(W.M.n(j)+1/2)-gammaln(W.PM.n(j)+1/2)...
               -(W.M.n(j)-W.PM.n(j))*psi(W.M.n(j)+1/2);
    mKL_BKc(j)=(W.M.n(j)+1/2)/W.M.c(j)*(W.M.c(j)-W.PM.c(j))...
              -(W.PM.n(j)+1/2)*log(W.M.c(j)/W.PM.c(j));
    mKL_BKv(j)=...
        +(W.M.n(j)+1/2)/W.M.c(j)*(...
        -W.PM.v(j)*(W.M.mu(j)-W.PM.mu(j))^2) ...
        -0.5*log(W.M.v(j)/W.PM.v(j))...
        -W.PM.v(j)/W.M.v(j)/2+1/2;
end
KL_BK=-sum(mKL_BKv+mKL_BKc+mKL_BKn); % this is the KL divergence
W.estF.BKterms=[sum(mKL_BKn) sum(mKL_BKc) sum(mKL_BKv)];
%% KL divergence of emission parameters for c(t) (same as VB4)
mKL_BKv=zeros(1,W.Nc);
mKL_BKc=zeros(1,W.Nc);
mKL_BKn=zeros(1,W.Nc);
for k=2:W.Nc
    mKL_BKn(k)=gammaln(W.Mc.n(k)+1/2)-gammaln(W.PMc.n(k)+1/2)...
               -(W.Mc.n(k)-W.PMc.n(k))*psi(W.Mc.n(k)+1/2);
    mKL_BKc(k)=(W.Mc.n(k)+1/2)/W.Mc.c(k)*(W.Mc.c(k)-W.PMc.c(k))...
              -(W.PMc.n(k)+1/2)*log(W.Mc.c(k)/W.PMc.c(k));
    mKL_BKv(k)=...
        +(W.Mc.n(k)+1/2)/W.Mc.c(k)*(...
        -W.PMc.v(k)*(W.Mc.mu(k)-W.PMc.mu(k))^2) ...
        -0.5*log(W.Mc.v(k)/W.PMc.v(k))...
        -W.PMc.v(k)/W.Mc.v(k)/2+1/2;
end
KL_BKc=sum(-mKL_BKv-mKL_BKc-mKL_BKn);
W.estF.BcKcterms=-[sum(mKL_BKn) sum(mKL_BKc) sum(mKL_BKv)];
%% assembly of the free energy (same as VB5)
W.estF.Fterms=[lnZ -KL_BK -KL_A -KL_pi -KL_BKc -KL_Ac -KL_R -KL_pih];
W.F=sum(W.estF.Fterms);
if(~isfinite(W.F))
    error('VB3_VBEMiter:F_not_finite','Nan/Inf generated in F')
end
%% do estimates 
%% light-weight estimates (every time)
%% state maps
W.est.cMap=cMap;
W.est.sMap=sMap;
W.est.Ms=Ms;
W.est.Mc=Mc;
%% occupation probabilities
W.est.Paverage=sum(qt,1);
W.est.Paverage=W.est.Paverage/sum(W.est.Paverage);

W.est.cAverage=W.est.Paverage*W.est.Mc;
W.est.sAverage=W.est.Paverage*W.est.Ms;

W.est.sVisible=sum(qt(1:end,1:W.N),1); % average p(s) excluding dirt states
W.est.sVisible=W.est.sVisible/sum(W.est.sVisible);
%% transition rates and dwell times: s(t)
W.est.Q=Q;
W.est.lnQ=lnQ;
W.est.lnQss=lnQss;
W.est.lnQsc=lnQsc;
W.est.lnQcc=lnQcc;
W.est.A=rowNormalize(W.M.wA-W.PM.wA);
W.est.dA=zeros(size(W.M.wA));
a0=sum(W.M.wA,2);
for k=1:W.N
    W.est.dA(k,:)=sqrt(W.est.A(k,:).*(1-W.est.A(k,:))/(1+a0(k)));
end
%W.est.tD=diag(W.est.A)./(1-diag(W.est.A).^2)/(X.fSample0/X.downSampling);
W.est.tD=1./(1-diag(W.est.A))/(X.fSample0/X.downSampling);
%% transition rates and dwell times: c(t)
W.est.Ac=rowNormalize(W.Mc.wA-W.PMc.wA);
W.est.dAc=zeros(size(W.Mc.wA));
ah0=sum(W.Mc.wA,2);
for k=1:W.N
    W.est.dAc(k,:)=sqrt(W.est.Ac(k,:).*(1-W.est.Ac(k,:))/(1+ah0(k)));
end

W.est.Rc=rowNormalize(W.Mc.wR-W.PMc.wR);
W.est.Rc(1,:)=0; % c=1 does not have an unsticking rate
W.est.dRc=zeros(size(W.Mc.wR));
rh0=sum(W.Mc.wR,2);
for k=1:W.Nc
    W.est.dRc(k,:)=sqrt(W.est.Rc(k,:).*(1-W.est.Rc(k,:))/(1+rh0(k)));
end
%W.est.tStick=W.est.Ac(:,1)./(1-W.est.Ac(:,1).^2)/(X.fSample0/X.downSampling);
%W.est.tUnstick=diag(W.est.Rc)./(1-diag(W.est.Rc).^2)/(X.fSample0/X.downSampling);
W.est.tStick=1./(1-W.est.Ac(:,1))/(X.fSample0/X.downSampling);
W.est.tUnstick=1./(1-diag(W.est.Rc))/(X.fSample0/X.downSampling);
W.est.tUnstick(1)=0; % c=1 does not have an unsticking time

W.est.csRates=1./W.est.tUnstick;
W.est.csRares(1)=0; % c=1 does not have an unsticking rate
kMc=-log(W.est.Ac(:,1));
W.est.scRates=zeros(W.N,W.Nc);
W.est.scRates(:,1)=-X.fSample0/X.downSampling*kMc;
for k=1:W.N
    W.est.scRates(k,2:end)=kMc(k)*W.est.Ac(k,1)*(1-W.est.Ac(k,2:end));
end
W.est.scRates=X.fSample0/X.downSampling*W.est.scRates;

% emission parameters
W.est.sKaverage=W.M.mu;
W.est.sBaverage=(W.M.n+1/2)./W.M.c;
W.est.sKstd=sqrt(W.M.c/2./W.M.v./(W.M.n-0.5)); % sqrt(Var(K))
W.est.sBstd=sqrt(W.M.n+1/2)./W.M.c;            % sqrt(Var(B))
W.est.sKBcov=W.est.sKaverage.*W.est.sBaverage; % <K*B>

W.est.sRMS=sqrt(1./W.est.sBaverage./(1-W.est.sKaverage.^2));
W.est.sTC=-1./log(W.M.mu)/(X.fSample0/X.downSampling); % real time units
W.est.sNC=-1./log(W.M.mu);                              % sampling time units

W.est.cKaverage=[0 W.Mc.mu(2:end)];
W.est.cBaverage=[0 ... %sum(W.est.sBaverage.*W.est.sAverage) ... % better plot if Bc(1)>0
    (W.Mc.n(2:end)+1/2)./W.Mc.c(2:end)];
W.est.cKstd=[0 sqrt(W.Mc.c(2:end)/2./W.Mc.v(2:end)./(W.Mc.n(2:end)-0.5))];
W.est.cBstd=[0 sqrt(W.Mc.n(2:end)+1/2)./W.Mc.c(2:end)];
W.est.cKBcov=W.est.cKaverage.*W.est.cBaverage; % <K*B>

W.est.cRMS=[0 sqrt(1./W.est.cBaverage(2:end)./(1-W.est.cKaverage(2:end).^2))];
W.est.cTC=[0 -1./log(W.Mc.mu(2:end))]/(X.fSample0/X.downSampling); % real time units
W.est.cNC=-1./log(W.Mc.mu);                              % sampling time units
%% potential troublemakers and memory or processor hogs
if(do_estimates)
    try
        W.est.sRates=logm(W.est.A)*(X.fSample0/X.downSampling);
        W.est.sRates_error='none';
    catch me
        W.est.sRates=0*W.est.A;
        W.est.sRates_error=me;
    end
     W.est2.qt=qt;
     W.est2.qstLMax=qstLMax;
     W.est2.alpha=alpha;
     W.est2.beta=beta;
     W.est2.qst=qst;
     W.est2.lnqst=lnqst;
     [~,W.est2.zMaxP]=max(qt,[],2);
     W.est2.sMaxP=uint8(sMap(W.est2.zMaxP));
     W.est2.cMaxP=uint8(cMap(W.est2.zMaxP));
    
     % remnants of some tests to track down Viterbi weirdness: it looks like
     % to weak prioor is the problem, but I should really reimplement a
     % Viterbi algorithm that uses log(Q) and log(qst) instead, to avoid
     % problems with non-normalizeable distributions.
     % W.est2.zViterbi=uint8(sViterbi(Q,qst)); % Viterbi path
     % W.est2.zViterbi=uint8(VBviterbi(Q,qst)); % Viterbi path, faster algo     
     % [W.est2.zViterbi,s_lin,tw,t0]=VBlogViterbi(lnQ,lnqst,Q,qst);
     % if(~isempty(t0))
     %    clf, hold on
     %    plot(W.est2.zViterbi,'b')
     %   plot(s_lin+0.02,'m--')
     %    plot(W.est2.zMaxP+0.04,'k')
     %    plot(tw,s_lin(tw),'bx','markersize',10)
     %    keyboard
     % end
     %W.est2.zViterbi=uint8(VBlogViterbi(lnQ,lnqst,Q,qst));
     W.est2.zViterbi=uint8(VBviterbi_log(lnQ,lnqst)); % faster!     
     W.est2.sViterbi=uint8(sMap(W.est2.zViterbi));
     W.est2.cViterbi=uint8(cMap(W.est2.zViterbi));
end
%W.minimalStorage=@minimalStorage;
% ML: this function turns out to be very bulky to store, and was therefore
% removed from the structure. The internal function must still be present
% however, in order to be able to load old data. 
% Instead, added a few lines to make sure it is removed from every
% structure that is touched by this iteration again.
if(isfield(W,'minimalStorage'))
    W=rmfield(W,'minimalStorage');
end

    function WM=minimalStorage() % remove bulky fields
        WM=W;
        if(isfield(W,'est2'))
            WM=rmfield(W,'est2');
        end
        if(isfield(W,'minimalStorage'))
            W=rmfield(W,'minimalStorage');
        end
    end
end

