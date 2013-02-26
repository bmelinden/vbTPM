function W1=VB7_scmap_prior(W0,N,Nc)
% W1=VB7_scmap_prior(W0,N,Nc)
% generate prior parameters for a mixed model with N genuine and Nc-1
% spurious states, based on the prior parameters for a genuine state-only model
% with N+Nc-1 states.
%
% mixed model: 
% p(st,ct|st-1,ct-1)=A(st-1,st)*(delta(ct-1,1)*Ac(st-1,ct) 
%       + (1-delta(ct-1,1))*Rc(ct-1,ct)

w=VB7_priorParent(W0,N+Nc-1,1); % construct genuine-state only prior

W1.priorParameterOptions=w.priorParameterOptions;
W1.N=N;
W1.Nc=Nc;

% initial state
W1.PM.wPi=w.PM.wPi(1:N);
W1.PMc.wPi=[sum(w.PM.wPi(1:N)) w.PM.wPi(N+1:end)];

% transitions
W1.PM.wA=w.PM.wA(1:N,1:N); % genuine -> genuine pseudocounts
W1.PMc.wA=[sum(w.PM.wA(1:N,1:N),2) w.PM.wA(1:N,N+1:end)];  % gen. -> [gen. spurious]
W1.PMc.wR= [0 zeros(1,Nc-1);
    sum(w.PM.wA(N+1:end,1:N),2) w.PM.wA(N+1:end,N+1:end)]; % spurios -> [gen. spurious]

% emission variables
W1.PM.n =w.PM.n(1:N);
W1.PMc.n=[0 w.PM.n(N+1:end)];
W1.PM.c =w.PM.c(1:N);
W1.PMc.c=[0 w.PM.c(N+1:end)];
W1.PM.v =w.PM.v(1:N);
W1.PMc.v=[0 w.PM.v(N+1:end)];
W1.PM.mu =w.PM.mu(1:N);
W1.PMc.mu=[0 w.PM.mu(N+1:end)];

