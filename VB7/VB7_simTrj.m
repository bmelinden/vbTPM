function [s,c,x1,x2]=VB7_simTrj(A,K,B,Asc,Rc,Kc,Bc,xyc,T,Twu)
% function [s,c,x1,x2]=VB7_simTrj(A,K,B,Asc,Rc,Kc,Bc,xyc,T,Twu)
% A(i,j)=p(s_t=j|s_{t-1}=i), so the evolution of the hidden states is
% given by
% p(t)=p(t-1)*A, where p(t) is a row vector
% Twu is a warm-up time, i.e., Twu time points are simulated before
% starting on the actual trace.
%
% Asc is the transition probabilities for the dirt events:
% Asc(i,j)=p(c_t=j| c_{t-1}=1, s_{t-1}=i)
%
% Rc is the transition probabilities for ending dirt events:
% Rc(i,j)=p(c_t=j| c_{t-1}=i, i>1).
% The constraint i>1 means that the first row of Rc is ignored. In the
% VB4-VB6 algorithms, spurious states are assumed non-interconvertable,
% meaning that only the diagonal and the first column of Rc is nonzero. VB7
% can infer spurious state interconversion as well.
%
% Kc, Bc are emission parameters of the 'dirt' process c(t): 
% when c(t)==1, the emission is determined by s(t). Otherwise, c(t)
% determines the output. The first elements of Kc,Bc will therefore not be
% used in the simulation. xyc are a vector where each column is the
% tethering point for a dirt state (default zero).
%
% x1,x2 are noise generated from the same hidden state trajectory, but with
% different methods and random numbers.
% x1 : simply use the Kc, Bc, and xyc parameters.
% x2 : simulate off-center sticking points by having the sticking events
%      diffuse around the last non-stick position preeding the sticking
%      event (more realistic).

%% change log
% M.L. 2011-01-06 : started on a VB5 version, with two different ways to
%                   simuate sticking events.
% M.L. 2011-01-16 : started on a VB6 version, with one more way to simuate
%                   sticking events (specifying a mean for each state).
% M.L. 2011-02-02 : changed name to VB7 (nothing else changed)

%% parameter handling
if(~exist('xyc','var') || isempty(xyc)); xyc=zeros(2,length(Kc)); end
if(~exist('Twu','var') || isempty(Twu)); Twu=ceil(5*max(-1./log(K))); end
Twu=ceil(Twu);
T=T+Twu;

N=length(K);
Nc=length(Kc);
if(Nc==1)
    Rc=1;
    Asc=ones(N,1);
end

nuc=(xyc*diag(1-Kc))';
%% start of actual code
%% generate state trajectories
s=zeros(T,1,'uint8');
s(1)=1;
cumA=cumsum(A,2);
c=zeros(T,1,'uint8');
c(1)=1;
cumAsc=cumsum(Asc,2);
cumRc=cumsum(Rc,2);
for t=2:T
    ras=rand;
    rac=rand;
    s(t)=find(ras<cumA(s(t-1),:),1);
    if(c(t-1)==1)
        c(t)=find(rac<cumAsc(s(t-1),:),1);
    else
        c(t)=find(rac<cumRc(c(t-1),:),1);        
    end
end
% generate noise trajectory
x1=zeros(T,2);
x2=zeros(T,2);
bb=1./sqrt(2*B);
hbb=1./sqrt(2*Bc);
x1(1,:)=randn(1,2)*bb(1);
x2(1,:)=randn(1,2)*bb(1);
mu=[0 0];
for t=2:T
    ww=randn(1,2);
    if(c(t)==1)
        x1(t,:)=K(s(t))*x1(t-1,:)+bb(s(t))*ww;
        x2(t,:)=K(s(t))*x2(t-1,:)+bb(s(t))*ww;
    else
        x1(t,:)=Kc(c(t))*x1(t-1,:)+nuc(c(t),:)+hbb(c(t))*ww;
        if(c(t-1)==1) % set new sticking point at beginning of event
            mu=x2(t-1,:);
        end
        x2(t,:)=Kc(c(t))*x2(t-1,:)+(1-Kc(c(t)))*mu+hbb(c(t))*ww;
    end
end
%% generate true VBE output statistics
T=T-Twu;
s=s(Twu+1:end,:);
c=c(Twu+1:end,:);
x1=x1(Twu+1:end,:);
x2=x2(Twu+1:end,:);

