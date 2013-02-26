function [lns,s,S,t_warn,t0]=VBlogViterbi(lnQ,lnqst,Q,qst) 
% [lns,s,S,t_warn,t0]=VBlogViterbi(lnQ,lnqst,Q,qst)  
% most likely trajectory by the Viterbi algorithm, using log of transition
% matrix lnQ and emission likelihood lnqst in parallel with the linear 
% representation (which use Q, qst). 
%
% Thi file contains three different versions of the Viterbi algorithm, and
% was created for debugging purposes only (to check that they are indeed
% equivalent). lns,s and S are the three Viterbi paths, 
% t_warn/t0 are time points when some/all path probabilities are zero,
% probably due to numerical underflow in the algorithm Q,qst as opposed to
% lnQ=log(Q) and lnqst=log(qst). 
%
% M.L. 2011-12-21

[T,N]=size(lnqst);
% linear representation
pt=zeros(T,N);
maxPrevious=zeros(T,N);
pt(1,:)=qst(1,:)/sum(qst(1,:)); % initial probability
pp=zeros(1,N);
% log representation (in parallel)
lnpt=zeros(T,N);
lnMaxPrevious=zeros(T,N);
lnpt(1,:)=lnqst(1,:)-mean(lnqst(1,:)); % initial probability, not normalized
lnpp=zeros(1,N);
t_warn=[];
t0=[];
% algorithm version without large lnpt array
lnPP=zeros(1,N);
lnP0=zeros(1,N);
lnP1=zeros(1,N);
MaxPrev=zeros(T,N,'uint16'); % state variable
lnP1=lnqst(1,:)-mean(lnqst(1,:)); % initial probability, not normalized

for tV=2:T
    lnP0=lnP1;
    lnP1=zeros(1,N);
    for jV=1:N      % current target state
        for kV=1:N  % previous state
            lnPP(kV)=lnP0(kV)     +lnQ(kV,jV)+lnqst(tV,jV); % probability of most likely path that ends with kV -> jV            
            lnpp(kV)=lnpt(tV-1,kV)+lnQ(kV,jV)+lnqst(tV,jV); % probability of most likely path that ends with kV -> jV
            pp(kV)  =pt(tV-1,kV)*Q(kV,jV)*qst(tV,jV); % probability of most likely path that ends with kV -> jV
        end
        if(sum(pp==0)>0)
            %disp(['zero encounter: pp(' int2str(tV) ') = [' num2str(pp,3) '].'])
            %disp(['log repr.   : lnpp(' int2str(tV) ') = [' num2str(lnpp,3) '].'])
            if(isempty(t0) || t0(end)~=tV)
                t0(end+1)=tV;
            end
            if(sum(pp)==0)
                warning('linear repr. not normalizable')
                if(isempty(t_warn) || t_warn(end)~=tV)
                    t_warn(end+1)=tV;
                end
            disp('----------')
            end
            
        end
        % probability of previous state before ending up at jV.
        [lnP1(jV),         MaxPrev(tV,jV)]=max(lnPP);         
        [lnpt(tV,jV),lnMaxPrevious(tV,jV)]=max(lnpp); 
        [pt(tV,jV),maxPrevious(tV,jV)]=max(pp); 
    end
    % rescale, since we only need to keep track of relative probabilities.
    lnpt(tV,:)=lnpt(tV,:)-mean(lnpt(tV,:)); 
    pt(tV,:)=pt(tV,:)/sum(pt(tV,:));    
    lnP1=lnP1-mean(lnP1); % rescale to avoid numerical under- or overflow.

    % debug
    %c1=sum(abs(lnpp-lnPP));
    %c2=sum(abs(lnpt(tV,:)-lnP1));
    %c3=sum(abs(lnMaxPrevious(tV,:)-double(MaxPrev(tV,:))));    
    %if( c1>0 || c1>0 || c3>0 )
    %keyboard
    %end
end
lns=zeros(T,1);
[~,lns(T)]=max(lnpt(T,:));  
s=zeros(T,1);
[~,s(T)]=max(pt(T,:));  

S=zeros(T,1,'uint16');
[~,S(T)]=max(lnP1);  

for tV=T-1:-1:1
    S(tV)=MaxPrev(tV+1,S(tV+1));
    lns(tV)=lnMaxPrevious(tV+1,lns(tV+1));
    s(tV)=maxPrevious(tV+1,s(tV+1));
end
disp(sum(abs(lns-double(S))))
end

