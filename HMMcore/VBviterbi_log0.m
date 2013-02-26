function S=VBlogViterbi(lnQ,lnqst) 
% S=VBlogViterbi(lnQ,lnqst) 
% most likely trajectory by the Viterbi algorithm, using log of transition
% matrix lnQ and emission likelihood lnqst. This is the matlab version of
% VBviterbi_log.c, created for debugging purposes, but very much slower.

[T,N]=size(lnqst);
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
        end
        % probability of previous state before ending up at jV.
        [lnP1(jV),         MaxPrev(tV,jV)]=max(lnPP);         
    end
    lnP1=lnP1-mean(lnP1); % rescale to avoid numerical under- or overflow.
end
S=zeros(T,1,'uint16');
[~,S(T)]=max(lnP1);  
for tV=T-1:-1:1
    S(tV)=MaxPrev(tV+1,S(tV+1));
end
end

