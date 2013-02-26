function W1=VB7_GSconversion(W0,GtoS,StoG)
%% W1=VB7_GSconversion(W0,GtoS,StoG)
%
% function to move states between the 'genuine' and 'spurious' pools, while
% keeping the parameters as unchanged as possible. The idea is to move
% states back and forth between the two categories, to see which gives a
% better fit. Only the M, Mh, and prior fields are kept, the E-step and
% estimates fields are removed.
% 
% W0    : parent VB7 structure
% GtoS,StoG: indices of genuine/spurious states that are to be changed. For
% now, only one state can be changed at a time.

%% change-log
% M.L. 2010-12-27   : started on VB4 converter
% M.L. 2011-01-05   : started on VB5 converter
% M.L. 2011-02-03   : translated to VB7 converter
% M.L. 2011-02-21   : added handling of multiple conversions (by converting
%                     one state at a time)
% M.L. 2011-03-15   : caught handling of no conversion, and of converting
%                     all genuine states to spurious
% M.L. 2012-02-22   : changed how transition counts are converted in
%                     M.wA during genuine->spurious conversion 
% M.L. 2012-03-15   : switched to VB7_scmap_prior instead
%% action
if(~exist('StoG'))
    StoG=[];
end

if( isempty(GtoS) &&  isempty(StoG))
    W1=W0;
else
    if(length(GtoS)==1 && isempty(StoG) && W0.N>1 ) % convert a genuine state to spurious
        clear StoG;
        N=W0.N-1;
        Nc=W0.Nc+1;
        %W1=VB7_priorParent(W0,N,Nc); % new prior parameters
        W1=VB7_scmap_prior(W0,N,Nc);
        if(~isfield(W1,'PM'))
            warning('VB7_convertStates: incoming state did not come with correct prior parameters. Using defaults.')
            W1=VB7_priorParent(W1,N,Nc);
            W1=VB7_scmap_prior(W1,N,Nc);
        end
        % insert spurious state
        W1.Mc.wPi=[W0.Mc.wPi W0.M.wPi(GtoS)];%%
        W1.Mc.n=  [W0.Mc.n   W0.M.n(GtoS)  ];
        W1.Mc.c=  [W0.Mc.c   W0.M.c(GtoS)  ];
        W1.Mc.mu= [W0.Mc.mu  W0.M.mu(GtoS) ];
        W1.Mc.v=  [W0.Mc.v   W0.M.v(GtoS)  ];
        % remove genuine state
        sk=[1:GtoS-1 GtoS+1:W0.N];
        W1.M.wPi=W0.M.wPi(sk);
        W1.M.n  =  W0.M.n(sk);
        W1.M.c  =  W0.M.c(sk);
        W1.M.mu = W0.M.mu(sk);
        W1.M.v  =  W0.M.v(sk);
        
        % transfer observed transitions
        % observed transitions before conversion
        wA=W0.M.wA-W0.PM.wA;     
        wAc=W0.Mc.wA-W0.PMc.wA;
        
        % move transitions from G->G, S1->S1 to G,S1 -> Snew
        %W1.Mc.wA=[wAc(:,1)-wA(:,GtoS) wAc(:,2:end) wA(:,GtoS)]; % can lead
        %to negative transition counts in wAc(:,1)
        W1.Mc.wA=[sum(wA,2)-wA(:,GtoS) wAc(:,2:end) wA(:,GtoS)]; % hopefully better
        W1.Mc.wA=W1.Mc.wA(sk,:);     % remove old genuine state
        W1.Mc.wA=W1.Mc.wA+W1.PMc.wA; % add prior back in
        
        wR=W0.Mc.wR-W0.PMc.wR; % observed unsticking before conversion
        W1.Mc.wR=[wR wR(:,1)/W0.N; wAc(GtoS,:) wA(GtoS,GtoS)];
        W1.Mc.wR=W1.Mc.wR+W1.PMc.wR;
        
        % try to compensate for transitions via the removed state
        toR=wA(sk,GtoS);
        frR=wA(GtoS,sk);
        W1.M.wA=wA(sk,sk)+W1.PM.wA+0*(toR*ones(size(frR))+ones(size(toR))*frR);
    elseif(length(GtoS)>1 && isempty(StoG) && W0.N>1) % convert several genuine states to spurious
        clear StoG;
        GtoS=-sort(-GtoS); % sort in decreasing order
        W1=W0;
        for k=1:length(GtoS)
            W1=VB7_GSconversion(W1,GtoS(k));
        end
        
    elseif(isempty(GtoS)   && length(StoG)==1 && StoG~=1 ) % convert a spurious state to genuine
        clear GtoS;
        warning('VB7_convertToDirtstate: conversion to genuine state not tested!');
        N=W0.N+1;
        Nc=W0.Nc-1;
        %W1=VB7_priorParent(W0,N,Nc); % new prior parameters
        W1=VB7_scmap_prior(W0,N,Nc);
        if(~isfield(W1,'PM'))
            warning('VB7_convertStates: incoming state did not come with correct prior parameters. Using defaults.')
            W1=VB7_priorParent(W1,N,Nc);
            W1=VB7_scmap_prior(W1,N,Nc);
        end
        % insert new genuine state
        W1.M.wPi=[W0.M.wPi W0.Mc.wPi(StoG)];
        W1.M.n=  [W0.M.n   W0.Mc.n(StoG)  ];
        W1.M.c=  [W0.M.c   W0.Mc.c(StoG)  ];
        W1.M.mu= [W0.M.mu  W0.Mc.mu(StoG) ];
        W1.M.v=  [W0.M.v   W0.Mc.v(StoG)  ];
        % remove spurious state
        sk=[1:StoG-1 StoG+1:W0.Nc]; % spurious states to keep
        W1.Mc.wPi=W0.Mc.wPi(sk);
        W1.Mc.n = W0.Mc.n(sk);
        W1.Mc.c = W0.Mc.c(sk);
        W1.Mc.mu= W0.Mc.mu(sk);
        W1.Mc.v = W0.Mc.v(sk);
        
        % move transition points around
        wA =W0.M.wA -W0.PM.wA;
        wAc=W0.Mc.wA-W0.PMc.wA;
        wR =W0.Mc.wR-W0.PMc.wR;
        
        W1.Mc.wA=[wAc(:,sk); wR(StoG,sk)];
        W1.Mc.wA=W1.Mc.wA+W1.PMc.wA;
        
        W1.Mc.wR=wR(sk,sk)+W1.Mc.wR;
        
        W1.M.wA=[wA                                 wAc(:,StoG);
            wR(StoG,1)/(W1.N-1)*ones(1,W1.N-1) wR(StoG,StoG)];
        W1.M.wA=W1.M.wA+W1.PM.wA;
    else
        disp(['VB7_convertStates : GtoS  = [ ' int2str(GtoS(1:end)) ' ]'])
        disp(['VB7_convertStates : StoG = [ ' int2str(StoG(1:end)) ' ]'])
        error('VB7_convertStates : spurious state 1 canot be converted, and at least one genuine state must be left')
    end
end