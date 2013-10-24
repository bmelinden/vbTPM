%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_VBEMiter_nomex.m, VBEM iteration without mex files, in the vbTPM package
% =========================================================================
% 
% Copyright (C) 2013 Martin Lindén
% 
% E-mail: bmelinden@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or any later
% version.   
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% % Public License for more details.
% 
% Additional permission under GNU GPL version 3 section 7
% 
% If you modify this Program, or any covered work, by linking or combining it
% with Matlab or any Matlab toolbox, the licensors of this Program grant you 
% additional permission to convey the resulting work.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%% start of actual code
function [Wfun,init_parameters]=VB7_initialGuess_KBregion(x,tau0,tauS,wKBT,KBGrange,KBbins,KBSrange,tSigma,fSample,downSample,fig)
% [Wfun,init_parameters]=VB7_initialGuess_KBregion(x,tau0,tauS,wKBT,KBGrange,KBbins,KBSrange,tSigma,fSample,downSample,fig)
%
% returns function handle Wfun(W) that generates initial conditions for a
% model W with existing W.N size field, and prior fields. 
% The emission parameters are partly (2/3) placed within the box KBGrange
% and the spplied trajectory, and partly (1/3) distributed evenly in the
% box KBSrange. Spurious states are evenly distributed in KBSrange.
% This initial guess also adds 50% sticking states near K=[0.5 1] and 'high' B-values;
%
% parameters:
% x         : trajectory for which to make initial huess
% tSigma    : width of Gaussian filter to estimate K,B values 
%             (defailt=4 s)
% KBGrange =[Kmin Kmax Bmin Bmax] region in which to search for initial
%             values of genuine states. Default [0 1 0 2e-4]
% KBSrange =[Kmin Kmax Bmin Bmax] region in which to place initial guess
%             for spurious states. Default [0 1 0 1e-3]
% KBbins    : number of bins along each axis, default 100. Also accepts a
%             two-vector [Kbins Bbins]
% 
% wKBT      : strength (=number of counts) of initial guess. Default 10000.
% tau0      : mean dwell time [s] in each state. Default 2s.
% tauS      : mean time [s] between spurious events.  Default 10 s.
% fSample   : sample freqyency [Hz] of x
% downSample: planned downSample factor for analysis (rescales dwell times
%             to pseudo counts).
% fig       : figure number or axis handle. Optional. If given, a K,B
%             intensity plot and initial guess points will be generated.
%             An empty of invalid value will open a new window.
%

%% change-log
% M.L. 2010-10-19 : started, from VBE1_initialGuess_KBlist.m
% M.L. 2010-10-27 : fixed bug that sometimes gave |W.M.mu|>1 and W.M.b<0
%                   (which leads to errors)
% M.L. 2010-11-03 : added explicit downsample parameter to avoid confusion.
% M.L. 2011-02-02 : translated to VB7
% M.L. 2011-02-03 : different regions for genuine and spurious state
%                   placements, KBGrange, KBSrange.

%% default parameter values 
if(~exist('tau0','var')   ||   isempty(tau0));      tau0=2;  end
if(~exist('tauS','var')   ||   isempty(tauS));      tauS=10;  end
if(~exist('wKBT','var')   ||   isempty(wKBT));      wKBT=10000;  end
if(~exist('KBGrange','var')||   isempty(KBGrange));KBGrange=[0 1 0 2e-4];  end
if(~exist('KBSrange','var')||   isempty(KBSrange));KBSrange=[0 1 0 1e-3];  end
if(~exist('KBbins','var') ||   isempty(KBbins));  KBbins=[500 500]; end
if(length(KBbins)==1); KBbins=KBbins*[1 1];end
if(~exist('fSample','var') ||  isempty(fSample));   fSample=30;  end
if(~exist('downSample','var') ||  isempty(downSample));   downSample=1;  end
if(~exist('tSigma','var') ||   isempty(tSigma));    tSigma=2;  end

% checks
if(KBSrange(4)<=KBSrange(3) || KBSrange(2)<=KBSrange(1))
   error('inconsistent KBSrange. Intervals must be >0') 
end
if(KBGrange(4)<=KBGrange(3) || KBGrange(2)<=KBGrange(1))
   error('inconsistent KBGrange. Intervals must be >0') 
end

%% compute K,B histogram
[~,Ktrj,Btrj]=RMSKBgaussfilter(x,tSigma,fSample);


KrangeFraction=sum( (Ktrj>=KBGrange(1)).*(Ktrj<=KBGrange(2)))/length(Ktrj);
BrangeFraction=sum( (Btrj>=KBGrange(3)).*(Btrj<=KBGrange(4)))/length(Btrj);
disp(['VB7_initialGuess_KBregion emission parameter fractions:' ...
    num2str([KrangeFraction BrangeFraction],3)])

[~,Krange,Brange,Z,dK,dB]=XYhistogram(Ktrj,Btrj,KBbins,KBGrange);
[Kic,Bic]=THRregion(0.9,Krange,Brange,Z);
if(length(Kic)<30 || length(Bic) <30)
   warning('VB7_initialGuess_KBregion: less than 30 initial guess parameters detected. Consider changing the KBGrange or KBbins parameters.')
end
% remove degeneracies
dKic=linspace(dK/length(Kic),dK,length(Kic))';
dBic=linspace(dB/length(Kic),dB,length(Kic))';
Kic=Kic+dKic(randperm(length(Kic)));
Bic=Bic+dBic(randperm(length(Kic)));

if(exist('fig','var')) % draw the KB region
    try
        figure(fig);
    catch
        try
            axes(fig);
        catch
            figure
        end
    end
    XYhistogram(Ktrj,Btrj,KBbins,KBGrange,fig);
    hold on
    pp=plot(Kic,Bic,'wo');
    set(pp,'markersize',3,'markerface','w')
    disp(['number of points in the initial guess pool: ' int2str(length(Kic))])
    pause(0.1)
end
%% create parameter object
init_parameters.tau0=tau0;
init_parameters.tauS=tauS;
init_parameters.wKBT=wKBT;
init_parameters.KBGrange=KBGrange;
init_parameters.KBSrange=KBSrange;
init_parameters.KBbins=KBbins;
init_parameters.tSigma=tSigma;
init_parameters.fSample=fSample;
init_parameters.downSample=downSample;
init_parameters.dK=dK;
init_parameters.dB=dB;
init_parameters.Kic=Kic;
init_parameters.Bic=Bic;

Wfun=@Winit;
%% actual initial condition function
    function W=Winit(W)
        N=W.N; % size of this model
        Nc=W.Nc;
        
        % initial state
        W.M.wPi=10*ones(1,N);
        W.Mc.wPi=10*ones(1,Nc);
        
        %% dwell times and transition rates        
        %% guess transition matrices
        tss0=tau0; % s
        tcc0=tau0;  % s
        tsc0=tauS; % s
        
        % generate transition rates
        tss=tss0.*exp(0.2*randn(N,1));
        tsc=tsc0.*exp(0.2*randn(N,1));
        tcc=tcc0.*exp(0.2*randn(Nc,1));
        
        Pss=1-1./tss/fSample*downSample;
        Psc=1-1./tsc/fSample*downSample;
        Pcc=1-1./tcc/fSample*downSample;
        % guard against impossible parameters
        Pss(Pss<=0)=0.1;
        Psc(Psc<=0)=0.1;
        Pcc(Pcc<=0)=0.1;
        
        Ass=2+rand(N,N);
        Asc=2+rand(N,Nc);
        Acc=2+rand(Nc,Nc);
        
        Ass=Ass-diag(diag(Ass));
        Asc(:,1)=0;
        for k=1:N
            Ass(k,:)=(1-Pss(k))*Ass(k,:)/sum(Ass(k,:));
            Ass(k,k)=Pss(k);
            Asc(k,:)=[Psc(k) (1-Psc(k))*Asc(k,2:end)/sum(Asc(k,2:end))];
        end
        Acc=Acc-diag(diag(Acc));
        Acc(1,:)=0;
        for k=2:Nc
            Acc(k,:)=(1-Pcc(k))*Acc(k,:)/sum(Acc(k,:));
            Acc(k,k)=Pcc(k);
        end
        %% adjust prior strengths and deal with degenerate cases
        if(N>1)
            W.M.wA=Ass*wKBT;
            if(Nc>1)
                W.Mc.wA=Asc*wKBT;
                W.Mc.wR=Acc*wKBT;
            else
                W.Mc.wA=wKBT*ones(N,1);
                W.Mc.wR=0;
            end
        else
            W.M.wA=wKBT;
            if(Nc>1)
                W.Mc.wA=Asc*wKBT;
                W.Mc.wR=Acc*wKBT;                
            else
                W.Mc.wA=wKBT;
                W.Mc.wR=0;
            end
        end
        
        % initial guess for noise parameters
        NL=ceil(2*N/3);
        if(length(Kic)>NL)
            n0=randsample(length(Kic),NL,0); % sample without replacement
        else
           error('VB7_initialGuess_KBregion: too few emission parameter points found. Adjusting the range of initial emission parameters for (KBGrange) might help.')
        end
        K0=Kic(n0);
        B0=Bic(n0);
        
        % candidates for spurious states
        NS=N-NL;
        K0S=KBSrange(1)+rand(NS,1)*diff(KBSrange([1 2]));
        B0S=KBSrange(3)+(1e-5+rand(NS,1))*diff(KBSrange([3 4]));

        K0=[K0;K0S]';
        B0=[B0;B0S]';
        
        K0=min(K0,0.98); % make sure we start inside the stability limits
        K0=max(K0,-0.98);%
        
        W.M.n=wKBT*ones(1,N);
        W.M.c=W.M.n./B0;
        W.M.mu=K0;
        W.M.v=wKBT*ones(1,N);
        
        % actual spurious states
        K0c=[0 KBSrange(1)+rand(1,Nc-1)*diff(KBSrange([1 2]))];
        B0c=[0 KBSrange(3)+(1e-5+rand(1,Nc-1))*diff(KBSrange([3 4]))];

        K0c=min(K0c,0.98); % make sure we start inside the stability limits
        K0c=max(K0c,-0.98);%
        
        W.Mc.n=wKBT*[0 ones(1,Nc-1)];
        W.Mc.c=[0 W.Mc.n(2:end)./B0c(2:end)];
        W.Mc.mu=K0c;
        W.Mc.v=wKBT*[0 ones(1,Nc-1)];
        
        
    end
end

