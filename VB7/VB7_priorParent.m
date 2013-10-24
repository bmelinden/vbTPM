% use 1: create an empty model object from scratch
% Wparent=VB7_priorParent(fSample,ds,fPi,tD,tA,tDc0,fB,B0,fBc,Bc0,K0,Kstd,Kc0,KcStd,KBscaling)
% creates a VB7 parent structure with prior distribution parameters only,
% from which new objects with prior parameter fields can be created.
%
% fSample : sampling frequency, used to convert dwell
%           times to pseudo counts (i.e., K,B values not affected). 
%           Default 30 Hz.
% ds      : downsampling factor (default 1 =  no downsampling).
% fPi     : strength of prior on starting state (default=5).
% tD,tA   : Pseudocounts in the transition matrices are parameterized by an
%           expected dwell time tD [s] and a total strength tA [s] per row
%           (using fSample to set the time scale). Default tD=1 s, tA=1 s.
% tDc0    : expected dwell time of the 'genuine' state (e.g., time between
%           spurious events). Default 10 s.
% fB,  B0 : strength and expected value of B prior for s (default: fB=1,B0=1e-4).
% fBc,Bc0 : strength and expected value of B prior for c (default:fBc=1,Bc0=1e-4).
% Kstd,K0 : total variance and expectation value of K prior for s (default
%           Kstd=0.3, K0=0.5).
% KcStd,Kc0: standard deviation and expectation value of K prior for c 
%           (default KcStd=1, Kc0=0.8). 
% KBscaling : determines how the prior strength of the emission
%             parameters scales with number of states.
%           1: total strength held constant
%           2: individual state strength held constant (default)
%
% use 2: create models with prior parameters and specified number of
% states, based on parameters in another model object (possibly onme
% created by use 1 above).
%
% W=VB7_priorParent(Wparent,N,Nc)
% creates a new VB7 object with N states and Nc dirt states, with prior
% parameters according to the options in Wparent.priorParameterOptions. The
% created object inherits the prior options, and can in turn be used to
% generate new onjects.

%% change-log
% M.L. 2011-01-05   : started translation to VB5, added a downSample parameter
% M.L. 2011-02-02   : started translation to VB7, removed zeros in Rc prior
% M.L. 2011-02-14   : made PM.c, PMc.c consistent with constant <B> value,
%                     to fix bug with inducing splits of unpopulated
%                     states.
% M.L. 2011-02-25   : new parameter KBscaling to make the emission priors
%                     indepent of model size. (I think this is sensible,
%                     since the number of states is primarily dependent on
%                     the number of sticking events in the trajectory).
% M.L. 2012-03-15   : default KBscaling=2, not 1.
% M.L. 2012-05-02   : changed how wPi scales with number of states (equal
%                     weight on both spurious and genuine states).

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_priorParent.m, prior creator for the vbTPM package
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
function W=VB7_priorParent(fSample,ds,fPi,tD,tA,tDc0,fB,B0,fBc,Bc0,K0,Kstd,Kc0,KcStd,KBscaling)

%% determine if the call is to create a new parent or a new VB7 object
% use 1: VB7_priorParent(fSample,ds,fPi,...)
% use 2: VB7_priorParent(Wparent,N ,Nc)
if(nargin==3 && isfield(fSample,'priorParameterOptions'))
    WP=fSample; % rename to reduce confusion...
    N=ds;
    Nc=fPi;
    W=Winit(WP,N,Nc);
else
%% handle parameters
if(~exist('tD','var') || isempty(tD)); tD=1; end
if(~exist('tA','var') || isempty(tA)); tA=5; end
if(~exist('tDc0','var') || isempty(tDc0)); tDc0=10; end
if(sum(size(tDc0))>2)
    warning('VB7_priorparent: ignoring all but first element of tDc0');
    tDc0=tDc0(1);
end
if(~exist('fPi','var') || isempty(fPi)); fPi=5; end
if(~exist('fB','var') || isempty(fB)); fB=1; end
if(~exist('B0','var') || isempty(B0)); B0=1e-4; end
if(~exist('fBc','var') || isempty(fBc)); fBc=1; end
if(~exist('Bc0','var') || isempty(Bc0)); Bc0=1e-4; end
if(~exist('Kstd','var') || isempty(Kstd)); Kstd=0.3; end
if(~exist('K0','var') || isempty(K0)); K0=0.5; end
if(~exist('KcStd','var') || isempty(KcStd)); KcStd=1; end
if(~exist('Kc0','var') || isempty(Kc0)); Kc0=0.8; end
if(~exist('fSample','var') || isempty(fSample) || ~isreal(fSample));fSample=30; end
if(~exist('ds','var') || isempty(ds)); ds=1; end   
if(~exist('KBscaling','var') || isempty(KBscaling)); KBscaling=2; end   
%% add prior parameter field
pargs.fPi=fPi;
pargs.tA=tA;pargs.tD=tD;
pargs.tDc0=tDc0;
pargs.fB=fB;pargs.B0=B0;
pargs.fBc=fBc;pargs.Bc0=Bc0;
pargs.K0=K0;pargs.Kstd=Kstd;
pargs.Kc0=Kc0;pargs.KcStd=KcStd;
pargs.fSample=fSample;
pargs.downSample=ds;
pargs.KBscaling=KBscaling;

W.priorParameterOptions=pargs;
end
%% code to initialize new VB7 object
    function W=Winit(WP,N,Nc)
     
        W.priorParameterOptions=WP.priorParameterOptions;
        W.N=N;
        W.Nc=Nc;
        
        opt=W.priorParameterOptions;
        % initial state
        wPi0=opt.fPi/(N+Nc-1);
        W.PM.wPi = ones(1,N)*wPi0;   % prior on initial state distribution
        W.PMc.wPi= ones(1,Nc)*wPi0;
        W.PMc.wPi(1)=sum(W.PM.wPi);
        % transition rates
        if(N>1)
            Q=-1/opt.tD*eye(N)+1/opt.tD/(N-1)*(ones(N,N)-eye(N)); % rate matrix [1/s]
            W.PM.wA = opt.tA*opt.fSample/opt.downSample*expm(Q/opt.fSample*opt.downSample);
        else
            W.PM.wA=opt.tA*opt.fSample/opt.downSample;
        end
        
         if(Nc>0)
            pC=exp(-1/opt.tDc0/(opt.fSample/opt.downSample)); % probability to stay genuine
            W.PMc.wA=opt.tA*opt.fSample/opt.downSample*[pC*ones(N,1) (1-pC)/(Nc-1)*ones(N,Nc-1)];            
            pS=exp(-1/opt.tD/(opt.fSample/opt.downSample)); % probability to stay stuck
            W.PMc.wR=opt.tA*(opt.fSample/opt.downSample)*...
                [zeros(1,Nc); (1-pS)/(Nc-1)*([ones(Nc-1,1) 1-eye(Nc-1,Nc-1)])+[zeros(Nc-1,1) pS*eye(Nc-1,Nc-1)]];
        else
            W.PMc.wA=opt.tA*(opt.fSample/opt.downSample)*ones(Nc,1);
            W.PMc.wR=0;
         end 
        % emission model
        if(isfield(opt,'KBscaling') && opt.KBscaling==2)
            W.PM.n  = opt.fB*ones(1,N);         % prior on dispersion coeff. Bj
            W.PM.c  = (W.PM.n+0.5)/opt.B0;
            W.PM.v  = 1/2/opt.Kstd^2*ones(1,N); % prior for decay coeff. Kj
            W.PM.mu = opt.K0.*ones(1,N);
            % c=1 does not have emission parameters:
            W.PMc.n  = [0 opt.fBc*ones(1,Nc-1)];       % prior on dispersion coeff. Bj
            W.PMc.c  = [0 (W.PMc.n(2:end)+0.5)/opt.Bc0];
            W.PMc.v  = [0 1/2/opt.KcStd^2*ones(1,Nc-1)];% prior for decay coeff. Kj
            W.PMc.mu = [0 opt.Kc0.*ones(1,Nc-1)];
        else
            W.PM.n  = opt.fB/N*ones(1,N);         % prior on dispersion coeff. Bj
            %W.PM.c  = opt.fB/N/opt.B0.*ones(1,N);
            W.PM.c  = (W.PM.n+0.5)/opt.B0;
            W.PM.v  = 1/2/opt.Kstd^2*ones(1,N)/N; % prior for decay coeff. Kj
            W.PM.mu = opt.K0.*ones(1,N);
            % c=1 does not have emission parameters:
            W.PMc.n  = [0 opt.fBc/(Nc-1)*ones(1,Nc-1)];       % prior on dispersion coeff. Bj
            %W.PMc.c  = [0 opt.fBc/(Nc-1)/opt.Bc0.*ones(1,Nc-1)];
            W.PMc.c  = [0 (W.PMc.n(2:end)+0.5)/opt.Bc0];
            W.PMc.v  = [0 1/2/opt.KcStd^2/(Nc-1)*ones(1,Nc-1)];% prior for decay coeff. Kj
            W.PMc.mu = [0 opt.Kc0.*ones(1,Nc-1)];
        end
    end
end
