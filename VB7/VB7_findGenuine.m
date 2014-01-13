% [gCal,sCal,gTrj,sTrj]=VB7_findGenuine(Wcal,Wtrj,RMS0,fig)
%
% classifies the states of converged VB7 structures Wtrj (looping trace)
% and Wcal (calibration trace) into genuine and spurious, based on the
% following rules:
% 1) The calibration has one genuine state, which is the one with highest
%    total occupation time. 
% 2) The looping state closest to the calibration genuine state is the
%    genuine 'unlooped' state.
% 3) Genuine states satisfy K <= Ku, B>= Bu, where Ku,Ku are the
%    emission parameters of the unlooped state. (Since the unlooped state
%    corresponds to the longest tether, and K and B decreases and increases
%    with decreasing tether length respectively). 
% 4) Genuine states have a minimum RMS value RMS0< 1/B/(1-K^2). (Default is
%    RMS0=80 nm). 
%
% fig: visalize the classification in figure or figure handle fig
% (optional). fig=[] opens a new figure.
%
% output: gTrj,gCal are the geunie states, and sTrj,sCal are the spurious
% states of Wtrj and Wcal respectively.
% if Wtrj is not given, only the calibration structure is analyzed.

%% change-=log
% M.L. 2011-02-16   : started
% M.L. 2011-08-17   : added the option to visialize the classification
% M.L. 2011-09-13   : added the requirement K>0 for genuine states

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_findGenuine.m, automated classification of genuine/spurious states, 
% in the vbTPM package
% =========================================================================
% 
% Copyright (C) 2014 Martin Lindén
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
function [gCal,sCal,gTrj,sTrj]=VB7_findGenuine(Wcal,Wtrj,RMS0,fig)

%% default parameters3
if(~exist('RMS0','var') || isempty(RMS0)); RMS0=80; end
%% selection criteria
% find most occupied state in calibration trace
[~,gCal]=max(Wcal.est.sAverage); % genuine state
Kcal=Wcal.est.sKaverage(gCal);
Bcal=Wcal.est.sBaverage(gCal);
dKcal=Wcal.est.sKstd(gCal);
dBcal=Wcal.est.sBstd(gCal);

sCal=setdiff(1:Wcal.N,gCal); % spurious calbration states

if(exist('Wtrj','var') && ~isempty(Wtrj))
    K=Wtrj.est.sKaverage;
    B=Wtrj.est.sBaverage;
    dK=Wtrj.est.sKstd;
    dB=Wtrj.est.sBstd;
    
    % 1) identify unlooped state of looping trace, as 'closest' to the
    % calibration genuine state
    [~,sU]=min((K/Kcal-1).^2+(B/Bcal-1).^2);
    Ku=K(sU);
    Bu=B(sU);
    dKu=dK(sU);
    dBu=dB(sU);
    
    % lower right corner of genuine parameter area
    Kmax=max(Ku,Kcal+3*dKcal);%max(Ku+3*dKu,Kcal+3*dKcal);
    Bmin=min(Bu,Bcal-3*dBcal);%min(Bu-3*dBu,Bcal-3*dBcal);    
    
    % 2) find states on the right side of the calibration trace, and above
    % the RMS threshold.
    gTrj2=find( (K<=Kmax).*(K>0).*(B>=Bmin).*(Wtrj.est.sRMS>RMS0));
    gTrj=gTrj2;
    
    if(0) % debug: compare with old criterion
        % apply rule 3 and 4
        gTrj1=find( (0<=(Ku-K)).*(B-Bu>=0).*(Wtrj.est.sRMS>RMS0));
        
        
        gTrj12=union(gTrj1,gTrj2);
        if(length(gTrj1)==length(gTrj12) && length(gTrj2)==length(gTrj12))
            gTrj=gTrj1;
        else
            %% plot rule 2
            figure(10)
            clf
            hold on
            errorbarxyML(K,B,dK,dB,[],[],'k.','k')
            hold on
            errorbarxyML(Kcal,Bcal,dKcal,dBcal,[],[],'go','g')
            hold on
            pp=plot(Kcal,Bcal,'ko');set(pp,'markersize',8,'markerface','g')
            plot([ 0 Kmax Kmax],[Bmin Bmin max(B)],'--r')
            plot(Kmax,Bmin,'r*')
            
            %% plot difference
            pp=plot(K(gTrj1),B(gTrj1),'ro');set(pp,'markersize',10,'linew',2)
            pp=plot(K(gTrj2),B(gTrj2),'bx');set(pp,'markersize',10,'linew',2)
            pp=plot(Ku,Bu,'kd');set(pp,'markersize',10,'linew',2)
            
            
            keyboard
        end
        
        
        %% plot rule 2
        figure(10)
        clf
        pp=plot(Kcal,Bcal,'ko');set(pp,'markersize',8,'markerface','g')
        hold on
        errorbarxyML(K,B,dK,dB,[],[],'k.','k')
        hold on
        errorbarxyML(Kcal,Bcal,dKcal,dBcal,[],[],'go','g')
        hold on
        
        plot([ 0 Kmax Kmax],[Bmin Bmin max(B)],'--r')
        KK=linspace(0,0.99,100);
        plot(KK,1./RMS0^2./(1-KK.^2),'r--')
        plot(Kmax,Bmin,'r*')
        set(gca,'yscale','log')
        axis([0 1 3e-5 5e-4])
    end
    
    sTrj=setdiff(1:Wtrj.N,gTrj);
else
    gTrj=[];
    sTrj=[];
end
%% visualize results
if(exist('fig','var'))
    try
        figure(fig)
        hold on
    catch me
        if(strcmp(me.identifier,'MATLAB:Figure:IntegerHandle'))
            % then fig was empty, or an invalid figure handle; open a
            % new figure instead
            fig=figure();
            clf
        elseif(strcmp(me.identifier,'MATLAB:Figure:HandleInUse'))
            % the fig was a handle to an existing figure, which we use
            axes(fig)
        end
    end
    x=linspace(0,Kmax,100);
    maxB=1/RMS0^2./(1-Kmax.^2);
    set(gca,'xlim',[-0.05 1],'box','on',...
        'ylim',[0 1.1*max([1.2*maxB Wtrj.est.sBaverage Wcal.est.sBaverage])])
    
    pp=fill([x Kmax 0 0],[1/RMS0^2./(1-x.^2) Bmin Bmin 1/RMS0^2],'g');
    hold on
    set(pp,'facecolor',0.9*[1 1 1],'edgecolor','none','LineStyle','-')
    leg={'genuine range'};
    % 1) display calibration spurious states
    if(~isempty(sCal))
        spurious_cal=plot(Wcal.est.sKaverage(sCal),Wcal.est.sBaverage(sCal),'+');
        set(spurious_cal,'markeredge',0*[1 1 1])
        leg{end+1}='cal. spurious';
    end
    % 2) display trajectory spurious states
    if(~isempty(sTrj))        
        pp=plot(Wtrj.est.sKaverage(sTrj),Wtrj.est.sBaverage(sTrj),'x');
        set(pp,'markeredge',0*[1 1 1])
        leg{end+1}='trj. spurious';
    end
    if(~isempty(gCal))
        pp=plot(Wcal.est.sKaverage(gCal),Wcal.est.sBaverage(gCal),'ko');
        set(pp,'markeredge','k','markerface','r')
        leg{end+1}='cal. genuine';
    end
    if(~isempty(gTrj))
        pp=plot(Wtrj.est.sKaverage(gTrj),Wtrj.est.sBaverage(gTrj),'kd');
        set(pp,'markeredge','k','markerface','b')
        leg{end+1}='trj. genuine';
    end
    legend(leg,2)
    xlabel('K')
    ylabel('B [nm^{-2}]')
    %axis tight
    %    foo=get(gca,'ylim');
    %maxB=2^(ceil(log2(max(foo(2),Bmin)))); % round upwards
        
end



