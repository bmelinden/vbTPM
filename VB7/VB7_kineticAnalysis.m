% [res,R]=VB7_kineticAnalysis(runinputfile,loadfile,doplots,resort)
% 
% Create factorial HMM models with spurious states in a class of its own.
% example
% [res,R]=VB7_kineticAnalysis('RI_ds3run_100pM_tA01','ML_ds3_0100pM.mat')
%

% runinputfile          VB7 runinput file
% loadfile              VB7 file with state classification (e.g., from
%                       the VB7_batch_postprocess GUI).
% doplots=true/false    plot some stuff
% resort=true/false     reorder thestates such that the spurious states
%                       come last.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_kineticAnalysis.m, automated generation of factorial models based on 
% state classifications, part of the vbTPM package
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
function [res,R]=VB7_kineticAnalysis(runinputfile,loadfile,doplots,resort)

if(~exist('doplots','var') || isempty(doplots)); doplots=true;end
if(~exist('resort','var') || isempty(resort)); resort=false;end

%% load analyzed models and prune fields
savefile=[loadfile '_VB7kin.mat'];
if(exist(savefile) && ~resort)
    R=load(savefile);
elseif(exist(loadfile))
    R=load(loadfile);
else
    error(['files not found: ' loadfile ', ' savefile])
end

% check if sorting and state conversion is already done
if(~isfield(R,'isSorted') || R.isSorted==false)
    for k=1:length(R.trj)
        for b=1:length(R.trj{k})
            if(R.includetrj{k}(b))
                tic
                %% sort models, extract Viterbi paths
                W=R.trj{k}{b}.Wtrj;
                s=R.trjstates{k}{b};                
                [W,s,~]=sortModel(W,s);
                
                % normal model
                dat=VB7_getTrjData(runinputfile,k,b);
                dat.dat=VB7_preprocess(dat.x,dat.downSample,dat.fSample);               
                R.RMS{k}{b}=dat.RMS;
                W=VB7iterator(W,dat.dat);
                R.trjstates{k}{b}=s;
                R.trj{k}{b}.Wtrj=W;
                W=VB7_VBEMiter(W,dat.dat,'estimate');

                spr=find(R.trjstates{k}{b}==1); % spurious states
                [sVit,tdVit]=getDwellTRJ(W.est2.sViterbi,spr); % viterbi path with spurious states removed
                R.trj_vit_s{k}{b}=sVit;
                R.trj_vit_td{k}{b}=tdVit;
                [sML,tdML]=getDwellTRJ(W.est2.sMaxP);
                R.trj_sMaxP_s{k}{b}=sML;
                R.trj_sMaxP_td{k}{b}=tdML;
                %h.viterbi_t{k}{b}=t0Vit;
                % plot: 
                % stairs(cumsum([0 R.trj_vit_td{k}{b}]),R.trj_vit_s{k}{b}([1:end end]),'r-')
                
                % factorial model                
                srm=find(s==1);
                Wgs=VB7_GSconversion(W,srm);
                Wgs=VB7iterator(Wgs,dat.dat);
                R.trj{k}{b}.Wgs=Wgs;

                Wgs=VB7_VBEMiter(Wgs,dat.dat,'estimate');
                
                [sVit,tdVit]=getDwellTRJ(Wgs.est2.sViterbi);
                R.gs_vit_s{k}{b}=sVit;
                R.gs_vit_td{k}{b}=tdVit;
                [sML,tdML]=getDwellTRJ(Wgs.est2.sMaxP);
                R.gs_sMaxP_s{k}{b}=sML;
                R.gs_sMaxP_td{k}{b}=tdML;
                
                disp(['-----' int2str([k b]) ' --- ' num2str(toc) ' s.'])
            end
        end
    end
    R.isSorted=true;
    save(savefile,'-struct','R');
    disp(['saved to ' savefile])
else
    disp([loadfile ' already sorted.'])
end
%% construct aggregated models and histograms
pstates=2:length(R.statelabels);         % states to use for property histograms
res.pstates=pstates;
%res.pstatelabels={R.statelabels{pstates}};
res.pstatelabels=R.statelabels(pstates);

histP  =[];
histPgs=[];
histtD  =[];
histtDgs=[];
tDwell  =cell(1,0);
tDwellgs=cell(0);

n=0;
for k=1:length(R.trj)
    for b=1:length(R.trj{k})
        if(R.includetrj{k}(b))
            W  =R.trj{k}{b}.Wtrj;
            Wgs=R.trj{k}{b}.Wgs;
            n=n+1;
                        
            for m=1:length(pstates)
                % find state index
                ind=find(R.trjstates{k}{b}==pstates(m));
                
                if(~isempty(ind))
                    % fill out probability plots.
                    histP(n,m)=sum(W.est.Paverage(ind));
                    histPgs(n,m)=sum(Wgs.est.Paverage(ind));
                    
                    % mean life time, weighted average
                    histtD(n,m)  =sum(  W.est.Paverage(ind).*W.est.tD(ind)')/histP(n,m);
                    histtDgs(n,m)=sum(Wgs.est.Paverage(ind).*Wgs.est.tD(ind)')/histPgs(n,m);
                    
                    % build dwell time hiostograms
                    % not finished yet!
                    %for j=ind
                    %    s =find(R.trj_vit_td{k}{b}==j);
                    %    tD=R.trj_vit_td{k}{b}(s); % dwell times in this state
                    %    tDwell{m}=[tDwell{m} tD]  % collect dwell time list for this state
                    %    
                    %    s =find(R.gs_vit_td{k}{b}==j);
                    %    tD=R.gs_vit_td{k}{b}(s); % dwell times in this state
                    %    tDwellgs{m}=[tDwellgs{m} tD]  % collect dwell time list for this state
                    %end
                end
            end            
            histP(n,:)  =histP(n,:)/sum(histP(n,:));
            histPgs(n,:)=histPgs(n,:)/sum(histPgs(n,:));
        end
    end
end

%% plot stuff and add to result structure
%% looping probability plots
res.histP  =histP;
res.histPgs=histPgs;
if(doplots) 
    
    % unlooped + non-loopers
    figure(1),    clf
    subplot(2,1,1)
    
    hist(sum(res.histP(:,[1 5:7]),2),0.025:0.05:1) % unlooped states
    %hist(sum(res.histP(:,[2:4]),2),0.025:0.05:1) % looped 
    xlabel('p(s=U,NL)')
    title('unlooped probability')
    subplot(2,1,2)
    hist(sum(res.histPgs(:,[1 5:7]),2),0.025:0.05:1) % unlooped states
    %hist(sum(res.histPgs(:,[2:4]),2),0.025:0.05:1) % looped 
    title('unlooped probability, gsHMM')
    xlabel('p(s=U,NL)')

    % looping probabilities 
    figure(2),clf,hold on
    mar={'-k.','-rd','-bo'};    
    for m=1:3
        [a,b]=hist(histP(:,1+m),0.025:0.05:1);
        plot(b,a,mar{m})
        leg{m}=R.statelabels{pstates(m+1)};
    end
    xlabel('state probability')
    legend(leg)
    clear leg;
    box on
    
end
%% mean dwell times
res.histtD=histtD;
res.histtDgs=histtDgs;
if(doplots)
    tMax=1.5*max(max(histtD));
    figure(3),clf
    col='krbgmkrbgmkrbgm';
    mar='o>s<d^po>s<d^po>s<d^p';    
    subplot(2,1,1)
    hold on
    for m=1:size(histtD,2)
        [a,b]=hist(histtD(:,m),linspace(0,tMax,50));
        plot(b,a,['-' mar(m) col(m)])
        leg{m}=R.statelabels{pstates(m)};
    end
    legend(leg)
    clear leg;
    box on
    title('mean dwell time histograms, HMM')
    xlabel('<\tau_j>')
    subplot(2,1,2)
    hold on
    for m=1:size(histtDgs,2)
        [a,b]=hist(histtDgs(:,m),linspace(0,tMax,50));
        plot(b,a,['-' mar(m) col(m)])
    end
    box on
    title('mean dwell time histograms, parallel HMM')
    xlabel('<\tau_j>')
end
%% dwell time distributions: Viterbi vs model
res.tDwell=tDwell;
res.tDwellgs=tDwellgs;

if(doplots)
    
end
end
%% misc functions
function [W,s,ind]=sortModel(W,s)
% reorder the states such that spurious states come last

if(W.Nc==1)
sm=2*max(s);
s(s==1)=sm;     % put spurious states last
[s,ind]=sort(s);
sf={'PM','M','E','est'};
for f=1:length(sf)
    fn=fieldnames(W.(sf{f}));
    for g=1:length(fn)
        V=W.(sf{f}).(fn{g});        
        [a,b]=size(V);
        ia=1;ib=1;
        if(a==W.N)
            ia=ind;
        end
        if(b==W.N)
            ib=ind;
        end
        V=V(ia,ib);
        W.(sf{f}).(fn{g})=V;        
    end
end
s(s==sm)=1;
end
end
