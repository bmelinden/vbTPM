function s=VB7_inspectStates(W,X,s,KBax)
% s=VB7_inspectStates(W,X,s,KBax)
%
% a little script that plots coordinates X around each occurrence of each
% state. Only 'genuine' states (s) are considered.
%
% W     : converged VBE object
% X     : x,y coordinates to plot
% s     : state sequence to use (default: Viterbi sequence)      
% KBax  : axes in which to indicate current states (not cleared).
%         If omitted, a new figure is closed, a KB scatter plot is added,
%         and the figure is closed when finished. The trajectories are
%         always plotted in a new figure and closed afterwards.
% 

% M.L. 2011-08-22   : updated to conform with VB7 
% M.L. 2011.09-01   : removed excessive text output
% M.L. 2012-01-19   : handle downsampled trajectories, close new windows on
%                     quit.

%% non-adjustable parameters
RMS=gaussFilter(X,3); % filter width
flank=20;             % number of points before and after each 2-dwell
%% input parameters
dSample=W.priorParameterOptions.downSample;
fSample=W.priorParameterOptions.fSample;
if(~exist('s','var') || isempty(s))
    %disp('computing Viterbi path')
    if(~isfield(W,'est2'))
        dat=VB7_preprocess(X,dSample,fSample);
        W=VB7_VBEMiter(W,dat,'estimate');
        clear dat;
    end
    s=W.est2.sViterbi;
end
if(dSample>1) % upsample state path
    s0=s;
    s=zeros(1,1+(length(s(2:end))*dSample));
    s(1)=s0(1);
    for mm=1:dSample
        s(1+mm:dSample:end)=s0(2:end);
    end
    if(length(s)<length(RMS))
        if(length(RMS)-length(s)<dSample )
            s(end:length(RMS))=s(end);
        else
            disp('could not upsamle correctly, or inconsistent state sequence length')
            keyboard
        end
    end
    clear s0;
end
fig1=figure();
%% plot data and states
[STATE,DWELL,Dstart,Dend]=getDwellTRJ(s);
sWatch=STATE(1);
t=(1:size(X,1))/30;

figure(fig1),clf
f1=subplot(3,1,1);
hold on
fs=plot(t(s==sWatch),RMS(s==sWatch),'go');
plot(t,W.est.sRMS(s),'-r.')
plot(t,RMS,'k')
[m,n]=sort(W.est.sRMS);
%set(gca,'ygrid','on')%,'ytick',m,'yticklabel',int2str(n')
for mm=1:length(n)
    plot(t([1 end]),m(mm)*[1 1],'k:')
end
box on

ind=find(STATE==sWatch); 
    
f2=subplot(3,1,2);
hold on
plot(t,X(:,1),'-r');
plot(t,X(:,2)+0,'-b');
box on

f3=subplot(3,1,3);hold on
set(f3,'ylim',[0.5 W.N+0.5],'ytick',1:W.N,'ygrid','on')
plot(t,s,'-k')
fss=plot(t(s==sWatch),s(s==sWatch),'go');
box on

% plot emission parameters
if(exist('KBax','var') && ~isempty(KBax))
    axes(KBax)
    closeKBax=false;
    KBnextplot=get(KBax,'nextplot');
    hold on
else
    figure();
    KBax=gca;
    closeKBax=true;
    hold on
    plot(W.est.sKaverage,W.est.sBaverage,'k.')
end
pem=plot(W.est.sKaverage(sWatch),W.est.sBaverage(sWatch),'r+','markersize',10);
%xlabel('K')
%ylabel('B')

%% watch transition state events

keepWatching=true;
globalview=false;
k=1;
dt=t(2)-t(1);
while(keepWatching)
    if(isempty(ind))
        ind=[2];
        disp(['state ' int2str(sWatch) ' not occupied!'])
    end
    figure(fig1)
    % find time indices to watch
    ti=max(1,Dstart(ind(k))-flank):min(Dend(ind(k))+flank,length(s));
    
    i0=Dstart(ind(k));
    t0=t(i0);
    i1=Dend(ind(k));
    t1=t(i1);
        
    %set([f1 f2 f3],'xtick',sort(union([t0-dt/2 t1+dt/2],t(ti([1 end])))),...
    set([f1 f2 f3],'xtick',[t0-dt/2 t1+dt/2],'xgrid','on')
    set([f1 f2 f3],'xlim',t(ti([1 end])))
    set([f1 f2],'ylim',[-600 600])
    set([f1 f2 f3],'xticklabel',num2str(get(gca,'xtick')',5))
    if(globalview)
        set([f1 f2 f3],'xlim',t([Dstart(1) Dend(end)]))
    end
    axes(f1)
    box on
    ylabel('RMS [nm]')
    xlabel('t [s]')

    axes(f2)
    box on
    ylabel('coordinate [nm]')
    xlabel('t [s]')

    axes(f3)
    box on
    ylabel('state')
    xlabel('t [s]')

    if(ind(k)>1 && ind(k)<length(s) && k<length(STATE))
        axes(f1)
        if(ind(k)>1 && ind(k)<length(STATE))
            title(int2str(STATE(ind(k)+[-1 0 1])))
        elseif(ind(k)==1 && ind(k)<length(STATE))
            title(int2str(STATE(ind(k)+[0 1])))
        elseif(ind(k)>1 && ind(k)==length(STATE))
            title(int2str(STATE(ind(k)+[-1 0])))
        else
            title('cannot plot state sequence. Only one state?')
        end            
    end
    
    todo=menu('choose action',['next ' int2str(sWatch)],['previous ' int2str(sWatch)],...
              'next dwell','previous dwell',...
              'state +1','state -1',...
              'flank +1','flank -1','global/local','quit');
    switch todo
        case 1
            k=min(k+1,length(ind));
        case 2
            k=max(1,k-1);
        case 3
            newDwell(1);
        case 4
            newDwell(-1);
        case 5
            newState(+1);
        case 6
            newState(-1);
        case 7
            flank=ceil(1+flank*1.1);
        case 8
            flank=floor(flank/1.1-1);            
        case 9
            globalview=~globalview;
        case 10
            keepWatching=false; 
    end
    %disp(int2str([k length(ind)]))
end
% close windows that have been opened by the code
close(fig1);
delete(pem);
if(closeKBax)
    close(get(KBax,'parent'))
else
    axes(KBax)
    set(KBax,'nextplot',KBnextplot); % reset hold status
end

    % move focus one (ds= +-1) dwell time forwrad or backward
    function newDwell(ds) 
        oldInd=ind(k);
        newInd=ind(k)+ds;
        newInd=max(1,newInd);newInd=min(length(STATE),newInd);
        %disp([sWatch STATE(newInd)])
        
        sWatch=STATE(newInd);
        ind=find(STATE==sWatch); 
        set(fs,'xdata',t(s==sWatch),'ydata',RMS(s==sWatch));
        set(fss,'xdata',t(s==sWatch),'ydata',s(s==sWatch));
        
        set(pem,'xdata',W.est.sKaverage(sWatch),'ydata',W.est.sBaverage(sWatch));        
        k=find(ind==newInd);
    end
    function newState(ds) % shift to another state
        sWatch=sWatch+ds;
        if(sWatch>W.N)
            sWatch=1;
        elseif(sWatch<1)
            sWatch=W.N;
        end
        oldIndAll=ind;
        oldInd=ind(k);
        ind=find(STATE==sWatch);
        if(isempty(ind)) % the new state is unoccupied
            %disp(['state ' int2str(sWatch) ' not occupied!'])
            ind=oldIndAll;
            newState(ds);% try the next state in line
            return
        end
        k=find(ind>oldInd,1);
        if(isempty(k))
            k=length(ind);
        end
        set(fs,'xdata',t(s==sWatch),'ydata',RMS(s==sWatch));
        set(fss,'xdata',t(s==sWatch),'ydata',s(s==sWatch));
        set(pem,'xdata',W.est.sKaverage(sWatch),'ydata',W.est.sBaverage(sWatch));        
        
    end
end
