function s=inspectStates(W,X,s,fig)
% inspectStates(W,X,s,fig)
%
% a little script that plots coordinates X around each occurrence of state
% sWatch.
% W     : converged VBE object
% X     : x,y coordinates to plot
% s     : state sequence to use (default: Viterbi sequence)      
% fig   : figure in which to plot (default: current)

%% prepare parameters
if(~exist('s','var') || isempty(s))
    disp('computing Viterbi path')
    tic
    s=sViterbi(W.Q,W.qst);
    toc
end
if(exist('fig','var'))
    figure(fig)
end
ca=gca;
cf=gcf;
%% non-adjustable parameters
RMS=gaussFilter(X,3); % filter width
flank=20;             % number of points before and after each 2-dwell

%% plot data and states
[STATE,DWELL,Dstart,Dend]=getDwellTRJ(s);
sWatch=STATE(1);

t=(1:size(X,1))/30;
axes(ca),clf
f1=subplot(3,1,1);
hold on
plot(t,RMS,'k')
[m,n]=sort(W.RMSaverage);
fs=plot(t(s==sWatch),RMS(s==sWatch),'go');
plot(t,W.RMSaverage(s),'-r.')
[m,n]=sort(W.RMSaverage);
set(gca,'ytick',m,'yticklabel',int2str(n'),'ygrid','on')

ind=find(STATE==sWatch); 
f2=subplot(3,1,2);
hold on
plot(t,X(:,1),'-rs');
plot(t,X(:,2)+0,'-bs');

f3=subplot(3,1,3);hold on
set(f3,'ylim',[0.5 W.N+0.5],'ytick',1:W.N,'ygrid','on')
plot(t,s,'-k')
fss=plot(t(s==sWatch),s(s==sWatch),'go');
%% watch transition state events

keepWatching=true;
k=1;
dt=t(2)-t(1);
while(keepWatching)
    figure(cf)
    % find indices to watch
    ti=max(1,Dstart(ind(k))-flank):min(Dend(ind(k))+flank,length(s));
    
    i0=Dstart(ind(k));
    t0=t(i0);
    i1=Dend(ind(k));
    t1=t(i1);
        
    set([f1 f2 f3],'xtick',sort(union([t0-dt/2 t1+dt/2],t(ti([1 end])))),...
        'xgrid','on','xlim',t(ti([1 end])))
    set([f1 f2],'ylim',[-600 600])
    set([f1 f2 f3],'xticklabel',num2str(get(gca,'xtick')',5))
    box on
    axes(f1)
    ylabel('RMS [nm]')
    xlabel('t [s]')

    axes(f2)
    ylabel('coordinate [nm]')
    xlabel('t [s]')

    axes(f3)
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
              'flank +1','flank -1','quit');
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
            keepWatching=false; 
    end
    disp(int2str([k length(ind)]))
end
    % move focus one (ds= +-1) dwell time forwrad or backward
    function newDwell(ds) 
        oldInd=ind(k);
        newInd=ind(k)+ds;
        newInd=max(1,newInd);newInd=min(length(STATE),newInd);
        disp([sWatch STATE(newInd)])
        
        sWatch=STATE(newInd);
        ind=find(STATE==sWatch); 
        set(fs,'xdata',t(s==sWatch),'ydata',RMS(s==sWatch));
        set(fss,'xdata',t(s==sWatch),'ydata',s(s==sWatch));
                
        k=find(ind==newInd);
    end
    function newState(ds) % shift to another state
        sWatch=sWatch+ds;
        if(sWatch>W.N)
            sWatch=1;
        elseif(sWatch<1)
            sWatch=W.N;
        end
        oldInd=ind(k);
        ind=find(STATE==sWatch);
        k=find(ind>oldInd,1);
        if(isempty(k))
            k=length(ind);
        end
        set(fs,'xdata',t(s==sWatch),'ydata',RMS(s==sWatch));
        set(fss,'xdata',t(s==sWatch),'ydata',s(s==sWatch));
        
    end
end
