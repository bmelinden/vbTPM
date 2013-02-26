function xout=inspectLacData()
% xout=inspectLacData()
% function to read position data from lacdata.mat files, drift-correct and
% reorganize it, and also apply manual cuts to remove sticking etc.
%
% input is read from lacdata.mat in current directory, and the resulting
% selection is written to beadSelection.mat for later use.
% A cell vector of position traces is also returned.
%
% M.L. 2010-10-28

%% change-log

%% code
tic
disp('loading data...')
D=load('lacdata.mat');
toc
disp('concatenating data...')
xraw=TPM_RP_concatenate(D.lacX,D.lacY);
clear D;
toc
Fcut=0.05;  % Hz
Fsample=30; % Hz
tSigma=4;   % s (gaussian filter width)
x=cell(size(xraw));
RMS=x;
K=x;B=x;
disp('applying drift correction...')
for k=1:length(xraw)
    for j=1:length(xraw{k})
        x{k}{j}=BWdriftcorrect(xraw{k}{j},Fcut,Fsample);
        [RMS{k}{j},K{k}{j},B{k}{j}]=VBE1gaussFilter(x{k}{j},tSigma,Fsample);
    end
end
toc
file=exist('beadSelection.mat','file');
c=2;
if(file)
    c=menu('continue earlier selection?','yes','no');
    if(c==1)
        load beadSelection.mat
    end
end
if(c==2 || ~file )
    bindex=[];
    tinterval=[];
    for k=1:length(xraw)
        for j=1:length(xraw{k})
            bindex(end+1,:)=[k j];
            tinterval(end+1,:)=[1 length(K{k}{j})];
        end
    end
end

figure(1)
clf

t=(1:size(x{1}{1}))/Fsample;
y=150*rand(size(t));
subplot(2,1,1)
hold on
pp=plot(t,y);
ppi=plot(t,y,'r');
hold off
%set(gca,'ylim',[50 200])
xlabel('t [s]')
ylabel('RMS [nm]')

subplot(2,2,3)
[a,b]=hist(y,0:1:300);
ph=bar(b,a);
xlabel('RMS [nm]')

subplot(2,2,4)
[a,b]=hist3([K{1}{1} B{1}{1}],{0:0.005:1, (0:0.005:1)*1e-4});
pc=pcolor(b{1},b{2},a');
axis tight;
ylabel('B [nm^{-2}]');
xlabel('K');
set(pc,'edgecol','non')

if(~exist('n','var'))
    n=1;
    rmsBin=1;
    tSigma=4;
end
again=true;
while(again)
    k=bindex(n,1);
    j=bindex(n,2);
    ti=tinterval(n,1):tinterval(n,2);
    
    y=RMS{k}{j};
    t=(1:length(RMS{k}{j}))/Fsample;
    set(pp,'xdata',t,'ydata',y);    
    set(ppi,'xdata',t(ti),'ydata',y(ti));
    
    [a,b]=hist(RMS{k}{j}(ti),0:rmsBin:300); 
    set(ph,'xdata',b(2:end-1),'ydata',a(2:end-1)/sum(a(2:end-1)));    
    
    [a,b]=hist3([K{k}{j}(ti) B{k}{j}(ti)],{0:0.005:1, (0:0.005:1)*1e-4});
    set(pc,'cdata',a');
    pause(0.1)
    subplot(2,1,1)
    title([ int2str(bindex(n,:)) ', window: ' num2str(tSigma,3) ', bin: ' num2str(rmsBin,3)]) 
    
    c=menu('changes',...
        'quit','next','previous',...
        '+ window','- window','+ bin','- bin','discard trj',...
        'select time interval');
    if(c==1) % save and quit
        again=false;
        save beadSelection.mat bindex tinterval n tSigma rmsBin
    end
    if(c==2); % next trajectory
        n=n+1;
        n=n+(1-n)*(n>size(bindex,1));
    end
    if(c==3); % previous trajectory
        n=n-1;
        n=n+(size(bindex,1)-n)*(n<1);
    end
    if(c==4);tSigma=tSigma*sqrt(2);end
    if(c==5);tSigma=tSigma/sqrt(2);end
    if(c==6);rmsBin=rmsBin*sqrt(2);end    
    if(c==7);rmsBin=rmsBin/sqrt(2);end    
    if(c==8);
        bindex=[bindex(1:n-1,:); bindex(n+1:end,:)];
        tinterval=[tinterval(1:n-1,:); tinterval(n+1:end,:)];        
        n=max(1,min(n,size(bindex,1)));    
    end
    if(c==9);
       subplot(2,1,1)        
       t12=t(tinterval(n,:));
       y12=get(gca,'ylim');
       h=imrect(gca,[t12(1) y12(1) t12(2) y12(2)] );
       title('double-clik in rectangle to mark selection');
       L=wait(h);
       delete(h);
       tinterval(n,:)=[find(t>L(1),1,'first') find(t<(L(1)+L(3)),1,'last')];
       tinterval(n,2)=min(tinterval(n,2),length(K{k}{j}));
    end
    if( (c>=2 &&c<=5) || c==8 )
        k=bindex(n,1);
        j=bindex(n,2);
        [RMS{k}{j},K{k}{j},B{k}{j}]=VBE1gaussFilter(x{k}{j},tSigma,Fsample);
    end
end
close

xout=cell(1,size(bindex,1));
for m=1:size(bindex,1)
    k=bindex(m,1);
    j=bindex(m,2);
    ti=tinterval(m,1):tinterval(m,2);    
    xout{m}=x{k}{j}(ti,:);
end
