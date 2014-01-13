% s=manualSegmentML(file)
%
% function to do 3-segmentation on a data file file that contains the fields
% RMS2, RMS4, RMS8  : RMS traces at various windo sizes
% x                 : xy coordinates
%
% and possibly s, an old segmentation that you can continue working on.
%
% how it works: 
% 1) start the program with the appropriate file
% 2) if an earlier segmentation is on that file, you can reset it or modify it.
% 3) in working mode, you are given a choice of states to add. Add a new
% time interval of s=state by dragging the rectangle to enclose it. 
% 4) you can switch between seeing RMS traces and raw xy coordinates (the
% latter is good for spotting sticking events).
% 5) DONE ends the program, and asks about writing the resulting
% segmentation to the file. (If you ask no, the file is not changed).
%
%
% M.L. 2010-11-11

%% change-log
% M.L. 2010-11-16   : new rectangles start where the last one left off

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% manualSegmentML.m, part of the vbTPM package
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
function s=manualSegmentML(file)

%% parameters
fSample=30;
sLabels={'unlooped','mid loop','short loop','undecided','stuck','tracking error'};
%file='PUC306_50pMSJLacI_withRMS';
%% load data
disp(['reading data from ' file ]);
load(file);
%% initialize variables
t=(1:length(RMS2))/fSample;
RMSview=true;
resetSeq=~exist('s','var');
if(~resetSeq) % ask if reset anyway
c=menu('reset sequence?','yes','no');
switch c
    case 1
        resetSeq=true;
end
end
if(resetSeq)
    s=4*ones(size(t));
    s(1:5)=[1 2 3 5 6];
    RMS0=zeros(1,6);
end
%% draw segmentation & data
figure(1);
clf
hold on
f1=gca;
seg=zeros(size(RMS0));
for k=1:length(RMS0)
    ind=find(s==k);
    if(~isempty(ind))
        RMS0(k)=(mean(RMS2(ind))+mean(RMS4(ind))+mean(RMS8(ind)))/3;
    end
    RMS0(4)=0; % undecided
    SH(k)=plot(t([1 end]),[1 1]*RMS0(k),'k--');
    if(~isempty(ind))
        seg(k)=plot(t(ind),RMS0(k)*ones(size(ind)),'b.');
    else
        seg(k)=plot(0,0,'b.');        
    end
end
set(seg(4),'color','k') % undecided
set(seg(5),'color','k') % stuck
set(seg(6),'color','k') % tracker error

pr=plot(t,RMS2,'k');
pr(2)=plot(t,RMS4,'r');
pr(3)=plot(t,RMS8,'g');
segL=plot(t,RMS0(s),'b');
box on
set(f1,'xlim',t(end)*[-0.05 1.05])
RMSyrange=get(f1,'ylim');
XYyrange=max(RMSyrange)*[-1.5 1.5];
tLeft=-Inf;
%% segmentation
morechanges=true;
while(morechanges)
    for k=1:length(RMS0)
        ind=find(s==k);
        if(~isempty(ind))
            RMS0(k)=(mean(RMS2(ind))+mean(RMS4(ind))+mean(RMS8(ind)))/3;
        end
        RMS0(4)=0; % undecided
        set(SH(k),'ydata',[1 1]*RMS0(k));
        set(seg(k),'xdata',t(ind),'ydata',RMS0(k)*ones(size(ind)));
    end    
    set(segL,'ydata',RMS0(s));
    
    c=menu('new state:',sLabels{1},sLabels{2},sLabels{3},sLabels{4},sLabels{5},sLabels{6},...
        'DONE','RMS <-> xy');
    %disp(c)
    if(c==7)
         morechanges=false;
    elseif(c<=6) % modify state assignments
        y12=get(f1,'ylim');y12=0.9*y12+0.1*y12([2 1]);
       
        % old starting rectangle
        xRange=get(f1,'xlim');        
        t12=0.8*xRange+0.2*xRange([2 1]);
        
        if(tLeft>=xRange(1) &&  tLeft<=xRange(2)) % then continue from latest dwell
            if(exist('tLeft','var')) % then modify x-range to fit previous rectangle
                t12(1)=tLeft;
                t12(2)=max(t12(2),tLeft+abs(diff(t12)));
            end
        end
                
       h=imrect(gca,[t12(1) y12(1) diff(t12) diff(y12)] );
       title(['double-clik in rectangle to select new ' sLabels{c} ' points']);
       L=wait(h);
       
       tInterval=[L(1) L(1)+L(3)];
       ind=find((t>tInterval(1)).*(t<tInterval(2)));
       s(ind)=c;       
       delete(h);
       
       tLeft=tInterval(2); % left side of next rectangle
       
    elseif(c==8) % switch between displaying RMS and x,y
        RMSview=~RMSview;
        if(RMSview) % reset RMSview
           set(pr(1),'ydata',RMS2);
           set(pr(2),'ydata',RMS4);
           set(pr(3),'xdata',t,'ydata',RMS8);
           XYyrange=get(f1,'ylim');
           set(f1,'ylim',RMSyrange);           
           %y12=get(f1,'ylim');
           %set(f1,'ylim',[0 max(abs(y12))]);
        else % switch to xy view
           set(pr(1),'ydata',x(:,1));
           set(pr(2),'ydata',x(:,2));
           set(pr(3),'xdata',t(1),'ydata',0);     
           RMSyrange=get(f1,'ylim');
           set(f1,'ylim',XYyrange);           
           %y12=get(f1,'ylim');
           %set(f1,'ylim',[-1 1]*max(y12))
        end
        
    end
end
%% write result to file?
c=menu('add segmentation to file?','yes','no');
if(c==1)
    save(file,'s','RMS0','sLabels','-append')
end



