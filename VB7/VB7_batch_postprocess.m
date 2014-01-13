% gfig=VB7_batch_postprocess(runinputfile,xtrastates,savefile)
%
% GUI for manual state classification and analysis
% runinputfile  : VB7 runinput file (optional).
% xtrastates    : cell vector of extra state labels, in excess of
%                 'spurious','U', and 'L'. Default: {'M','B'}
%                give 0, false, or 'none' for no extra states.
% savefile      : filename where the analysis is to be stored (optional).
%                 It is possible to load and continue an old analysis, but
%                 the savefile must have been created with the same
%                 runinput file.

% M.L. 2011-10-21   : started GUI for postprocessing of single data set
%                      (all beads anayzed in a single runinput file).
% M.L. 2011-11-12   : % Removed unused states from standard analysis.
% M.L. 2011-12-02   : dealing with missing trajectories 
% M.L. 2012-01-23   : added dropdown menus to make it possible to switch
%                     what is plotted on the axes.
% M.L. 2012-02-02   : added filename to save analysis to as an argument, to
%                     make it easier to analyze the same data in different
%                     ways.
% M.L. 2012-03-13   : added fat calibration RMS line to RMS trajectory
%                     plot. Started on some toggle function for the RMS
%                     trajectory plot.
% M.L. 2012-03-14   : added optional conversion to parallell HMM models.
%                     Not completed, and deactivated by default.
% M.L. 2012-05-04   : made viterbi/sMaxP paths extend all the way to the
%                     last time point
% M.L. 2013-04-19   : added non-looping state as a default extra state
% M.L. 2013-04-25   : better GUI stability by disallowing callback
%                     interruption and qeueing 
% to do-list:
% - parallellize viterbi path calculations
% 
% - add command line options to get rid of popup menus
% - fix analysis file name: load and save different file names in GUI
% - converge dual HMM and save dwell sequences as part of analysis
%
% - display Ptot and dwell time for all states in the classification panel,
%   or change what is plotted: K,Tc,B,log(B),RMS,<Tdwell>,...
%   or have a 3D emission plot?
% - bigger font?
% - button to reset to default analysis for single state
% - plot calibration line in KB plot?
% - RMS histograms below KB plot?
% - help text
% - sort the states (genuine first, in decreasing RMS order), and make sure
%   the sorting information is not lost when saving/reloading
% - clustering algorithm on K,B,Tdwell?

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_batch_postprocess.m, postprocessing manager for the vbTPM package.
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
function gfig=VB7_batch_postprocess(runinputfile,xtrastates,savefile)
VB7_printGPL('VB7_batch_postprocess.m')

%% input parameters
% runinput file
%runinputfile='VB7runinput_PUC306_1pM_ds1';
if(~exist('runinputfile','var') || isempty(runinputfile))
    [runinputfile,runinputdir]=uigetfile('*.m','Select the runinput file for state classification');
else
    [runinputfile,runinputdir]=dirFileSplit(runinputfile);
end
if(strcmp(runinputfile(end-1),'.'))
    runinputfile=runinputfile(1:end-2);
end
startdir=pwd;
cd(runinputdir);
runinputdir=pwd;

h.tSigma=4;
h.RMSmin=80; % nm
h.runinputfile=runinputfile;
%% internal parameters
h.trjRMS_path_color{1}=[0.7 0.5 0]; % color of viterbi path 
h.trjRMS_path_color{2}=[0.7 0.7 0.7]; % color of seq. of mosty likely states

h.xyLabel={'K'};
h.xyFunction={@(Wt,Wc,k,cg)(Wt.est.sKaverage(k))};
h.xyLabel{end+1}='B [nm^{-2}]';
h.xyFunction{end+1}=@(Wt,Wc,k,cg)(Wt.est.sBaverage(k));

h.xyLabel{end+1}='K-K_{cal.}';
h.xyFunction{end+1}=@(Wt,Wc,k,cg)(Wt.est.sKaverage(k)-Wc.est.sKaverage(cg));
h.xyLabel{end+1}='B-B_{cal.} [nm^{-2}]';
h.xyFunction{end+1}=@(Wt,Wc,k,cg)(Wt.est.sBaverage(k)-Wc.est.sBaverage(cg));

h.xyLabel{end+1}='RMS [nm]';
h.xyFunction{end+1}=@(Wt,Wc,k,cg)(Wt.est.sRMS(k));
h.xyLabel{end+1}='RMS-RMS_{cal.} [nm]';
h.xyFunction{end+1}=@(Wt,Wc,k,cg)(Wt.est.sRMS(k)-Wc.est.sRMS(cg));

h.xyLabel{end+1}='corr. time: \tau_c [s]';
h.xyFunction{end+1}=@(Wt,Wc,k,cg)(Wt.est.sTC(k));
h.xyLabel{end+1}='corr. time: \tau_c-\tau_{c,cal.} [s]';
h.xyFunction{end+1}=@(Wt,Wc,k,cg)(Wt.est.sTC(k)-Wc.est.sTC(cg));

% The physically interesting dwell time is one where the spurious states
% are removed from the calculation, which means it has to be recomputed
% every time a state classification changes!
% h.xyLabel{7}='dwell time: <\tau>  [s]';
% h.xyFunction{7}=@(Wt,Wc,k,cg)(Wt.est.tD(k));
%% initialize internal state
% initial choices for x- and y- axes in scatter plot
h.Xf=1; % K-values
h.Yf=2; % B-values
% read HMM results and recompute some stuff
[HMMresultsfile,cal,trj]=VB7_batch_manage(runinputfile,'collect','save',true);
eval(runinputfile);
% filename for analysis file
if(exist('savefile','var') && ~isempty(savefile))
    if(length(savefile)>4)
        if(strcmp(savefile(end-3:end),'.mat')==0)
        savefile=[savefile '.mat'];
        end
    else
        savefile=[savefile '.mat'];        
    end
    HMManalysisfile_wpath=savefile;
else
    HMManalysisfile_wpath=[HMMresultsfile(1:end-4) '_analyzed.mat'];    
end
disp(['Saving results in ' HMManalysisfile_wpath ])

[HMManalysisfile,HMManalysispath]=dirFileSplit(HMManalysisfile_wpath);
h.HMManalysisfile=HMManalysisfile;
h.HMManalysispath=HMManalysispath;
% look for existing analysis file
continueoldanalysis=false;
if(exist(HMManalysisfile_wpath,'file'))
    disp('continue existing analysis? (choose from the menu window that just appeared)')
    choice=menu('continue existing analysis?','yes','no');
   if(choice==1)
      continueoldanalysis=true;
   end
   clear choice;
end

% debug-attempt
if(continueoldanalysis)
    R=load(HMManalysisfile_wpath);
    trj=R.trj;
    cal=R.cal;
end
    
for k=1:length(trj)
%%%%%%%for k=1:3 % debug trick to draw faster
    for b=1:length(trj{k})
        if(~isempty(trj{k}{b}) && ~isempty(cal{k}))
            tic
            W=trj{k}{b};
            % load data, using the path in the runinput file (which can be
            % changed by editing that file)
            filename=W.trjDataFile;
            filename=[source_path filesep filename(find(filename==filesep,1,'last')+1:end)];            
            foo=load(filename);
            x=foo.(W.looping_xyfield);
            if(W.driftcorrection)
                x=BWdriftcorrect(x,W.fCut,W.Wtrj.priorParameterOptions.fSample);
            end
            % extract sampling frequency
            fSample=W.Wtrj.priorParameterOptions.fSample;%/W.Wtrj.priorParameterOptions.downSample;
            if(isfield(h,'fSample') && h.fSample ~=fSample)
                error('All trajectories not analyzed with the same sampling frequency and/or downsampling.')
            end
            h.fSample=fSample;
            h.downSample=W.Wtrj.priorParameterOptions.downSample;
            if(1) %pre-compute RMS trajectory
                h.trjRMS{k}{b}=RMSKBgaussfilter(x,h.tSigma,h.fSample);
            end
            if(1)% pre-compute Viterbi path etc
                X=VB7_preprocess(x,W.Wtrj.priorParameterOptions.downSample);
                % compute all estimates
                W.Wtrj=VB7_VBEMiter(W.Wtrj,X,'estimate');
                %h.viterbi{k}{b}=uint8(W.Wtrj.est2.sViterbi);
                [h.viterbi_s{k}{b},~,h.viterbi_t{k}{b}]=getDwellTRJ(W.Wtrj.est2.sViterbi);
                h.viterbi_s{k}{b}(end+1)=h.viterbi_s{k}{b}(end);
                h.viterbi_t{k}{b}(end+1)=W.Wtrj.T;
                %h.viterbi{k}{b}=uint8(W.Wtrj.est2.sViterbi);                
                %h.sMaxP{k}{b}=uint8(W.Wtrj.est2.sMaxP);
                [h.sMaxP_s{k}{b},~,h.sMaxP_t{k}{b}]=getDwellTRJ(W.Wtrj.est2.sMaxP);
                h.sMaxP_s{k}{b}(end+1)=h.sMaxP_s{k}{b}(end);
                h.sMaxP_t{k}{b}(end+1)=W.Wtrj.T;
            end
            disp(['precomputed for trj ' int2str([k b]) ', ' num2str(toc) ' s.'])
        end
    end
end
clear foo X;


% initiaization
if(continueoldanalysis)
    R=load(HMManalysisfile_wpath);
else
    R=doDefaultAnalysis(trj,cal,h.RMSmin);
    % xtrastates & plot symbols
    if(~exist('xtrastates','var') || isempty(xtrastates) || (islogical(xtrastates) && extrastates==true))
        xtrastates={'M','B','NL'};
    elseif((isnumeric(xtrastates) && xtrastates==0) || (islogical(xtrastates) && extrastates==false) || sum(strcmp(xtrastates,'none'))>0)
        xtrastates={};
    end
    R.statelabels={'sp.','U','L'};
    for k=1:length(xtrastates)
        R.statelabels{end+1}=xtrastates{k};
    end
end
foo1='o^s>d<p';foo2='krgmc';
h.statesym={'.k'};
for k=2:length(R.statelabels)
    h.statesym{k}=[foo1(mod(k-1,length(foo1))+1) foo2(mod(k-1,length(foo2))+1)];
end
clear foo1 foo2
% make index map
h.kb=[];
h.Nmax=1;
for k=1:length(trj)
    for b=1:length(trj{k})
       h.kb(end+1,1:2)=[k b]; 
       if(R.hastrj{k}(b))
          h.Nmax=max(h.Nmax,trj{k}{b}.Wtrj.N); 
       end
    end
end
h.k=0;h.b=0;
% make sure to start at an existing trajectory
for m=1:size(h.kb,1)
    if(R.hastrj{h.kb(m,1)}(h.kb(m,2)))
        h.k=h.kb(m,1);
        h.b=h.kb(m,2);
        break
    end
end
if(h.k==0)
    disp('No converged trajectories found')
    return
end
%% lay out GUI figure
gfig=figure(1); % GUI figure
clf
set(gfig,'toolbar','auto')
% compute figure size
pos=get(gfig,'position');
set(0,'units','pixels');
spos=get(0,'screensize');

figx=pos(1);
figy=pos(2);
figw=max(pos(3),700); figh=pos(4); % present size

KBsize=400; % size of KB plot panel

rowh=35;
colw=25;
col1=30;
M=length(R.statelabels);
SCheight=25+(h.Nmax+3)*rowh;
SCwidth=10+col1+colw*(M+1);

figw=max(figw,KBsize+SCwidth+40);
figh=max(figh,max(KBsize+20+rowh,SCheight)+50+3*rowh);


%set(gfig,'pos',[pos(1:2) figw figh])
scrw=spos(3);scrh=spos(4);
set(gfig,'position',[max(0,min(figx,scrw-figw)) max(0,min(figy,scrh-figh-110)) figw figh])
clear spos scrw scrh
%% file name panel
h.fileNamePanel=uipanel('title','file names','parent',gfig,...
                'position',[10/figw 1-(20+3*rowh)/figh 1-20/figw (20+3*rowh)/figh]);
h.fileNameRI=uicontrol('style','text','string',['runiput file: ' runinputfile],'parent',h.fileNamePanel,...
          'position',[10 2*rowh (figw-20)/2-20 rowh],'HorizontalAlignment','left');
h.fileNameAN=uicontrol('style','text','string',['save file: ' HMManalysisfile],'parent',h.fileNamePanel,...
          'position',[10+(figw-20)/2 2*rowh (figw-20)/2-20 rowh],'HorizontalAlignment','left');   
h.saveAnalysis=uicontrol('style','pushbutton','string','save analysis','parent',h.fileNamePanel,...
           'position',[10+(figw-20)/2 5+1.1*rowh 100 0.8*rowh],'callback',@saveAnalysis);
h.fileNamePrev=uicontrol('style','pushbutton','string','prev.','parent',h.fileNamePanel,...
          'position',[10+(figw-20)/2 5+0.1*rowh 45 0.8*rowh],'callback',@prevTrj);
h.fileNameNext=uicontrol('style','pushbutton','string','next','parent',h.fileNamePanel,...
          'position',[65+(figw-20)/2 5+0.1*rowh 45 0.8*rowh],'callback',@nextTrj);
% information about this particular trajectory    
h.fileNameCAL=uicontrol('style','text','parent',h.fileNamePanel,...
          'string',['calibration: ' dirFileSplit(cal{h.k}.calDataFile)],...
          'position',[10 1*rowh (figw-20)/2-20 rowh],...
          'HorizontalAlignment','left');
h.fileNameTRJ=uicontrol('style','text','parent',h.fileNamePanel,...
          'string',['looping: ' dirFileSplit(trj{h.k}{h.b}.trjDataFile)],...
          'position',[10 0 (figw-20)/2-20 rowh],'HorizontalAlignment','left');
h.fileNameNUM=uicontrol('style','text','parent',h.fileNamePanel,...
          'position',[115+(figw-20)/2 0 150 rowh],...
          'string',['bead : ' int2str(h.k) '. subset : ' int2str(b) '.'],...
          'HorizontalAlignment','left');
%% trajectory action button group
h.KBactions=uipanel('title','','parent',gfig,'units','normalized',...
                    'position',[10/figw 1-(40+4*rowh)/figh KBsize/figw rowh/figh]);
h.excludeTRJ=uicontrol('style','checkbox','string','exclude trajectory',...
                       'position',[10 0.1*rowh 150 0.8*rowh],'parent',h.KBactions,...
                       'callback',@excludeTrj);
h.inspectTRJ=uicontrol('style','pushbutton','parent',h.KBactions,...
                       'position',[KBsize-10-100 0.1*rowh 100 0.8*rowh],...
                       'string','inspect','callback',@inspectTrj);                   
%% emission plot panel and RMS trajectory
h.KBpanel=uipanel('title','emission parameters','parent',gfig,'units','normalized',...
                'position',[10/figw 1-(60+4*rowh+KBsize)/figh KBsize/figw KBsize/figh]);
h.XYaxes=axes('parent',h.KBpanel,'position',[0.2 0.25 0.77 0.67],'box','on');



h.Xmenupanel=uipanel('title','x-axis','parent',h.KBpanel,'units','normalized',...
            'position',[0.02 5/KBsize 0.47 1.5*rowh/KBsize]);
h.Ymenupanel=uipanel('title','y-axis','parent',h.KBpanel,'units','normalized',...
            'position',[0.51 5/KBsize 0.47 1.5*rowh/KBsize]);

h.Xmenu=uicontrol('style','popupmenu','parent',h.Xmenupanel,'units','normalized',...
                  'position',[0.05 0.05 0.9 0.9],'callback',@setXYscatter,...
                  'string',h.xyLabel,'value',h.Xf);    
h.Ymenu=uicontrol('style','popupmenu','parent',h.Ymenupanel,'units','normalized',...
                  'position',[0.05 0.05 0.9 0.9],'callback',@setXYscatter,...
                  'string',h.xyLabel,'value',h.Yf);        
%% RMS plot control panel
h.RMSpanel=uipanel('title','RMS plot features','parent',gfig,'units','normalized',...
                    'position',[10/figw 1-(110+4.5*rowh+KBsize)/figh KBsize/figw 1.4*rowh/figh]);
h.RMS_rms=uicontrol('style','checkbox','string','RMS trj','Value',true,...
                       'position',[10 0.1*rowh 100 0.8*rowh],'parent',h.RMSpanel,...
                       'callback',@toggleRMS);
h.RMS_viterbi=uicontrol('style','checkbox','string','Viterbi path','Value',true,...
                       'position',[90 0.1*rowh 100 0.8*rowh],'parent',h.RMSpanel,...
                       'callback',@toggleRMS);
h.RMS_Pmax=uicontrol('style','checkbox','string','argmax p(s_t)','Value',true,...
                       'position',[200 0.1*rowh 120 0.8*rowh],'parent',h.RMSpanel,...
                       'callback',@toggleRMS);
h.RMS_RMSlevels=uicontrol('style','checkbox','string','levels','Value',true,...
                       'position',[320 0.1*rowh 70 0.8*rowh],'parent',h.RMSpanel,...
                       'callback',@toggleRMS);

% plot RMS
h.RMSfig=figure(2); 
clf
set(h.RMSfig,'toolbar','auto')
h=plotNewRMS(h,trj,cal,R);
figure(h.RMSfig)
axis tight

% plot states
h.KBdots=[];
h.KBtext=[];
h.backgroundToggle=zeros(1,length(R.statelabels)+1);


h=plotAllStates(h,R,trj,cal);
h=plotNewStates(h,trj,cal,R);
h=updateKBsymbols(h,trj,R);

% see http://www.mathworks.se/help/techdoc/creating_guis/f15-998661.html#f15-1000880
% for some useful hints on axes
%% state classification panel
h.SCpanel=uipanel('title','state classification','parent',gfig,'units','normalized',...
    'position',[(KBsize+30)/figw 1-(40+3*rowh+SCheight)/figh SCwidth/figw SCheight/figh]);
h.SCtext=[];
for k=1:h.Nmax
    % button group for each state
    h.SCstateGroup(k)=uibuttongroup('visible','on',...
        'SelectionChangeFcn',@stateChange,'parent',h.SCpanel,...
        'position',[5/SCwidth (rowh*(h.Nmax-k)+5)/SCheight 1-10/SCwidth 0.9*rowh/SCheight],...
        'bordertype','line','tag',['state ' int2str(k)],'userdata',0);
    % state number to the left in each row
    h.SCtext(end+1)=uicontrol('style','text','string',int2str(k),...
        'parent',h.SCstateGroup(k),...
        'position',[10 0.15*rowh 0.8*col1 rowh/2],...
        'HorizontalAlignment','right');
    for m=1:M % row of radio buttons for each state
        h.SCstateButton(k,m)=uicontrol('style','radio',...
            'position',[col1+colw*(m-0.5) rowh/4 colw rowh/2],...
            'parent',h.SCstateGroup(k),...
            'tag',int2str([k m]),'userdata',[k m]);
    end
end

% scatter plot toggle buttons
h.SCtextpanel=uipanel('visible','on','parent',h.SCpanel,...
    'position',[5/SCwidth rowh*(h.Nmax+1.5)/SCheight 1-10/SCwidth 2.3*rowh/SCheight],...
    'bordertype','line','title','scatter plot');
% make room for state label header
%set(h.SCstateGroup(1),'position',get(h.SCstateGroup(1),'position')+[0 0 0 rowh/SCheight]);
for m=1:M    
    h.SCscatterToggle(m)=uicontrol('style','togglebutton','string',R.statelabels{m},...
        'position',[col1+colw*(m-0.5) 0.1*rowh colw 0.7*rowh],...%[colw*(m-0.5) 2 colw rowh/2],...
        'parent',h.SCtextpanel(1),'callback',@scatterPlotToggle);    
    h.SCstateLabels(m)=uicontrol('style','text','string',R.statelabels{m},...
        'position',[col1+colw*(m-0.5) 0.9*rowh colw rowh/2],...%[colw*(m-0.5) 2 colw rowh/2],...
        'parent',h.SCstateGroup(1)); 
end
h.SCscatterToggle(M+1)=uicontrol('style','togglebutton','string','cal.',...
    'position',[10 0.1*rowh col1 0.7*rowh],...%[col1+colw*(m-0.5) rowh/4 colw rowh/2],...%[colw*(m-0.5) 2 colw rowh/2],...
    'parent',h.SCtextpanel(1),'callback',@scatterPlotToggle);

% polygon selection button
h.SCpolygonConvertButton=uicontrol('style','pushbutton','parent',h.SCtextpanel,...
    'position',[10 0.9*rowh SCwidth-30 0.7*rowh],...
    'callback',@massConversion,'string','convert many');
% turn on and correctly label the needed number of buttons
updateClassificationButtons(h,R.trjstates{h.k}{h.b},R.includetrj{h.k}(h.b));            
%% initialize GUI data
% put data into the GUI data structure
myhandle=guihandles(gfig);
myhandle.h=h;
myhandle.R=R;
myhandle.trj=trj;
myhandle.cal=cal;
guidata(gfig,myhandle);

% make all callback functions uninterruptible
cp=get(gfig,'Children'); % all panels
for np=1:length(cp)
    cc=get(cp,'Children'); % all panel component
    for nc=1:length(cc)
        set(cc{nc},'interruptible','off')
        set(cc{nc},'busyaction','cancel')
    end
end
end
%% callback functions
function toggleRMS(hObject,eventdata)
% replot RMS plot with new toggle values
%disp('Death and horror: RMS trace toggle not implemented!')
myhandle=guidata(gcbo);
h=myhandle.h;
R=myhandle.R;
trj=myhandle.trj;

h=toggleRMS_doit(h);
h=updateKBsymbols(h,trj,R);
h.trjRMS_path_color{1}=[0.7 0.5 0]; % color of viterbi path 

myhandle.h=h;
guidata(gcbo,myhandle);

end
function setXYscatter(hObject,eventdata)
% changes the property being plotted on the scatter plot x axis, after
% selection has changed in the corresponding dropdown menu.
myhandle=guidata(gcbo);
h=myhandle.h;
R=myhandle.R;
trj=myhandle.trj;
cal=myhandle.cal;


xRange=get(h.XYaxes,'xlim');
yRange=get(h.XYaxes,'ylim');

if(hObject==h.Xmenu)
    h.Xf=get(hObject,'value'); % new function number
elseif(hObject==h.Ymenu)
    h.Yf=get(hObject,'value'); % new function number
else
    error('setXYscatter called by other than Xmenu or Ymenu!')
end

c=get(h.XYaxes,'children');
delete(c);                 % delete everything in the scatter plot
clear c;
h.KBdots=[];h.KBtext=[];   % throw away handles to deleted objects

% redraw scatter plot from scratch
h=plotAllStates(h,R,trj,cal);
h=plotNewStates(h,trj,cal,R);
h=updateKBsymbols(h,trj,R);

axis tight
if(hObject==h.Xmenu)
set(h.XYaxes,'ylim',yRange); % put in the y-range we had in the beginning
else
set(h.XYaxes,'xlim',xRange); % put in the y-range we had in the beginning    
end

myhandle.h=h;myhandle.R=R;
guidata(gcbo,myhandle); % check in GUI data again
end
function saveAnalysis(hObject,eventdata)
% save the present sta.te of analysis
myhandle=guidata(gcbo);
h=myhandle.h;
R=myhandle.R;
trj=myhandle.trj;
cal=myhandle.cal;
%% do some standard analysis
M=length(R.statelabels);
% emission helper variables, see ML19, eqs. C33-36
clear Wtrj;
Wtrj.N=M;
Wtrj.E.wA=zeros(M,M);
Wtrj.E.ds_1=zeros(1,M);
Wtrj.E.C=zeros(1,M);Wtrj.E.V=zeros(1,M);
Wtrj.E.M=zeros(1,M);Wtrj.E.U=zeros(1,M);
Wtrj.T=0;

Wcal.E.C=0;Wcal.E.V=0;
Wcal.E.M=0;Wcal.E.U=0;Wcal.T=0;

for k=1:length(trj)
    if(sum(R.includetrj{k})>0) % then include this calibration curve
        W=cal{k}.Wcal;
        Wcal.E.C=Wcal.E.C+W.E.C(R.calGenuine(k));
        Wcal.E.V=Wcal.E.V+W.E.V(R.calGenuine(k));
        Wcal.E.M=Wcal.E.M+W.E.M(R.calGenuine(k));
        Wcal.E.U=Wcal.E.U+W.E.U(R.calGenuine(k));
        Wcal.T=Wcal.T+W.T;
    end
    for b=1:length(trj{k})
        if(R.includetrj{k}(b))
            
            W=trj{k}{b}.Wtrj;
            ind=R.trjstates{k}{b};
            N=W.N;
            for i=1:N
                si=ind(i);
                Wtrj.T=Wtrj.T+W.T;
                % emission parameters
                Wtrj.E.C(si)=Wtrj.E.C(si)+W.E.C(i);
                Wtrj.E.V(si)=Wtrj.E.V(si)+W.E.V(i);
                Wtrj.E.M(si)=Wtrj.E.M(si)+W.E.M(i);
                Wtrj.E.U(si)=Wtrj.E.U(si)+W.E.U(i);
                % summed transition matrix (no prior)
                Wtrj.E.ds_1(si)=Wtrj.E.ds_1(si)+W.E.ds_1(i);
                for j=1:N
                    sj=ind(j);
                    Wtrj.E.wA(si,sj)=Wtrj.E.wA(si,sj)+W.E.wA(i,j); % see eq C18
                end
            end
        end
    end
end
f=W.priorParameterOptions.fSample;     % original sampling freq.
ds=W.priorParameterOptions.downSample; % downsampling factor

% remove unpopulated states
ind=find(diag(Wtrj.E.wA)>0); % states to keep in the analysis
if(ismember(1,ind))
    ind=[setdiff(ind,1);1];   % put spurious states last, exclude unoccupied states
end
R.statelabels_Wxxx={R.statelabels{ind}};


M=length(ind);
Wtrj.N=M;
Wtrj.E.ds_1= Wtrj.E.ds_1(ind);
Wtrj.E.wA  = Wtrj.E.wA(ind,ind);
Wtrj.E.M   = Wtrj.E.M(ind);
Wtrj.E.C   = Wtrj.E.C(ind);
Wtrj.E.V   = Wtrj.E.V(ind);
Wtrj.E.U   = Wtrj.E.U(ind);

Wtrj.M.wA=Wtrj.E.wA;
Wtrj.M.wPi=Wtrj.E.ds_1;
% summed emission parameters (without priors)
Wtrj.M.n=Wtrj.E.M;
Wtrj.M.c=Wtrj.E.C-Wtrj.E.U.^2./Wtrj.E.V;
Wtrj.M.v=Wtrj.E.V;
Wtrj.M.mu=Wtrj.E.U./Wtrj.E.V;

Wcal.M.n=Wcal.E.M;
Wcal.M.c=Wcal.E.C-Wcal.E.U.^2./Wcal.E.V;
Wcal.M.v=Wcal.E.V;
Wcal.M.mu=Wcal.E.U./Wcal.E.V;

% interesting averages, transitions
Ttot=sum(Wtrj.M.wA,1); % in downsampled units
Wtrj.est.Paverage=Ttot/sum(Ttot); %relative probabilities of genuine states
% in upsampled units
Pcal=Wcal.M.n/Wcal.T/ds;
% exclude spurious states from rate calculations
A=rowNormalize(Wtrj.M.wA(1:end-1,1:end-1)); % accumulated transition rates, per downsampled time step
Q=logm(A)*f/ds;     % transition rates, per second
tDA=(1./(1-diag(A))*ds/f)';
tDQ=(-1./diag(Q))';
%R.A=A;
%R.Q=Q;


% emission parameters
K=Wtrj.M.mu;
dK=sqrt(Wtrj.M.c./(2*Wtrj.M.v.*(Wtrj.M.n-0.5)));
Kcal=Wcal.M.mu;
dKcal=sqrt(Wcal.M.c./(2*Wcal.M.v.*(Wcal.M.n-0.5)));
B=(Wtrj.M.n+0.5)./Wtrj.M.c;
dB=sqrt((Wtrj.M.n+0.5))./Wtrj.M.c;
Bcal=(Wcal.M.n+0.5)./Wcal.M.c;
dBcal=sqrt((Wcal.M.n+0.5))./Wcal.M.c;

RMS=sqrt(1./B./(1-K.^2));
RMScal=sqrt(1./Bcal./(1-Kcal.^2));

Tc=-ds/f./log(K);
TcCal=-ds/f/log(Kcal);
%% display analysis results
disp(' ')
disp('Analysis according to current state classification:')
disp('(Assuming units of s and nm.)')

for k=0:M+1;fprintf('----------');end;fprintf('\n');
fprintf('          ');
fprintf('%10s','cal.',R.statelabels_Wxxx{1:end})
fprintf('\n')

fprintf('%10s','Ptot')
fprintf('     %1.3f',[Pcal Wtrj.est.Paverage ])
fprintf('\n')

fprintf('%10s','<Tdwell>')
fprintf('%10s','---')
fprintf('%10.3f',tDA)
fprintf('%10s','---')

fprintf(' s (*)\n')
fprintf('%10s','<Tdwell>')
fprintf('%10s','---')
fprintf('%10.3f',tDQ)
fprintf('%10s','---')
fprintf(' s (*)\n')

fprintf('%10s','<K>')
fprintf('%10.3f',[Kcal K ])
fprintf('\n')
%fprintf('%10s','std[K]')       % should rather boot-strap the trajectories
%fprintf('%10.2e',[dKcal dK ])
%fprintf('\n')

Bexp=floor(mean(log10([Bcal B])));
fprintf('%10s','<B>')
fprintf('%10.3f',[Bcal B ]/10^(Bexp))
fprintf(' 1e%d/nm^2\n',Bexp)
%fprintf('%10s','std[B] ')      % should rather boot-strap the trajectories
%fprintf('%10.3f',[dBcal dB ]/1e-5)
%fprintf(' 1e-5/nm^2\n')

fprintf('%10s','RMS')
fprintf('%10.3f',[RMScal RMS ])
fprintf(' nm\n')

fprintf('%10s','Tcorr')
fprintf('%10.3f',[TcCal Tc ])
fprintf(' s\n')
for k=0:M+1;fprintf('----------');end;fprintf('\n');
disp('(*): Average dwell times computed by two different methods,')
disp('     not including the spurious states. ')
for k=0:M+1;fprintf('----------');end;fprintf('\n');

disp('Transition matrix (non-spurious states only) A: ')
fprintf('%20s',' ')
for mm=1:M-1
    fprintf('%12s',R.statelabels_Wxxx{mm})
end
fprintf('\n')
for mm=1:M-1
    fprintf('%20s',' [ ')
    fprintf('%12.2e',A(mm,:))
    fprintf(' ]\n')
end
fprintf('\n')

disp('Transition rate matrix (non-spurious states only) Q=logm(A)/dt: ')
fprintf('%20s',' ')
for mm=1:M-1
    fprintf('%12s',R.statelabels_Wxxx{mm})
end
fprintf('\n')
for mm=1:M-1
    fprintf('%20s',' [ ')
    fprintf('%12.2e',Q(mm,:))
    fprintf(' ]\n')
end
fprintf('\n')
%% save results
%R.Wtrj=Wtrj;
%R.Wcal=Wcal;
R.trj=trj;
R.cal=cal;
save([h.HMManalysispath filesep h.HMManalysisfile],...
    '-struct','R');
disp(['wrote analysis to '  [h.HMManalysispath filesep h.HMManalysisfile]]);
cc=menu('create factorial models?','yes','no');
if(cc==1)
    VB7_kineticAnalysis(h.runinputfile,[h.HMManalysispath filesep h.HMManalysisfile],false,true);
end
end
function nextTrj(hObject,eventdata)
% change to next trajectory
myhandle=guidata(gcbo);
h=myhandle.h;
R=myhandle.R;
trj=myhandle.trj;
cal=myhandle.cal;
% sort states in current models

%[a,b,ind]=sortModel(trj{h.k}{h.b}.Wtrj,R.trjstates{h.k}{h.b});
%trj{h.k}{h.b}.Wtrj=a;
%R.trjstates{h.k}{h.b}=b;
%h.viterbi{h.k}{h.b}=ind(h.viterbi{h.k}{h.b});
%h.sMaxP{h.k}{h.b}=ind(h.sMaxP{h.k}{h.b});

% move on to next existing trajectory
mMax=length([R.hastrj{:}]);
for mm=1:mMax
    h.b=h.b+1;
    % check if
    if(h.b>length(trj{h.k}))
        h.b=1;
        h.k=h.k+1;
    end
    if(h.k>length(trj))
        h.k=1;
    end
    if(R.hastrj{h.k}(h.b))
        break
    end
end
updateTrjText(h,cal,trj);
updateClassificationButtons(h,R.trjstates{h.k}{h.b},R.includetrj{h.k}(h.b));
h=plotNewRMS(h,trj,cal,R);
h=plotNewStates(h,trj,cal,R);
h=updateKBsymbols(h,trj,R);

myhandle.h=h;
myhandle.R=R;
myhandle.trj=trj;
guidata(gcbo,myhandle);
end
function prevTrj(hObject,eventdata)
% change to next trajectory
myhandle=guidata(gcbo);
h=myhandle.h;
R=myhandle.R;
trj=myhandle.trj;
cal=myhandle.cal;
% sort states in current models
%[a,b,ind]=sortModel(trj{h.k}{h.b}.Wtrj,R.trjstates{h.k}{h.b});
%trj{h.k}{h.b}.Wtrj=a;
%R.trjstates{h.k}{h.b}=b;
%h.viterbi{h.k}{h.b}=ind(h.viterbi{h.k}{h.b});
%h.sMaxP{h.k}{h.b}=ind(h.sMaxP{h.k}{h.b});

% move on to next previous existing trajetory
mMax=length([R.hastrj{:}]);
for mm=1:mMax    
    h.b=h.b-1;
    % check if
    if(h.b==0)
        h.k=h.k-1;
        if(h.k==0)
            h.k=length(trj);
        end
        h.b=length(trj{h.k});
    end
    if(R.hastrj{h.k}(h.b))
        break
    end
end

updateTrjText(h,cal,trj);
updateClassificationButtons(h,R.trjstates{h.k}{h.b},R.includetrj{h.k}(h.b));
h=plotNewRMS(h,trj,cal,R);
h=plotNewStates(h,trj,cal,R);
h=updateKBsymbols(h,trj,R);

myhandle.h=h;
myhandle.R=R;
myhandle.trj=trj;
guidata(gcbo,myhandle);
end
function scatterPlotToggle(hObject,source)
myhandle=guidata(gcbo);
h=myhandle.h;
R=myhandle.R;

ss=find(hObject==h.SCscatterToggle);
if(ss<=length(R.statelabels))
    h.backgroundToggle(ss)=~h.backgroundToggle(ss);
else
    h.backgroundToggle(end)=~h.backgroundToggle(end);
end
    
h=updateBGscatter(h,R,myhandle.trj);

myhandle.h=h;
myhandle.R=R;
guidata(gcbo,myhandle);
end
function stateChange(hObject,source)
% react to change of state classification
myhandle=guidata(gcbo);
h=myhandle.h;
R=myhandle.R;
trj=myhandle.trj;
if(R.includetrj{h.k}(h.b))
    changedState=find(h.SCstateGroup==hObject);
    newLabel=find(source.NewValue==h.SCstateButton(changedState,:));
    %set(h.SCstateButton(changedState,:),'value',0);
    %set(h.SCstateButton(changedState,newLabel),'value',1);
    R.trjstates{h.k}{h.b}(changedState)=newLabel;
    h=updateKBsymbols(h,myhandle.trj,R);
else
    updateClassificationButtons(h,R.trjstates{h.k}{h.b},R.includetrj{h.k}(h.b));
    menu('cannot change state in excluded trajectory!','OK')
end
h=updateBGscatter(h,R,trj);

myhandle.h=h;
myhandle.R=R;
myhandle.trj=trj;
guidata(gcbo,myhandle);
end
function excludeTrj(hObject,source)
% flag this trajectory to not include in analysis
myhandle=guidata(gcbo);
h=myhandle.h;
R=myhandle.R;

if(get(hObject,'value')==1) % then the trajectory should be marked to be excludeed
    % update flag
    R.includetrj{h.k}(h.b)=false;
    h=updateBGscatter(h,R,myhandle.trj);
    % set all states to spurious
    R.trjstates{h.k}{h.b}(1:end)=1;
else % marked to be included again
    % update flag
    if(R.hastrj{h.k}(h.b))
        R.includetrj{h.k}(h.b)=true;
        % make KB background white again
        set(h.XYaxes,'color','w')
        % perform default analysis to classify states
        [~,R.trjstates{h.k}{h.b}]=defaultAnalysis(...
            myhandle.cal{h.k}.Wcal,myhandle.trj{h.k}{h.b}.Wtrj,h.RMSmin);
        h=plotNewStates(h,myhandle.trj,myhandle.cal,R);
        h=updateKBsymbols(h,myhandle.trj,R);
    else
        menu('no converged model available for this trajectory','OK')
        
    end
end
h=updateKBsymbols(h,myhandle.trj,R);
updateClassificationButtons(h,R.trjstates{h.k}{h.b},R.includetrj{h.k}(h.b));

myhandle.h=h;
myhandle.R=R;
guidata(gcbo,myhandle);
end
function inspectTrj(hObject,source)
% branch out to the inspect states function
myhandle=guidata(gcbo);
h=myhandle.h;
R=myhandle.R;
W=myhandle.trj{h.k}{h.b}.Wtrj;
Wc=myhandle.cal{h.k}.Wcal;
% load data
a=load(myhandle.trj{h.k}{h.b}.trjDataFile);
if(isfield(myhandle.trj{h.k}{h.b},'looping_xyfield'))
    looping_xyfield=myhandle.trj{h.k}{h.b}.looping_xyfield;    
    fCut=myhandle.trj{h.k}{h.b}.fCut;
    driftcorrection=myhandle.trj{h.k}{h.b}.driftcorrection;
else
   eval(h.runinputfile)
end
fSample=W.priorParameterOptions.fSample;
ds=W.priorParameterOptions.downSample;
x=a.(looping_xyfield);
clear a
if(driftcorrection)
    x=BWdriftcorrect(x,fCut,fSample);
else
end
X=VB7_preprocess(x,ds,fSample);
W=VB7_VBEMiter(W,X,'estimate');

% new function is able to handle arbitrary scatter plots:
xsc=h.xyFunction{h.Xf}(W,Wc,1:W.N,R.calGenuine(h.k));
ysc=h.xyFunction{h.Yf}(W,Wc,1:W.N,R.calGenuine(h.k));
VB7_inspectStates(W,x,W.est2.sViterbi,h.XYaxes,xsc,ysc);

end
function massConversion(hObject,source)
% reclassify states based on interactive polygon
myhandle=guidata(gcbo);
h=myhandle.h;
R=myhandle.R;
trj=myhandle.trj;
cal=myhandle.cal;

% display instructions
disp('Convert many states at once:')
disp('1. Create a polygon by clicking in the emission plot.')
disp('   Double-click to place the last vertex')
disp('2. Move the vertices to fine-tune your selection')
disp('3. Select target state from the pop-up menu')
disp('All visible states will be converted to the target state.')
disp('')
% draw polygon
p=impoly(h.XYaxes,'closed',true);

% select new state or cancel
newstate=menu('target state:',...
    {R.statelabels{:},'cancel'});

% perform conversion
if(newstate<=length(R.statelabels)) % then update visible states
    % get polygon vertices
    pch=get(p,'children');
    XV=get(pch(end),'xdata');
    YV=get(pch(end),'ydata');
    % update states in active trajectory
    if(R.includetrj{h.k}(h.b))
        %K=trj{h.k}{h.b}.Wtrj.est.sKaverage;
        %B=trj{h.k}{h.b}.Wtrj.est.sBaverage;
        Wt=trj{h.k}{h.b}.Wtrj;
        Wc=cal{h.k}.Wcal;
        K=h.xyFunction{h.Xf}(Wt,Wc,1:Wt.N,R.calGenuine(h.k));
        B=h.xyFunction{h.Yf}(Wt,Wc,1:Wt.N,R.calGenuine(h.k));
        inp=inpolygon(K,B,XV,YV);
        R.trjstates{h.k}{h.b}(inp)=newstate;
    end
    % update states in background scatter plot
    vis=h.backgroundToggle; % visible states
    for k=1:length(trj)
        for b=1:length(trj{k})
            if(~(k==h.k && b==h.b) && R.includetrj{k}(b))                
                %inp=inpolygon(trj{k}{b}.Wtrj.est.sKaverage,trj{k}{b}.Wtrj.est.sBaverage,XV,YV);
                Wt=trj{k}{b}.Wtrj;
                Wc=cal{k}.Wcal;
                K=h.xyFunction{h.Xf}(Wt,Wc,1:Wt.N,R.calGenuine(k));
                B=h.xyFunction{h.Yf}(Wt,Wc,1:Wt.N,R.calGenuine(k));
                inp=inpolygon(K,B,XV,YV);
                
                for s=1:trj{k}{b}.Wtrj.N
                    if(inp(s) &&  vis(R.trjstates{k}{b}(s)))
                        R.trjstates{k}{b}(s)=newstate;
                    end
                end
                
            end
        end
    end
    %h.backgroundToggle(newstate)=1;
    h=updateBGscatter(h,R,trj);
    h=updateKBsymbols(h,trj,R);
    updateClassificationButtons(h,R.trjstates{h.k}{h.b},R.includetrj{h.k}(h.b));
else
    disp('Conversion canceled.')
end
delete(p);


myhandle.h=h;
myhandle.R=R;
guidata(gcbo,myhandle);
end
%% plot update functions
function h=toggleRMS_doit(h)
% update on/off status of graphs in the RMS plot
if(get(h.RMS_viterbi,'Value')==false)
   set(h.trjRMS_path(1),'linestyle','none');
else
   set(h.trjRMS_path(1),'linestyle','-');    
end
if(get(h.RMS_rms,'Value')==false)
   set(h.trjRMS_rms,'linestyle','none');
else
   set(h.trjRMS_rms,'linestyle','-');    
end
if(get(h.RMS_Pmax,'Value')==false)
   set(h.trjRMS_path(2),'linestyle','none');
else
   set(h.trjRMS_path(2),'linestyle','--');    
end
if(get(h.RMS_RMSlevels,'Value')==false)    
   set([h.trjRMS_line h.calRMS_line],'linestyle','none');
else
   set([h.trjRMS_line h.calRMS_line],'linestyle','--');
end

end
function h=plotAllStates(h,R,trj,cal)
% scatter plot of all states
axes(h.XYaxes);
hold on
h.trjplot=cell(size(trj));
h.calplot_gen=zeros(1,length(trj));
h.calplot_spr=[];
for k=1:length(trj)
    h.trjplot{k}=cell(size(trj{k}));
    for b=1:length(trj{k})
        if(R.hastrj{k}(b))
            %K=trj{k}{b}.Wtrj.est.sKaverage;
            %B=trj{k}{b}.Wtrj.est.sBaverage;
            K=h.xyFunction{h.Xf}(trj{k}{b}.Wtrj,cal{k}.Wcal,1:trj{k}{b}.Wtrj.N,R.calGenuine(k));
            B=h.xyFunction{h.Yf}(trj{k}{b}.Wtrj,cal{k}.Wcal,1:trj{k}{b}.Wtrj.N,R.calGenuine(k));
            h.trjplot{k}{b}=zeros(1,trj{k}{b}.Wtrj.N);
            for s=1:trj{k}{b}.Wtrj.N
                h.trjplot{k}{b}(s)=plot(K(s),B(s),'.');
            end
        end
    end
    if(0<sum(R.hastrj{k})) % only plot calibration if there are any trajectories
        %Kc=cal{k}.Wcal.est.sKaverage;
        %Bc=cal{k}.Wcal.est.sBaverage;
        g=R.calGenuine(k);
        Kc=h.xyFunction{h.Xf}(cal{k}.Wcal,cal{k}.Wcal,1:cal{k}.Wcal.N,g);
        Bc=h.xyFunction{h.Yf}(cal{k}.Wcal,cal{k}.Wcal,1:cal{k}.Wcal.N,g);
        h.calplot_gen(k)=plot(Kc(g),Bc(g),'.');
        foo=plot(Kc([1:g-1 g+1:end]),Bc([1:g-1 g+1:end]),'.');
        if(~isempty(foo))
            h.calplot_spr(k)=foo;
            %h.calplot_spr(end+1)=foo;
        end
    end
end
xlabel(h.xyLabel{h.Xf})
ylabel(h.xyLabel{h.Yf})    
h=updateBGscatter(h,R,trj);
end
function h=updateBGscatter(h,R,trj)
% update symbols on background scatter plot
M=length(R.statelabels);

% calibration genuine states
if(h.backgroundToggle(M+1))
    set(h.calplot_gen(h.calplot_gen>0),'marker','o','color','b');
else
    set(h.calplot_gen(h.calplot_gen>0),'marker','none' );
end
% calibration spurious states
if(h.backgroundToggle(1))
    set(h.calplot_spr(h.calplot_spr>0),'marker','.','color',0.5*[1 1 1]);    
else
    set(h.calplot_spr(h.calplot_spr>0),'marker','none');
end

for k=1:length(trj)
    for b=1:length(trj{k})
        if(R.includetrj{k}(b))
            for s=1:trj{k}{b}.Wtrj.N
                ss=R.trjstates{k}{b}(s);
                if(h.backgroundToggle(ss))
                    set(h.trjplot{k}{b}(s),'marker',h.statesym{ss}(1),...
                        'color',h.statesym{ss}(2),'markerface','none');
                    col=0.5*get(h.trjplot{k}{b}(s),'color')+0.5;
                    set(h.trjplot{k}{b}(s),'color',col);
                else
                    set(h.trjplot{k}{b}(s),'marker','none');
                end
            end
        else
            set(h.trjplot{k}{b},'marker','none');
        end
    end
end
end
function h=plotNewStates(h,trj,cal,R)
% produce a scatter plot of active states with text.

% scatterplot:
delete(h.KBdots);
delete(h.KBtext);
h.KBdots=[];
h.KBtext=[];
axes(h.XYaxes);
%K=trj{h.k}{h.b}.Wtrj.est.sKaverage;
%B=trj{h.k}{h.b}.Wtrj.est.sBaverage;
K=h.xyFunction{h.Xf}(trj{h.k}{h.b}.Wtrj,cal{h.k}.Wcal,1:trj{h.k}{h.b}.Wtrj.N,R.calGenuine(h.k));
B=h.xyFunction{h.Yf}(trj{h.k}{h.b}.Wtrj,cal{h.k}.Wcal,1:trj{h.k}{h.b}.Wtrj.N,R.calGenuine(h.k));
for k=1:length(K)
    h.KBdots(k)=plot(K(k),B(k),'.');
    h.KBtext(k)=text(K(k),B(k),[' ' int2str(k)],...
        'VerticalAlignment','bottom','HorizontalAlignment','left');
    hold on
end
% calibration state will remain untouched
Kc=h.xyFunction{h.Xf}(cal{h.k}.Wcal,cal{h.k}.Wcal,R.calGenuine(h.k),R.calGenuine(h.k));
Bc=h.xyFunction{h.Yf}(cal{h.k}.Wcal,cal{h.k}.Wcal,R.calGenuine(h.k),R.calGenuine(h.k));
%%Kc=cal{h.k}.Wcal.est.sKaverage(R.calGenuine(h.k));
%%Bc=cal{h.k}.Wcal.est.sBaverage(R.calGenuine(h.k));

h.KBdots(end+1)=plot(Kc,Bc,'ko');
h.KBtext(end+1)=text(Kc,Bc,'cal.',...
    'VerticalAlignment','top','HorizontalAlignment','left','color','b');
set(h.KBdots(end),'markerface','b')
end
function h=updateKBsymbols(h,trj,R)
% update KB plot symbols and text color to agree with current state
% classfication

% to do: also update RMS state lines

if(R.includetrj{h.k}(h.b))
    for k=1:trj{h.k}{h.b}.Wtrj.N
        s=R.trjstates{h.k}{h.b}(k);
        set(h.KBdots(k),'color','k',...%h.statesym{R.trjstates{h.k}{h.b}(k)}(2),...
            'markerface',h.statesym{s}(2),...
            'marker',h.statesym{s}(1));
        set(h.KBtext(k),'color',h.statesym{s}(2));
        set(h.trjRMS_line(k),'color',h.statesym{s}(2));
    end
    set(h.calRMS_line,'color','g');
    set([h.KBtext h.trjRMS_rms],'color','k');
    for nn=1:length(h.trjRMS_path)
        set(h.trjRMS_path(nn),'color',h.trjRMS_path_color{nn});
    end
    pause(0.1)
else
    delete(h.KBdots);
    h.KBdots=[];% indicate RMS for the different states    
    delete(h.KBtext);
    h.KBtext=[];
    
    set([h.calRMS_line h.trjRMS_line h.trjRMS_rms h.trjRMS_path],'color',get(h.KBpanel,'backgroundcolor'))
end

if(R.includetrj{h.k}(h.b)==false)
    set([h.XYaxes h.RMSaxes h.RMSfig],'color',get(h.KBpanel,'backgroundcolor'))
else
    set([h.XYaxes h.RMSaxes h.RMSfig],'color','w')
end
end
function h=plotNewRMS(h,trj,cal,R)
% replot the RMS trajectory
% plot the new RMS trajectory, remembering the y-limits

calRMS=cal{h.k}.Wcal.est.sRMS(R.calGenuine(h.k));
RMS=trj{h.k}{h.b}.Wtrj.est.sRMS;
T=length(h.trjRMS{h.k}{h.b});

figure(h.RMSfig)
ylim=get(gca,'ylim');
clf
hold on

% indicate calibration state RMS
h.calRMS_line=plot([1 T]/h.fSample,calRMS*[1 1],'g','linewidth',2);

% viterbi path (most likely path)
%h.trjRMS_path=plot((1:h.downSample:T)/h.fSample,RMS(h.viterbi{h.k}{h.b}),...
%    'color',h.trjRMS_path_color{1},'linewidth',2);
h.trjRMS_path=stairs(h.viterbi_t{h.k}{h.b}*h.downSample/h.fSample,RMS(h.viterbi_s{h.k}{h.b}),...
    'color','g','linewidth',2,'linestyle','-.');

% sequence of most likely states (different from viterbi path)
%h.trjRMS_path(2)=plot((1:h.downSample:T)/h.fSample,RMS(h.sMaxP{h.k}{h.b}),...
%    'color',h.trjRMS_path_color{2},'linewidth',1,'linestyle','--');
h.trjRMS_path(2)=stairs(h.sMaxP_t{h.k}{h.b}*h.downSample/h.fSample,RMS(h.sMaxP_s{h.k}{h.b}),...
    'color',h.trjRMS_path_color{1},'linewidth',2);

% plot RMS data
h.trjRMS_rms=plot((1:T)/h.fSample,h.trjRMS{h.k}{h.b},'k');

% indicate RMS for the different states
h.trjRMS_line=0*RMS;
for mm=1:length(RMS)
    h.trjRMS_line(mm)=plot([1 T]/h.fSample,RMS(mm)*[1 1],'k--');
    %set(h.trjRMS_line(mm),'linewidth',2)
end


h=toggleRMS_doit(h);
h.RMSaxes=gca;
set(gca,'ylim',ylim)
xlabel('t')
ylabel('RMS')
title('RMS plot of current trajectory')
legend('calibration RMS','Viterbi path','argmax p(s_t)','RMS data',2)
box on
end
%% panel update functions
function updateTrjText(h,cal,trj)
% update filenames and bead/trj numbers for the current trajectory
set(h.fileNameCAL,'string',['calibration: ' dirFileSplit(cal{h.k}.calDataFile)]);
set(h.fileNameTRJ,'string',['looping: ' dirFileSplit(trj{h.k}{h.b}.trjDataFile)]);
set(h.fileNameNUM,'string',['bead : ' int2str(h.k) '. subset : ' int2str(h.b) '.']);
end
function updateClassificationButtons(h,Wstates,includeTrj)
% set radio buttons in state panel according to states Wstates
for k=1:length(Wstates)
    set(h.SCstateButton(k,:),'value',0);
    set(h.SCstateButton(k,Wstates(k)),'value',1);
    set(h.SCstateGroup(k),'visible','on');
end
for k=1+length(Wstates):h.Nmax
    set(h.SCstateGroup(k),'visible','off');
end
% set exclude state indicator correctly
set(h.excludeTRJ,'value',~includeTrj);
% set scatter toggle buttins correctly
for ss=1:length(h.SCscatterToggle)
    set(h.SCscatterToggle(ss),'value',h.backgroundToggle(ss))
end

end
%% misc functions
function [gCal,trjstates]=defaultAnalysis(Wcal,Wtrj,RMS0)
% do a default state assignment based on RMS values and the default
% genuine/spurious criterion

[gCal,~,gTrj,~]=VB7_findGenuine(Wcal,Wtrj,RMS0);
trjstates=ones(1,Wtrj.N);   % start w sprurious trj states
trjstates(gTrj)=3;          % make all genune states looped

gRMS=Wtrj.est.sRMS(gTrj);
[~,gU]=max(gRMS);
trjstates(gTrj(gU))=2;      % the genuin state with largest RMS is unloped!
end
function R=doDefaultAnalysis(trj,cal,RMS0)
R.hastrj=cell(size(trj));   % indicator for existing analyzed trajectories
R.trjstates=cell(size(trj));% state number array
R.calGenuine=zeros(size(trj));
%calK=zeros(size(trj));
%calB=zeros(size(trj));
trjK=cell(size(trj));
trjB=cell(size(trj));
Nmax=1;
firstk=Inf;
firstb=Inf;
for k=1:length(trj)
    R.hastrj{k}   =true(1,length(trj{k}));
    R.trjstates{k}=cell(1,length(trj{k}));
    for b=1:length(trj{k})
        if(length(cal)>=k && ~isempty(cal{k}) && isfield(cal{k},'Wcal')&&~isempty(trj{k}{b}) && isfield(trj{k}{b},'Wtrj'))
            Wcal=cal{k}.Wcal;
            Wtrj=trj{k}{b}.Wtrj;
            Nmax=max(Nmax,Wtrj.N);
            if(k<=firstk && b<=firstb)
                firstk=k;firstb=b;
            end
        else % probably because something is not converged here
            R.hastrj{k}(b)=false; % remember that this trj is no good
            disp(['bead ' int2str(k) ', trj ' int2str(b) ' : no results found.'])
            continue
        end
        [R.calGenuine(k),R.trjstates{k}{b}]=defaultAnalysis(Wcal,Wtrj,RMS0);
        %calK(k)=Wcal.est.sKaverage(R.calGenuine(k));
        %calB(k)=Wcal.est.sBaverage(R.calGenuine(k));
        trjK{k}{b}=Wtrj.est.sKaverage;
        trjB{k}{b}=Wtrj.est.sBaverage;
    end
end
R.includetrj=R.hastrj;
end
function [fname,fdir]=dirFileSplit(ff)
% splits off the path to a file and its file name
ind=find(ff==filesep,1,'last');
if(isempty(ind))
    ind=0;
end
fdir=ff(1:ind-1);
fname=ff(ind+1:end);
if(isempty(fdir))
    fdir=['.' filesep];
end
end

