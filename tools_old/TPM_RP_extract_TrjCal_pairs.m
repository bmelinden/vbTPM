% [xTrj,xCal,xcorrTrj,xcorrCal,trjName,RMStrj]=...
%           TPM_RP_extract_TrjCal_pairs(prefix,corr,raw,visualcheck)
% 
% reads nolacdata.mat and lacdata.mat from ocal directiry, and outputs
% coordinate time series of all calibration and looping trajectories, one
% trajectory per file. The selection and cutting of time series is read
% from the analysis structure, so that the end result can be compared
% exactly with the results from histogram analysis.
%
% prefix is a string to identify this data sets from others, e.g.
% 'PUC306_100fMSJLacI' leads to filenames of the form
% PUC306_100fMSJLacI_bead41_trjA_corr.mat etc., instead of just 
% bead41_trjA_corr.mat etc., which is the default if no prefix is given.
% (an underscore is inserted automatically)
%
% corr,raw are flags to indicate if driftcorrected and raw coordinates are
% to be extracted respectively.
%
% visualcheck (optional) is a boolean that, if true, draws the raw
% coordinates and lets the user inspect the cuts that are in addition to
% the ones based on previous RMS analysis, as a safety check.

%% change-log
% M.L. 2011-02-03   : added possibility to give filename prefixes
% M.L. 2011-02-22   : keep track of indices to trajectories in original
%                     lacdata.mat files to make it easier to backtrack some
%                     trajectory.

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TPM_RP_extract_TrjCal_pairs.m, part of the vbTPM package
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
function [xTrj,xCal,xcorrTrj,xcorrCal,trjName,RMStrj]=TPM_RP_extract_TrjCal_pairs(prefix,corr,raw,visualcheck)

%% deal with input
source=pwd;
if(~strcmp(source(end),filesep))
    source=[source filesep];
end
if(~exist('prefix','var')|| isempty(prefix))
    prefix='';
end
if(~isempty(prefix) && ~strcmp(prefix(end),'_'))
    prefix=[prefix '_'];
end
if(~exist('visualcheck','var')|| isempty(visualcheck))
    visualcheck=false;
end
if(~exist('corr','var')|| isempty(corr))
    corr=false;
end
if(~exist('raw','var')|| isempty(raw))
    raw=false;
end
beaddir=[prefix 'XYdata'];

%% code

% load file structures
tic
disp('loading data...')
C=load([source 'nolacdata.mat']);
D=load([source 'lacdata.mat']);
toc

disp('concatenating data...')
tic
lacx=TPM_RP_concatenate(D.lacX,D.lacY);
laccorrx=TPM_RP_concatenate(D.laccorrX,D.laccorrY);
nolacx=TPM_RP_concatenate(C.nolacX,C.nolacY);
nolaccorrx=TPM_RP_concatenate(C.nolaccorrX,C.nolaccorrY);
toc

NF=length(lacx);   % number of coordinate batches
NC=length(nolacx); % number of calibration batches
if(NF~=2*ceil(NF/2)) % then N is an even integer
    C
    D
    error('Does not know how to handle odd number of lacX coordinate fields')
end
if(NF==NC)
   two_looping_rounds=false;
elseif(NF==2*NC)
   two_looping_rounds=true;
   disp('Found two looping trajectories per calibration')
else
    error('NF ~= NC or 2*NC. Cannot handle this case...')
end

disp('extracting analyzed trajectories...')
tic
xTrj=cell(0);
xCal=cell(0);
xcorrTrj=cell(0);
xcorrCal=cell(0);
trjName=cell(0);
RMStrj=cell(0);
trjRecord=cell(0);
calRecord=cell(0);
nBead=0;
nTrj=0;
tMaxRMS=[];
tMax=[];
disp('found the following combinations:')
for k1=1:NC
    if(two_looping_rounds)
        k2=k1+NF/2;
        disp([C.nolacnames{k1} ', ' D.alllacnames{k1} ', ' D.alllacnames{k2}])
    else
        disp([C.nolacnames{k1} ', ' D.alllacnames{k1}])        
    end
end
for k1=1:NC
    
    mRMS=0;
    for m=1:length(D.record{k1})
       if(D.record{k1}(m)==1) % then this bead was recorded for further analysis
        nBead=nBead+1;
        nTrj=nTrj+1;
        mRMS=mRMS+1;
        
        % calibration curve
        trjName{nBead,1}=C.nolacnames{k1};
        xCal{nBead,1}=nolacx{k1}{m};
        xcorrCal{nBead,1}=nolaccorrx{k1}{m};
        calRecord{nBead,1}=[k1 m];
        
        % looping trajectory        
        trjName{nBead,2}=D.alllacnames{k1};
        
        % cut according to previous RMS analysis
        tMax1=find(D.lacANAL_RMS{k1}(mRMS,:)~=0,1,'last');
        xTrj{nBead,1}=lacx{k1}{m}(1:tMax1,:);
        xcorrTrj{nBead,1}=laccorrx{k1}{m}(1:tMax1,:);  
        RMStrj{nBead,1}=D.lacANAL_RMS{k1}(mRMS,1:tMax1);
        trjRecord{nBead,1}=[k1 m];
        
        % cut for HMM analysis: remove portions where coordinates do not
        % change at all.         
        tMAxXY=min(find(xcorrTrj{nBead,1}(1:end-1,1)==xcorrTrj{nBead,1}(2:end,1),1),...
                   find(xcorrTrj{nBead,1}(1:end-1,2)==xcorrTrj{nBead,1}(2:end,1),1));
        tMax2=min([tMax1 tMAxXY]);
        
        if(visualcheck)
            figure(1),clf,hold on
            plot(xcorrTrj{nBead,1}(:,1),'k')
            plot(xcorrTrj{nBead,1}(:,2)-100,'b')
            plot(xcorrTrj{nBead,1}(1:tMax2,1)+500,'r')
            plot(xcorrTrj{nBead,1}(1:tMax2,2)+400,'g')
            pause(0.2)
            if(tMax1 ~= tMax2)
                pause
            end
        end
        
        xTrj{nBead,1}=lacx{k1}{m}(1:tMax2,:);
        xcorrTrj{nBead,1}=laccorrx{k1}{m}(1:tMax2,:);  
        RMStrj{nBead,1}=D.lacANAL_RMS{k1}(mRMS,1:tMax2);

        tMax(nBead,1)=tMax2;
        tMaxRMS(nBead,1)=tMax1;
        
        
       end
       
       if(two_looping_rounds)
           k2=k1+NF/2;
           mRMS2=0;
           if(length(D.record)>=k2 && ~isempty(D.record{k2}) && D.record{k2}(m)==1)
               % then this bead was recorded twice
               mRMS2=mRMS2+1;
               nTrj=nTrj+1;
               
               trjName{nBead,3}=D.alllacnames{k2};
               tMax1=find(D.lacANAL_RMS{k2}(mRMS2,:)~=0,1,'last'); 
               xTrj{nBead,2}=lacx{k2}{m}(1:tMax1,:);
               xcorrTrj{nBead,2}=laccorrx{k2}{m}(1:tMax1,:);
               RMStrj{nBead,2}=D.lacANAL_RMS{k2}(mRMS2,1:tMax1);
               trjRecord{nBead,2}=[k2 m];

               % cut for HMM analysis: remove portions where coordinates do not
               % change at all.
               tMAxXY=min(find(xcorrTrj{nBead,2}(1:end-1,1)==xcorrTrj{nBead,2}(2:end,1),1),...
                          find(xcorrTrj{nBead,2}(1:end-1,2)==xcorrTrj{nBead,2}(2:end,2),1));
               tMax2=min([tMax1 tMAxXY]);
               
               if(visualcheck)
                   figure(1),clf,hold on
                   plot(xcorrTrj{nBead,2}(:,1),'k')
                   plot(xcorrTrj{nBead,2}(:,2)-100,'b')
                   plot(xcorrTrj{nBead,2}(1:tMax2,1)+500,'r')
                   plot(xcorrTrj{nBead,2}(1:tMax2,2)+400,'g')
                   pause(0.2)
                   if(tMax1 ~= tMax2)
                       pause
                   end
               end
               
               xTrj{nBead,2}=lacx{k2}{m}(1:tMax2,:);
               xcorrTrj{nBead,2}=laccorrx{k2}{m}(1:tMax2,:);
               RMStrj{nBead,2}=D.lacANAL_RMS{k2}(mRMS2,1:tMax2);               

               tMax(nBead,2)=tMax2;
               tMaxRMS(nBead,2)=tMax1;
           end
       end              
    end
end
disp(['... found ' int2str(nTrj) ' trajectories from ' int2str(nBead) ' beads'])
toc

disp('writing files...')
mkdir(beaddir);
nDig=int2str(ceil(log10(nBead))); % number format, e.g., bead01* or bead0001*
for k=1:nBead

    tCutRMS=tMaxRMS(nBead,1);tCut=tMax(nBead,1);
    
    if(raw)
        file   =sprintf(['%sbead%0' nDig 'd_cal_raw.mat'],prefix  ,k);
        x=xCal{k,1};
        name=[trjName{k,1} ', raw bead positions'];        
        record=calRecord{k,1};
        save([beaddir filesep file],'x','name','tCutRMS','tCut','record')
        disp(['wrote : ' beaddir filesep file])
        
        file =sprintf(['%sbead%0' nDig 'd_trjA_raw.mat'],prefix,k);
        x=xTrj{k,1};
        name=trjName{k,2};
        record=trjRecord{k,1};
        save([beaddir filesep file],'x','name','tCutRMS','tCut','record')
        disp(['wrote : ' beaddir filesep file])
    end
    if(corr)
        file  =sprintf(['%sbead%0' nDig 'd_cal_corr.mat'],prefix,k);
        x=xcorrCal{k,1};
        name=[trjName{k,1} ', drift-corrected bead positions.'];
        record=calRecord{k,1};
        save([beaddir filesep file],'x','name','tCutRMS','tCut','record')
        disp(['wrote : ' beaddir filesep file])
        
        file=sprintf(['%sbead%0' nDig 'd_trjA_corr.mat'],prefix,k);
        x=xcorrTrj{k,1};
        name=[trjName{k,2} ', drift-corrected bead positions.'];
        record=trjRecord{k,2};
        save([beaddir filesep file],'x','name','tCutRMS','tCut','record')
        disp(['wrote : ' beaddir filesep file])
    end
    if(two_looping_rounds)
        if(size(xTrj,2)>=2 && ~isempty(xTrj{k,2}))
            tCutRMS=tMaxRMS(nBead,2);tCut=tMax(nBead,2);
            
            if(raw)
                file =sprintf(['%sbead%0' nDig 'd_trjB_raw.mat'],prefix,k);
                x=xTrj{k,2};
                name=[trjName{k,3} ', raw bead positions, second trj.'];
                record=trjRecord{k,2};
                save([beaddir filesep file],'x','name','tCutRMS','tCut','record')
                disp(['wrote : ' beaddir filesep file])
            end
            if(corr)
                file=sprintf(['%sbead%0' nDig 'd_trjB_corr.mat'],prefix,k);
                x=xcorrTrj{k,2};
                name=[trjName{k,3} ', drift-corrected bead positions, second trj.'];
                record=trjRecord{k,2};
                save([beaddir filesep file],'x','name','tCutRMS','tCut','record')
                disp(['wrote : ' beaddir filesep file])
            end
        end
    end
    
end
disp('-----------------------')
disp('file name cell vectors (always named corr, even if raw was written):')
disp('calibration_filename={...')
for k=1:nBead
    disp([char(39) sprintf(['%sbead%0' nDig 'd_cal_corr.mat'],prefix,k) char(39) ',...'])
end
disp('};')
disp(' ')
disp('looping_filename={...')
for k=1:nBead
    str=['{' char(39) ...
        sprintf(['%sbead%0' nDig 'd_trjA_corr.mat'],prefix,k)  ...
        char(39)];
    if(two_looping_rounds && size(xTrj,2)>=2 && ~isempty(xTrj{k,2}))
        str=[str ', ' char(39) ...
            sprintf(['%sbead%0' nDig 'd_trjB_corr.mat'],prefix,k) ...
            char(39) '},...'];
    else
        str=[str '},...'];
    end
    disp(str)
end
disp('};')

