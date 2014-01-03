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
%
% updated to deal with new data format as of Feb 2011. M.L.

%% change-log
% M.L. 2011-02-03   : added possibility to give filename prefixes
% M.L. 2011-02-22   : keep track of indices to trajectories in original
%                     lacdata.mat files to make it easier to backtrack some
%                     trajectory.
% M.L. 2011-02-22   : dealing with new data format, using cutpts fields to
%                     judge how much to save

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TPM_RP_extract_TrjCal_pairs_v2.m, part of the vbTPM package
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
lacx=cell(1,length(D.lacX));
laccorrx=cell(1,length(D.lacX));
for k=1:length(D.lacX)
    for b=1:size(D.lacX{k},1);
        tCut=D.cutpts{k}(b);
   
        lacx{k}{b}=[D.lacX{k}(b,1:tCut)' ...
                    D.lacY{k}(b,1:tCut)'];
        laccorrx{k}{b}=[D.laccorrX{k}(b,1:tCut)' ...
                        D.laccorrY{k}(b,1:tCut)'];
    end
end
nolacx=cell(1,length(C.nolacX));
nolaccorrx=cell(1,length(C.nolacX));
for k=1:length(C.nolacX)
    for b=1:size(C.nolacX{k},1);
        nolacx{k}{b}    =[C.nolacX{k}(b,:)'     C.nolacY{k}(b,:)'];
        nolaccorrx{k}{b}=[C.nolaccorrX{k}(b,:)' C.nolaccorrY{k}(b,:)'];
    end
end
toc

NF=length(lacx);   % number of coordinate batches
NC=length(nolacx); % number of calibration batches
disp(['Found ' int2str(NC) ' calibration rounds, and ' ...
    int2str(NF) ' looping rounds.'])
if(NF~=NC) 
    C
    D
    error('Does not know how to handle different of calibration and coordinate batches. Please improve the code...')
end

disp('writing to files...')
nDig='2';
mkdir(beaddir);
rcfile=cell(0);
rlfile=cell(0);
ccfile=cell(0);
clfile=cell(0);
for k=1:NF
    for b=1:size(D.lacX{k},1);
        if(~isempty(lacx{k}{b}))
            
            if(visualcheck)
                figure(1),clf,hold on
                plot(lacx{k}{b}(:,1)+100,'k')
                plot(lacx{k}{b}(:,2)-100,'r')
                plot(laccorrx{k}{b}(:,1)+100,'b')
                plot(laccorrx{k}{b}(:,2)-100,'g')
                
                pause(0.5)
            end
            
            if(raw)
                rcfile{k,b}   =sprintf(['%sbead%0' nDig 'darea%0' nDig 'd_cal_raw.mat'],prefix,k,b);
                x= nolacx{k}{b};
                name=[C.nolacnames{k} ', raw calibration trace'];
                index=[k b];
                save([beaddir filesep rcfile{k,b}],'x','name','index')
                disp(['wrote : ' beaddir filesep rcfile{k,b}])
                
                rlfile{k,b} =sprintf(['%sbead%0' nDig 'darea%0' nDig 'd_trj_raw.mat'],prefix,k,b);
                x=lacx{k}{b};
                name=[C.nolacnames{k} ', raw looping trace'];
                index=[k b];
                save([beaddir filesep rlfile{k,b}],'x','name','index')
                disp(['wrote : ' beaddir filesep rlfile{k,b}])
            end
           if(corr)
                ccfile{k,b} = sprintf(['%sbead%0' nDig 'darea%0' nDig 'd_cal_corr.mat'],prefix,k,b);
                x    = nolaccorrx{k}{b};
                name = [C.nolacnames{k} ', drift-corrected calibration trace'];
                index= [k b];
                save([beaddir filesep ccfile{k,b}],'x','name','index')
                disp(['wrote : ' beaddir filesep ccfile{k,b}])
                
                clfile =sprintf(['%sbead%0' nDig 'darea%0' nDig 'd_trj_corr.mat'],prefix,k,b);
                x=laccorrx{k}{b};
                name=[C.nolacnames{k} ', drift-corrected looping trace'];
                index=[k b];
                save([beaddir filesep clfile{k,b}],'x','name','index')
                disp(['wrote : ' beaddir filesep clfile{k,b}])
            end
            
            
        end
    end
end
disp('-----------------------')
disp('file name cell vectors:')
disp('calibration_filename={...')

if(raw)
    for k=1:NF
        for b=1:size(D.lacX{k},1);
            if(~isempty(lacx{k}{b}))
                disp([char(39) rcfile{k,b} char(39) ',...'])
            end
        end
    end
    disp('};')
    disp(' ')
    disp('looping_filename={...')
    for k=1:NF
        for b=1:size(D.lacX{k},1);
            if(~isempty(lacx{k}{b}))
                disp([ '{' char(39) rlfile{k,b} char(39) '},...'])
            end
        end
    end
    disp('};')
    disp(' ')
end
if(corr)
    for k=1:NF
        for b=1:size(D.lacX{k},1);
            if(~isempty(lacx{k}{b}))
                disp([char(39) ccfile{k,b} char(39) ',...'])
            end
        end
    end
    disp('};')
    disp(' ')
    disp('looping_filename={...')
    for k=1:NF
        for b=1:size(D.lacX{k},1);
            if(~isempty(lacx{k}{b}))
                disp([char(39) clfile{k,b} char(39) ',...'])
            end
        end
    end
    disp('};')
    disp(' ')
end

return

