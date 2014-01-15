% opt=VB7_getOptions(runinputfile)
%
% convert HMM runinput parameters from a runinput file into an options
% structure opt. In fact, all variables created by the command
% eval(runinputfile) are stored in the opt structure. 
%
% Additionally, the path part of runinputfile is added to the fields
% source_path and target_paths, so that the position of data and results
% are interpreted as given relative to the position of the runinput file.
%
% M.L. 2014-01-14

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_getOptions.m, access analysis options in the vbTPM package
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
function opt=VB7_getOptions(runinputfile)

% Split the runinput filename
[path_tmp, name_tmp, ext_tmp] = fileparts(runinputfile);
if(isempty(path_tmp))
    path_tmp='.';
    %else
    %warning('VB7_getOptions warning: runinput file not in the current folder,')
    % keep warning until bug is fixed!
end
if(ismember('.',name_tmp))
   error(['runinputfile : ' name_tmp '. It is currently not possible to handle runinputfile names containing a .'])
end

% read the raw options from the runinput file
oldFolder = cd(path_tmp);
% check that the runinput file is actually present in this folder
abspathfile=[oldFolder filesep path_tmp filesep name_tmp '.m'];
if(exist(abspathfile,'file'))
    eval(name_tmp)
else
    error(['VB7_getOptions: runinput file not found: ' abspathfile ' .'])    
end
cd(oldFolder);
clear oldFolder abspathfile; % forget what folder the options file happend to be called from
vv=whos;
opt=struct;

for m=1:length(vv)
    opt.(vv(m).name)=eval(vv(m).name);    
end
opt.localroot=pwd;

% modify target- and source paths to work relative to runinput file location
if(isfield(opt,'source_path'))
    opt.source_path=[path_tmp filesep opt.source_path];
else
    opt.source_path=[path_tmp filesep];
    warning(['VB7_getOptions: no source path specified in ' runinputfile '.'])
end
if(isfield(opt,'target_path'))
    opt.target_path=[path_tmp filesep opt.target_path];
else
    opt.target_path=[path_tmp filesep];
    warning(['VB7_getOptions: no target path specified in ' runinputfile '.'])
end
% add fileseparator at the end of the paths
if(~strcmp(opt.source_path(end),filesep))
    opt.source_path=[opt.source_path filesep];
end
if(~strcmp(opt.target_path(end),filesep))
    opt.target_path=[opt.target_path filesep];
end

% remove fields created in this file
opt=rmfield(opt,{'path_tmp','name_tmp','ext_tmp'});

% possibly check for missing fields and inserting fedault values

% make some sanity checks
