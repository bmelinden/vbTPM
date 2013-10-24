%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_VBEMiter_nomex.m, VBEM iteration without mex files, in the vbTPM package
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
function opt=VB7_getOptions(runinputfile)
% opt=VB7_getOptions(runinputfile)
%
% convert HMM runinput parameters from a runinput file into an options
% structure opt. In fact, all variables created by the command
% eval(runinputfile) are stored in the opt structure. 
% An (incomplete) sanity check of some parameter values is also performed.
%
% M.L. 2013-07-07

%% start of actual code
% Split the runinput filename
[path, name, ext] = fileparts(runinputfile);
if(isempty(path))
    path='.';
else
    warning('VB7_getOptions warning: runinput file not in the current folder,')
end
if(ismember('.',name))
   error(['runinputfile : ' name '. It is currently not possible to handle runinputfile names containing a .'])
end

% read the raw options from the runinput file
oldFolder = cd(path);
eval(name)
cd(oldFolder);
clear oldFolder; % forget what folder the options file happend to be called from
vv=whos;
opt=struct;

for m=1:length(vv)
    opt.(vv(m).name)=eval(vv(m).name);    
end
opt.localroot=pwd;

% possibly check for missing fields and inserting fedault values

% make some sanity checks
