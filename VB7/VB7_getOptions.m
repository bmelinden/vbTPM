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
