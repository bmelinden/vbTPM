% Add paths that are needed to run the vbTPM analysis 

dir0=pwd;
addpath(genpath([dir0 filesep '.' filesep 'VB7']))
addpath(genpath([dir0 filesep '.' filesep 'HMMcore']))
addpath(genpath([dir0 filesep '.' filesep 'tools']))
addpath([dir0 filesep '.'])
clear dir0
VB7_printGPL();
disp('Added local vbTPM paths.')

