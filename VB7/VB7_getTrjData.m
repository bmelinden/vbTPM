function [trjdat,caldat]=VB7_getTrjData(runinputfile,k,b)
% [trjdat,caldat]=VB7_getTrjData(runinputfile,k,b)
%
% retrieve looping and calibration data trajectory k,b for the analysis
% specified in runinput file. The objects trjdat,caldat have fields
%
% name      : name of data file
% data      : preprocessed data object, as used in the analysis
% x         : position trajectory
% RMS       : RMS trace
%
% M.L. 2012-02-23

% change-log
% ML 2012-06-01 : handle runinput filenames ending with .m
% ML 2013-05-11 : handle runinput files in different folders

% get settings for this runinput file
if(strcmp(runinputfile(end-1:end),'.m'))
    runinputfile=runinputfile(1:end-2);
end
[VB7tmp_path,VB7tmp_file]=fileparts(runinputfile);
VB7tmp_dir0=pwd;
cd(VB7tmp_path)
eval(VB7tmp_file);
cd(VB7tmp_dir0)
clear VB7tmp_path VB7tmp_file VB7tmp_dir0

% determine data file names and load data
caldat.name=[source_path calibration_filename{k}];
trjdat.name=[source_path looping_filename{k}{b}];
if(~strcmp(source_path(end),filesep))
    caldat.name=[source_path filesep calibration_filename{k}];
    trjdat.name=[source_path filesep looping_filename{k}{b}];
end
caldat=load(caldat.name);
trjdat=load(trjdat.name);

% rename the coordinate field to 'x'
if(~strcmp('x',calibration_xyfield))
    caldat.x=caldat.(calibration_xyfield);
    caldat=rmfield(caldat,calibration_xyfield);
end
if(~strcmp('x',looping_xyfield))
    trjdat.x=trj.(looping_xyfield);
    trjdat=rmfield(trjdat,looping_xyfield);
end

% driftcorrect by lowpass fileter, if specified in the runinputfile
if(driftcorrection)
    caldat.x=BWdriftcorrect(caldat.x,fCut,fSample);
    trjdat.x=BWdriftcorrect(trjdat.x,fCut,fSample);
    caldat.fCut=fCut;
    trjdat.fCut=fCut;
end

% preprocessed data
caldat.data=VB7_preprocess(caldat.x,downSample,fSample);
trjdat.data=VB7_preprocess(trjdat.x,downSample,fSample);

% RMS trace
caldat.RMS=RMSKBgaussfilter(caldat.x,tSigma,fSample);
trjdat.RMS=RMSKBgaussfilter(trjdat.x,tSigma,fSample);
caldat.tSigma=tSigma;
trjdat.tSigma=tSigma;

% useful information
trjdat.fSample=fSample;
trjdat.downSample=downSample;
caldat.fSample=fSample;
caldat.downSample=downSample;


