% script to go through all VB7 result files in this directory and make
% them a lot smaller by removing the reference to the internal function
% W.minimalStorage from the VB7 data structure. Mostly useful for rescuing 
% old analysis results.

 % change this line to a pattern suitable for 
f=dir('*.mat');
for k=1:length(f)
    disp(f(k).name)
    a=load(f(k).name);
    if(isfield(a,'Wtrj'))
        a.Wtrj=VB7_minimalStorage(a.Wtrj);
        save(f(k).name,'-struct','a')
    end
end
