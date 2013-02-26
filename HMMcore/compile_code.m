% Matlab script to compile all c files in this folder.
% M.L. 2011

mex -setup

ff=dir('*.c');

for k=1:length(ff)
    disp(ff(k).name)
    mex(ff(k).name,['-' computer('arch')],'-O')    
end