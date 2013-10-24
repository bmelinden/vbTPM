function VB7_printGPL(name)
disp(' ');
disp('=========================================================================');
if(exist('name','var'))
    disp([name ', part of vbTPM.'])
    disp(' ');
end
disp('vbTPM is a package for variational Bayes analysis of Tethered Particle ');
disp('Motion data. See http://sourceforge.net/projects/vbtpm/');
disp(' ' );
disp('Copyright (C) 2013 Martin Lindén, E-mail: bmelinden@gmail.com');
disp('=========================================================================');
disp(' ');
