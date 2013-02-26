function [Kmin,Bmax]=VB7_KBupperleft(RMS0,B0)
% [Kmin,Bmax]=VB7_KBupperleft(RMS0,B0)
%
% find a minimum K and maximum B from solving the equations
% B = B0 * (1-K)
% B = 1 /RMS0^2/(1-K^2)
%
% Can me used as an extra selection criterion for genuine states.

% M.L. 2011-09-12   : started

Keq=@(k)((1-k).*(1-k.^2)-1/RMS0^2/B0);

Kmin=fzero(Keq,0.5);
Bmax=B0*(1-Kmin^2);


