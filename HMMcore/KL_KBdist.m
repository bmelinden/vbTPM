% KL=KL_KBdist(n,c,mu,v,n0,c0,mu0,v0)
%
% compute elementwise Kullback-Leibler divergence for vbTPM bead motion 
% distrobutions q(K,B).
%
% KL(n,c,mu,v|n0,c0,mu0,v0) = int q*log(q/q0) dKdB

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KL_KBdist.m, part of the vbTPM package
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
function KL=KL_KBdist(n,c,mu,v,n0,c0,mu0,v0)

KL=-(n+0.5)./c.*(c-c0-v0.*(mu-mu0).^2)...
   +0.5*log(v./v0)+(n0+0.5).*log(c./c0)...
   -gammaln(n+0.5)+gammaln(n0+0.5)...
   +(n-n0).*psi(n+0.5)+v0./v/2-0.5;

