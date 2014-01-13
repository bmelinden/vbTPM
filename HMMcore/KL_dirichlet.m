% KL=KL_dirichlet(w,u)
% 
% compute the Kullback-Leibler divergence of two Dirichlet distributions
% with parameters w and u,
% KL = KL(w||u) = int p(x|w) ln ( p(x|w) / p(x|u) ) dx
% 
% In addition, elements with u(j)=0 are ignored, equivalent to 
% u=u(find(u~=0));w=w(find(u~=0));
%
% (In the VBHMM setting, w is the variational distribution, and u is the
% prior.)
%
% M.L. 2012-05-02

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KL_dirichlet.m, part of the vbTPM package
% =========================================================================
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
function KL=KL_dirichlet(w,u)

% find non-zero indices
ind=find(u~=0);
u=u(ind);
w=w(ind);

w0=sum(w);
u0=sum(u);

KL=gammaln(w0)-gammaln(u0)...
    +sum(gammaln(w)-gammaln(u)...
    +(w-u).*(psi(w)-psi(w0)));





