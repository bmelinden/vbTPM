% est=VB7_Mestimates(M,fSample,downSampling)
% Compute some estimates based on variational distribution
% parametrers (no hyper-parameters substracted).
%
% ML 2012-12-06

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_Mestimates.m, generate parameter estimates from variational 
% distributions in the vbTPM package.
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
function est=VB7_Mestimates(M,fSample,downSampling)

% occupation probabilities
est.Pocc=sum(M.wA,2)/sum(sum(M.wA,2));
% transition rates and dwell times: s(t)
est.A=rowNormalize(M.wA);
est.dA=zeros(size(M.wA));
a0=sum(M.wA,2);

for k=1:length(M.mu)
    est.dA(k,:)=sqrt(est.A(k,:).*(1-est.A(k,:))/(1+a0(k)));
end
est.tD=1./(1-diag(est.A))/(fSample/downSampling);

% emission parameters
est.Kaverage=M.mu;
est.Baverage=(M.n+1/2)./M.c;
est.Kstd=sqrt(M.c/2./M.v./(M.n-0.5)); % sqrt(Var(K))
est.Bstd=sqrt(M.n+1/2)./M.c;            % sqrt(Var(B))
est.KBcov=est.Kaverage.*est.Baverage; % <K*B>

est.RMS=sqrt(1./est.Baverage./(1-est.Kaverage.^2));
est.TC=-1./log(M.mu)/(fSample); % real time units; correction for downsampling not needed
est.NC=-1./log(M.mu);                              % sampling time units

