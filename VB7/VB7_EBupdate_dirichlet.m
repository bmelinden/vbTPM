% U=VB7_EBupdate_dirichlet(w,zpattern)
%
% Find Dirichlet distribution u that minimize the total KL divergences
% between each row of w and u, argmin_u sum_j KL_dirichlet(w(j,:)||u),
% i.e., each row of w is a Dirichlet distribution
%
% This means solving 
% psi(u) - psi(u0) = 1/M sum_j ( psi(w(j,:) - psi(w0(j)),
% with
% u0 = sum(u), w0(j)=sum(w(j,:))
%
% If present, zpattern specifices columns to ignore; only columns with
% zpattern>0 are included in the analysis. For columns with zpattern<=0,
% the corresponding entries of u are set to zero as well. By default, all
% columns are included.
% 
% ML 2012-05-02

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_EBupdate.m, empirical Bayes updates for Dirichlet distributions,
% in the vbTPM package
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
function U=VB7_EBupdate_dirichlet(w,zpattern)

opt=optimset('Display','off','TolFun',1e-15,'TolX',1e-15);
relTolF=1e-15;
maxIter=500;

N=size(w,2); % number of columns
M=size(w,1);
N0=N;
% ignore zero elements 
if(exist('zpattern','var')&& ~isempty(zpattern))
    ind=find(zpattern>0);
    N=length(ind);
else
    ind=1:N;
end
w=w(:,ind);
w0=sum(w,2);

% compute RHSs
wRHS=1/M*sum(psi(w)-psi(w0)*ones(1,N),1);

% solve it all in one go
lnU0=log(min(w,[],1)); % initial guess
lnFeq=@(x)(psi(exp(x))-psi(sum(exp(x)))-wRHS);
[lnU,~,Uflag,Uoutput]=fsolve(lnFeq,lnU0,opt);
U=exp(lnU);

if(Uflag<1)
    disp(['VB7_EBupdate_dirichelt fsolve output flag: ' int2str(Uflag)])
    Uoutput
    disp(['message : ' Uoutput.message])
    disp('--------------------')
end


%% alternative method: iteratively solve for each element of u
if(0) % does not seem to make a significant difference
    % initial guess: full numerical solution
    u=U;%min(w,[],1);
    Fu=0;
    for m=1:M
        Fu=Fu-KL_dirichlet(w(m,:),u);
    end
    
    for iter=1:maxIter
        u_old=u;
        
        for j=1:N
            u0j=sum(u([1:j-1 j+1:end]));
            lnfeq=@(lnuj)(psi(exp(lnuj))-psi(exp(lnuj)+u0j)-wRHS(j)); %scalar equation
            lnuj=fsolve(lnfeq,log(u(j)),opt);
            
            u(j)=exp(lnuj);
        end
        Fu_old=Fu;
        Fu=0;
        for m=1:M
            Fu=Fu-KL_dirichlet(w(m,:),u);
        end
        dFrel=abs((Fu-Fu_old)/Fu);
        %disp(num2str([dFrel u-u_old]))
        if(dFrel<relTolF)
            disp(['EB dirichlet converged after ' int2str(iter) ' iterations'])
            break
        end
    end
    if(iter>=maxIter)
        warning('EB_dirichlet reached maximum number of iterations')
    end
    disp(num2str(U./u-1))
    
    % reinsert columns corresponding to zero elements
    if(N0~=N)
        u2=zeros(1,N0);
        u2(ind)=u;
        u=u2;
    end
end

% reinsert columns corresponding to zero elements
if(N0~=N)
    u2=zeros(1,N0);
    u2(ind)=U;
    U=u2;
end

