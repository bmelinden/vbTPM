% data=VB7_preprocess(x,d,fSample,comment)
% 
% Preprocessing step including downsampling for the VB3-VB5,VB7 family of analysis
% code. 
% x is the position trajectory
% d is the downsampling factor (default 1= do downsampling)
% fSample is the native sampling frequency (only for later use). Default 30
% Hz.
% comment: optional string to attach to the data. Good for debugging.

%% change-log
% M.L. 2012-01-12   : added an optimized loop for the case d=1 (no
%                     downsampling).
% M.L. 2011-03-15   : added function to convert from downsampled to
%                     original time indices. Removed x1, x12t, since not
%                     used by the VB7 algorithm.
% M.L. 2011-02-03   : changed name to VB7_preprocess (same code as VB3)
% M.L. 2010-10-26   : saving the native sampling frequency as well
% M.L. 2010-10-20   : changed to averaging each block, to get more
%                     homogeneous expressions in the vbEM steps (basically
%                     multiply everything by the downSampling factor).
% M.L. 2010-10-19	: started

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VB7_preprocess.m, preprocessor for vbTPM 
% =========================================================================
% 
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
function data=VB7_preprocess(x,d,fSample,comment)

%% parameter check
if(~exist('d','var')|| isempty('d')); d=1; end
if(~exist('fSample','var')|| isempty('fSample')); fSample=30; end
if( d ~= round(d) || d<1 )
    disp(['d = ' num2str(d)])
    error('VB7_preprocess: downsampling factor d must be a positive integer')
end
[a,b]=size(x);
if(a<=2 || b~=2)
    if( a==2 && b>2) % needs transposing
        warning(['VB7_preprocess: size(x) = [' int2str(size(x)) ']. Column vector preferred.'])
        x=x';
    else % wrong format
        error(['VB7_preprocess: size(x) = [' int2str(size(x)) ']. Wrong data size.'])
    end
end
clear a b;
if(exist('comment','var') && ~isempty(comment)); data.comment=comment;end

%% start of actual code
data.downSampling=d;
data.fSample0=fSample;
%data.x1=x(1,:);
%data.x12=norm(data.x1)^2;
% truncate dta at the end, to make all blocks length d
L=1+floor((size(x,1)-1)/d); % number of full blocks
data.R2  =zeros(L,1); % x(t)^2
data.R2m1=zeros(L,1); % x(t-1)^2
data.X12=zeros(L,1); % x(t)x(t-1)
b1=2;
b2=b1+d-1;

data.tind(1,1)=1;
if(d==1) % then b2=b1, and we can avoid coding the innermost loop
    data=data;
    data.R2(2:end,1)  =sum(x(2:end,:).^2,2);
    data.R2m1(2:end,1)=sum(x(1:end-1,:).^2,2);
    data.X12(2:end,1) =sum(x(2:end,:).*x(1:end-1,:),2);
else % downsmapled version (slower)
    for k=2:L
        data.R2(k,1)  =sum(sum(x(b1:b2,:).^2    ,2))/d;
        data.R2m1(k,1)=sum(sum(x(b1-1:b2-1,:).^2,2))/d;
        data.X12(k,1)=sum(sum(x(b1:b2,:).*x(b1-1:b2-1,:),2))/d;
        b1=b1+d;
        b2=b2+d;
    end
end

% upsampling function
data.fUpSample=@upSampleIndices;
    function ti=upSampleIndices(ind)
        ti=zeros(1,length(ind)*d);
        n1=1;        
        for k=1:length(ind);
            if(ind(k)==1)
                n2=n1;
                ti(n1)=1;
            else
                n2=n1+d-1;
                ti(n1:n2)=(2+(ind(k)-2)*d):(2+(ind(k)-1)*d-1);
            end
            n1=n2+1;
        end
        ti=ti(1:n2);
end
end
