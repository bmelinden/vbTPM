function [h,K,B,RMS,Krange,Brange,KBn]=KBhistogram(X,tSigma,fSample,KBrange,NKB,fig)
% [h,K,B,RMS,Krange,Brange,KBn]=KBhistogram(X,tSigma,fSample,KBrange,NKB,fig)
%
% h : plot handle
% K,B,RMS . gaussian filtered estimates
% KBrange =[Kmin Kmax Bmin Bmax] range for the histogram, which is plotted
% in figure fig (or axes(fig)), if the argument is given.
% NBK = [NK NB] is teh number of lattice points in each dimension of the
% histogram
% 
% M.L. 2010-10-17

%% change log
% M.L. 2011-02-02 : improved help text

%% code
if(~exist('tSigma','var')|| isempty(tSigma));  tSigma=4;end
if(~exist('fSample','var')|| isempty(fSample));fSample=30;end
[RMS,K,B]=RMSKBgaussFilter(X,tSigma,fSample);

if(~exist('NKB','var')|| isempty(NKB));NKB=[100 100];end
if(length(NKB)==1); NKB=NKB*[1 1]; end
if(~exist('KBrange','var')|| isempty(KBrange));KBrange=[-1 1 0 max(B) ];end
   
Krange=linspace(KBrange(1),KBrange(2),NKB(1));
Brange=linspace(KBrange(3),KBrange(4),NKB(2));

[a,b]=hist3([K B],{Krange, Brange});
KBn=a';

if(exist('fig','var') && ~isempty(fig))
    try
        figure(fig);
    catch
        axes(fig);
    end
    h=pcolor(Krange,Brange,KBn/sum(sum(KBn)));
    axis(KBrange);
    ylabel('B');
    xlabel('K');
    set(h,'edgecol','non')
else
    h=[];
end

