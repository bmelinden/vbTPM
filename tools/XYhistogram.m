function [h,Xrange,Yrange,XYn,dx,dy]=XYhistogram(X,Y,XYbins,XYrange,fig)
% [h,Xrange,Yrange,XYn,dx,dy]=XYhistogram(X,Y,XYbins,XYrange,fig)
% h : plot handle
% Xrange,Yrange : bin center vectors for 2D histograms
% XYn           : histogram counts
% dx,dy         : bin widths
% 
% X,Y   : indata series
% XYbins=[nx ny]: number of bins in each direction (default 100).
% XYrange = [xmin xmax, ymin ymax] : range of histogram (default max/min of
%                                    X,Y in each direction).
% M.L. 2010-10-17

%% deal with parameters
if(~exist('XYbins','var')|| isempty(XYbins));XYbins=[100 100];end
if(length(XYbins)==1); XYbins=XYbins*[1 1]; end
if(~exist('XYrange','var')|| isempty(XYrange));XYrange=[min(X) max(X) min(Y) max(Y) ];end
X=reshape(X,length(X),1);
Y=reshape(Y,length(Y),1);
%% actual function   
Xrange=linspace(XYrange(1),XYrange(2),XYbins(1));
Yrange=linspace(XYrange(3),XYrange(4),XYbins(2));
[a,b]=hist3([X Y],{Xrange, Yrange});
XYn=a';
dx=diff(Xrange(1:2));
dy=diff(Yrange(1:2));

if(exist('fig','var') && ~isempty(fig))
    try
        figure(fig);
    catch
        axes(fig);
    end
    h=pcolor(Xrange,Yrange,XYn/sum(sum(XYn)));
    axis(XYrange);
    set(h,'edgecol','non')
else
    h=[];
end
