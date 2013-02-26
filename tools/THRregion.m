function [x,y,a]=THRregion(thr,Xrange,Yrange,Z)
% [x,y,a]=THRregion(thr,Xrange,Yrange,Z,dn);
%
% Thresholded selection from 2D histogram. THRregion finds all (x,y)-values
% such that Z(x,y)> a, where a is a threshold chosen such that 
% sum(Z(Z>a)) / sum(Z) ~ thr
% Xrange, Yrange are vectors of bin centers in the two dimensions, Z is an
% intensity matrix for each bin.
%  
% M.L. 2010-10-19

%% change-log
% M.L. 2011-02-04   : added linear term to f(a) to avoid flat regions.

if(thr<0 || thr >1)
    warning('THRregion called with thr outside [0,1],')
end
aMin=min(min(Z));
aMax=max(max(Z));

fMax=max(abs([thr 1-thr]))^2;
daf=0.01*fMax/(aMax-aMin);
f=@(a)((thr-sum(sum(Z(Z>a)))/sum(sum(Z)))^2+daf*(a-aMin));
a=fminbnd(f,aMin,aMax);

dx=Xrange(2)-Xrange(1);
dy=Yrange(2)-Yrange(1);
[x,y]=meshgrid(Xrange,Yrange);
x=x(Z>a)+dx/2;
y=y(Z>a)+dy/2;

