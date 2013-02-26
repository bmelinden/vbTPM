function [h,x,y,z]=PcolHHistogram(K,B,KBrange,NKB,fig)


if(~exist('NKB','var')|| isempty(NKB));NKB=[100 100];end
if(length(NKB)==1); NKB=NKB*[1 1]; end
if(~exist('KBrange','var')|| isempty(KBrange));KBrange=[-1 1 0 max(B) ];end
   


if(exist('fig','var'))
    try
        figure(fig);
    catch
        axes(fig);
    end
end

Krange=linspace(KBrange(1),KBrange(2),NKB(1));
Brange=linspace(KBrange(3),KBrange(4),NKB(2));

[a,b]=hist3([K B],{Krange, Brange});
x=b{1};
y=b{2};
z=a'/sum(sum(a));
h=pcolor(x,y,z);
axis(KBrange);
set(h,'edgecol','non')
