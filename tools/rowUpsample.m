function V=rowUpsample(v,ds,Tmax)
% upsamples the rows of the matrix v by a factor ds, and then truncates it
% to Tmax rows. The first row is not upsampled. This is meant to make the
% results of a downsampled VB3 analysis comparable to data that has not
% been downsampled. 
% M.L.2010-12-03

if(~exist('Tmax','var')|| isempty(Tmax))
    Tmax=1+ds*(size(v,1)-1);
end

V=zeros(Tmax,size(v,2));

V(1,:)=v(1,:);
for k=2:size(v,1)
   r1=1+(k-2)*ds;
   for m=1:ds
       V(r1+m,:)=v(k,:);
   end
end
V=V(1:Tmax,:);

