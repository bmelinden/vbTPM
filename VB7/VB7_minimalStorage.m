function W=VB7_minimalStorage(W) 
% remove bulky fields from the VB7 structure W to make storage to disk more
% efficient
if(isfield(W,'est2'))
    W=rmfield(W,'est2');
end
if(isfield(W,'minimalStorage'))
    W=rmfield(W,'minimalStorage');
end
