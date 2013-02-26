function ST=transposeStruct(S)
% ST=VB7_transpose(S)
% transpose S if S is a numeric type, or all numeric fields and subfields
% of S, if S is a structure. 

%M.:L. 2011-02-25

% basic case:
ST=S;
if(isnumeric(S))
    ST=S';
elseif(isstruct(S))
    fn=fieldnames(S);
    for k=1:length(fn)
       ST.(fn{k})=transposeStruct(S.(fn{k})); 
    end
else % in case S is neither a struct or a double, 
    ST=S;
end

