function A=rowNormalize(Q)
% normalize the rows af a transition matrix to enforce sum_j Q(i,j) = 1
A=Q;
for k=1:size(Q,1)
    A(k,:)=Q(k,:)/sum(Q(k,:));
end