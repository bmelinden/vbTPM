function x=TPM_RP_concatenate(lacX,lacY)
% x=TPM_RP_concatenate(lacX,lacY)
%
% function to transform raw data in the form 
%  lacX: {[M1 x120x B1 double]  ... [M4 x120x B4 double] }
%  lacY: {[M1 x120x B1 double]  ... [M4 x120x B4 double] }
% to 
%  x{i}{j} = [Mi*120 x 2 double], where i is the set index, and j is the bead
%
% This is the raw data format of xy-coordinate data saved by the RP group
% analysis software (as of now)
%
% M.L. 2010-06-04

x=cell(size(lacX));
for i=1:length(lacX) % loop over bead sets
    if(~isempty(lacX{i}))
        x{i}=cell(1,size(lacX{i},3));
        for j=1:length(x{i}); % loop over beads in set k
            x{i}{j}=zeros(120*size(lacX{i},1),2);
            for k=1:size(lacX{i},1)
                x{i}{j}(120*(k-1)+(1:120),1:2)=[lacX{i}(k,:,j)' lacY{i}(k,:,j)'];
            end
        end
    end
end
