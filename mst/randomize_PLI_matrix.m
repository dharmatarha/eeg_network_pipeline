function [M]=randomize_PLI_matrix(M)

[a,b]=size(M);
if a~=b
    error('Only square matrices please!');
end

%get upper triangle values
values=M(triu(M)~=0);
%randomize values
values=values(randperm(length(values)));
%fill upper triangle with randomized values
M(triu(M)~=0)=values;
%mirror upper triangle over diagonal
M=triu(M)+triu(M,1)';
