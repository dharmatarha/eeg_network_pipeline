function [MST]=get_MST(M)

matrix_inverse=1./M;

electrode=15;

%with the inverse PLI values - the MATLAB MST algorithm works with the
%smallest values, but we're searching for the largest values!!!, the inverse is suitable for tree construction
matrix_inverse(isinf(matrix_inverse))=0;
inverse_PLI_values = create_sparse_matrix (matrix_inverse);

%get the MST, needs bioinformatics toolbox
[ST_inverse,pred] = graphminspantree(inverse_PLI_values,'Method', 'Kruskal');


MST_matrix_inverse = full(ST_inverse);
MST=ST_inverse;
[row col] = find(MST_matrix_inverse~=0);
if length(col)==electrode-1
    for a=1:length(col)
        MST(row(a),col(a))=M(row(a),col(a));
    end
end

MST=MST+MST';

