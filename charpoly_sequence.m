function P = charpoly_sequence( M )
% CHARPOLY_SEQUENCE
% Compute characteristic polynomials of a sequence of matrices.
%
% M is either a cell array of square matrices
% or a (N,N,M) array, where each (:,:,k) is a matrix whose characteristic
% polynomial will be evaluated.

validateattributes(M, {'numeric','cell'}, {})

% if input matrix was not a cell, convert it to cell
if ~iscell(M)
    validateattributes(M, {'numeric'}, {})
    if  size(M,3) > 1
        M = num2cell(M, [1,2]);
    else
        M = {M};
    end
end


% assert all matrices are numeric and square
for k = 1:numel(M)
    validateattributes(M{k}, {'numeric'}, {})
end

% compute characteristic polynomials
P = squeeze(cellfun( @poly, M, 'uniformoutput',false ));

if isvector(P)
try
    P = cell2mat(P);
catch
end
end