function d = mindist(v)
% MINDIST
%
% Compute minimal distance between elements of v.

v = v(:);

if all( v )  == 0.
    d = 0.;
else
    % repeat vector to the matrix
    V = repmat( v, [1, numel(v)]);
    
    % pairwise distance between elements - 
    % need only nonzero elements from upper triangular part,
    % as matrix is symmetric, and diagonal always contains zeros
    D = nonzeros( triu( abs(V - V.'), 1 ) );

    % minimal off-diagonal distance 
    if isempty(D)
        d = 0;
    else
        d = min(D);
    end
        
end