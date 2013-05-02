function d = mindist(v)
% MINDIST
%
% Compute minimal distance between elements of v.

if ~iscolumn(v)
    v = v.';
end

if max(abs(v)) == 0.
    d = 0.;
else
    % repeat vector to the matrix
    V = repmat( v, [1, numel(v)]);
    
    % pairwise distance between elements
    D = abs(V - V.');
    
    % minimal off-diagonal distance (diagonal is obviously zero)
    d = min( D( eye(numel(v)) ~= 1 ) );
    
end