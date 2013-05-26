function d = mindist(v)
% MINDIST
%
% Compute minimal distance between elements of v.

Nel = numel(v);
v = v(:);

if max(v) == 0.
%if all( v == 0 )
    d = 0.;
else
    % repeat vector to the matrix
    V = repmat( v, [1, Nel]);
    
    % off-diagonal distances between elements (diagonal is obviously zero)
    D = triu( abs(V - V.'), 1);
    
    d = min( D(:) );
end