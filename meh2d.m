function [classes, compr, spectral] = meh2d( T, J )
% meh2d( T, tol , ... )
%
% inputs:
%
% T - integration period
% J - cell array of 2x2 Jacobian matrices (if a single matrix is to be
%     analyzed, pass it as {J})
%
%
% returns:
% classes - classification into mesohyperbolic types
% compr - numerical compressibility
% spectral.dets - determinants (if passed as input, this is just a copy)
% spectral.traces - traces (if passed as input, this is just a copy)

validateattributes(T, {'numeric'},{'scalar','nonnegative','finite'} );
validateattributes(J, {'cell'}, {});
cellfun( @(x)validateattributes(x, {'numeric'}, {}), J)

% inputs are Jacobian matrices
P = charpoly_sequence(J); % extract characteristic polynomials
% extract determinant and sum of minors
if iscell(P)
    %        tf = cellfun( @(x)(-x(2)), P );
    tf = cellfun( @(x)(-x(2)), P );
    df = cellfun( @(x)(x(3)), P );
else
    tf = -P(:,2);
    df = P(:,3);
end

% normalcy
nml = normalcy( J );

validateattributes( df, {'numeric'},{'real','finite'});
validateattributes( tf, {'numeric'},{'real','finite'});

assert( all( size(df) == size(tf) ), 'df and tf should have the same sizes');


%% computation of mesohyperbolicity

% incompressibility deviation = "compressibility"

compr = tf + df * T;
mh_flipping = df < 0;
mh_pure = df > 4/(T^2);
not_mh = ~( mh_flipping | mh_pure );

classes = nan(size(df));
classes(mh_flipping) = -1;
classes(mh_pure) = 1;
classes(not_mh) = 0;

spectral.dets = df;
spectral.traces = tf;
spectral.nml = nml;

function retval = normalcy( M )

retval = zeros(size(M));
for k = 1:numel(M)
    retval(k) = norm( M{k}*ctranspose(M{k}) -  ctranspose(M{k}) * M{k} );
end

