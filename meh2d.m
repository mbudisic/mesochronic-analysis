function [classes,compr] = meh2d( T, tol , varargin )
% meh3d( T, tol , df, mf )
%
% inputs: 
%
% T - integration period
% tol - tolerance of zero (set to 0 for default value)
% df - matrix/vector of determinants of Jacobians
% tf - matrix/vector of traces of Jacobians
%
% OR
%
% T - integration period
% tol - tolerance of zero (set to 0 for default value)
% J - cell matrix/vector of 3x3 Jacobian matrices
%

validateattributes(T, {'numeric'},{'scalar','nonnegative','finite'} );
validateattributes(tol, {'numeric'},{'scalar','nonnegative','finite'} );

if tol == 0
tol = 3e-3;
fprintf(1, 'Tolerance to zero %.1e\n',tol);
end

%% validate arguments
assert( numel(varargin)  <= 2, 'Too many input arguments' )
assert( numel(varargin)  > 0, 'Not enough input arguments' )

% inputs are Jacobian matrices
if numel(varargin) == 1 
    
    J = varargin{1};
    for k = 1:numel(J)
        validateattributes( J{k} , {'numeric'}, {'square', '2d', 'nrows', 3} ); % make sure inputs are 2d 3x3 matrices
    end
    P = squeeze( charpoly_sequence(J) ); % extract characteristic polynomials
    
    % extract determinant and sum of minors
    if iscell(P) 
%        tf = cellfun( @(x)(-x(2)), P );
       tf = cellfun( @(x)(-x(2)), P );
       df = cellfun( @(x)(x(3)), P );       
    else
       tf = -P(:,2);
       df = P(:,3);
    end

% inputs are df and mf - validate them
else
    df = varargin{1};
    tf = varargin{2};
end

validateattributes( df, {'numeric'},{'real','finite'});
validateattributes( tf, {'numeric'},{'real','finite'});

assert( all( size(df) == size(tf) ), 'df and tf should have the same sizes');


%% computation of mesohyperbolicity    

% incompressibility deviation = "compressibility"

compr = tf + df * T;
mh_flipping = df < 0;
mh_pure = df > (4/T^2);
not_mh = (df > 0) & (df < (4/T^2));

classes = nan(size(df));
classes(mh_flipping) = 1;
classes(mh_pure) = -1;
classes(not_mh) = 0;


