function [classes,compr] = meh2d( T, varargin )
% meh2d( T, tol , ... )
%
% inputs: 
%
% T - integration period
% df - matrix/vector of determinants of Jacobians
% tf - matrix/vector of traces of Jacobians
%
% OR
%
% T - integration period
% J - cell array of 2x2 Jacobian matrices (if a single matrix is to be
%     analyzed, pass it as {J})
%

validateattributes(T, {'numeric'},{'scalar','nonnegative','finite'} );

%% validate arguments
assert( numel(varargin)  <= 2, 'Too many input arguments (consult help)' )
assert( numel(varargin)  > 0, 'Not enough input arguments (consult help)' )

% inputs are Jacobian matrices
if numel(varargin) == 1
    J = varargin{1};
    validateattributes( J, {'cell'},{})
    for k = 1:numel(J)
        validateattributes( J{k} , {'numeric'}, {'square', '2d', 'nrows', 2} ); % make sure inputs are 2d 2x2 matrices
    end
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
mh_pure = df > 4/(T^2);
not_mh = ~( mh_flipping | mh_pure );

classes = nan(size(df));
classes(mh_flipping) = -1;
classes(mh_pure) = 1;
classes(not_mh) = 0;


