function [classes,quant,classes_ind] = meh2d( T, tol , varargin )
% meh3d( T, tol , df, mf )
%
% inputs: 
%
% T - integration period
% tol - tolerance of zero (set to 0 for default value)
% df - matrix/vector of determinants of Jacobians
% mf - matrix/vector of sum of minors of Jacobians
%
% OR
%
% T - integration period
% tol - tolerance of zero (set to 0 for default value)
% J - cell matrix/vector of 3x3 Jacobian matrices
%
%
%  classes - pseudocolor matrix of the same size as df
%            -2 -> EFocus
%            -1 -> ENode
%             0 -> Non-hyp
%             1 -> FNode
%             2 -> FFocus
%  quant   - structure of quantitiers used to determine mesohyperbolicity
%          .Cplus
%          .Cminus
%          .Delta
%          .Sigma
%          .T - averaging time
%
%  classes_ind - logical matrices of the same size as df
%              - indicator-functions for each of the classes 

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
    P = charpoly_sequence(J); % extract characteristic polynomials
    
    % extract determinant and sum of minors
    if iscell(P) 
%        tf = cellfun( @(x)(-x(2)), P );
       mf = cellfun( @(x)(x(3)), P );
       df = cellfun( @(x)(-x(4)), P );       
    else
       mf = P(:,3);
       df = -P(:,4);
    end

% inputs are df and mf - validate them
else
    df = varargin{1};
    mf = varargin{2};
end

validateattributes( df, {'numeric'},{'real','finite'});
validateattributes( mf, {'numeric'},{'real','finite'});

assert( all( size(df) == size(mf) ), 'df and mf should have the same sizes');


%% computation of mesohyperbolicity    

% Delta determines the real/complex character of eigenvalues
% Sigma determines stability of the two eigenvalues with matching stability
% (see writeup)

quant.Cplus = -T^3*df;
quant.Cminus = 3*T^3*df+2*T^2*mf-8;

% Evaluate Delta using Horner's scheme to reduce numerical errors.
% Coefficients listed in descending order:
Delta_coeff = { -4 * df.^4, ...
                -12 * df.^3 .* mf, ...
                -13 * df.^2 .* mf.^2, ...
                -6 * df .* mf .^3, ...
                -mf .* (mf.^3 - 18 * df.^2), ...
                18*df.*mf.^2, ...
                27 * df.^2 + 4 * mf.^3 };

Delta = zeros( size( df) );
for k = 1:length(Delta_coeff)
    Delta = (Delta * T) + Delta_coeff{k};
end
quant.Delta = Delta;
quant.Sigma = -df./( 3*T^3*df + 2*T^2*mf - 8. ) ;
quant.T = T;

% distance from nonhyperbolicity - scaled by T^3 to account for time
% factor
quant.nonhyp_distance = min( abs(quant.Cplus), abs(quant.Cminus) )/T^3;

% node/focus - dynamics of the pair of eigenvalues
% u/s stability of the pair of eigenvalues

% Cplus can avoid multiplying by T^3
classes_ind.hyp = quant.nonhyp_distance > tol;

classes_ind.node_f = (quant.Delta < 0) & (quant.Sigma < 0) & classes_ind.hyp;
classes_ind.node_e = (quant.Delta < 0) & (quant.Sigma > 0) & classes_ind.hyp;
classes_ind.focus_f = (quant.Delta > 0) & (quant.Sigma < 0) & classes_ind.hyp;
classes_ind.focus_e = (quant.Delta > 0) & (quant.Sigma > 0) & classes_ind.hyp;

% assign pseudocolors
classes = nan(size(df));
classes(classes_ind.focus_e) = -2;
classes(classes_ind.node_e) = -1;
classes(classes_ind.node_f) = 1;
classes(classes_ind.focus_f) = 2;
classes(~classes_ind.hyp) = 0;

fprintf(1,'MH F-Node # %d\n',sum(real(classes_ind.node_f(:))))
fprintf(1,'MH E-Node # %d\n',sum(real(classes_ind.node_e(:))))
fprintf(1,'MH F-Focus # %d\n',sum(real(classes_ind.focus_f(:))))
fprintf(1,'MH E-Focus # %d\n',sum(real(classes_ind.focus_e(:))))
fprintf(1,'MH Nonhyp # %d\n',sum(real(~classes_ind.hyp(:))))

summary  = @(x)prctile(x,[2.5 25 50 75 97.5]);

disp('Percentile summary: 2.5, 25, 50, 75, 97.5 %')
disp(['|Delta|: ' mat2str(summary(abs(quant.Delta(:))), 5 )])
disp(['|Sigma|: ' mat2str(summary(abs(quant.Sigma(:))), 5 )])
