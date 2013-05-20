function [classes, quants, spectral] = meh2d( T, J )
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
% 
% -- classes
%    Mesohyperbolicity class
%    -1 - hyperbolicity, orientation preserving
%    0  - ellipticity
%    1  - hyperbolicity, orientation reversing
%
% -- quants.compr
% Compressibility of mesochronic Jacobian J (theoretically identical to 0):
%   T * det J + tr J
%   This is a coarse measure of numerical error.
%
% -- quants.NonNml
% Non-normality of mesochronic Jacobian J: Frobenius norm of the commutator
%   || J* . J - J . J*||
%   When non-normality is zero, J has complete orthogonal basis of eigenvectors 
%   (it is *unitarily* diagonalizable)p.
%
% -- quants.Hyp
% Hyperbolicity of mesochronic Jacobian J: 
%   (T^2 * det J - 4) .* det J
%   When positive, J has a pair of real eigenvalues (flow map is hyperbolic for the integration time).
%
% -- quants.NonDefect
% Non-defectiveness of mesochronic Jacobian J: 
%   smallest distance between roots of the minimal polynomial of J
%   When zero, J is defective, i.e., it has an 
%   incomplete basis of eigenvectors (it is not diagonalizable)
% 
% -- quants.FTLE
%    Finite Time Lyapunov Exponent 
%    (1/T) * log norm(J)
%    norm is 2-norm, T is always positive 
%
% -- spectral.Dets - determinants of mc. Jacobians
% -- spectral.Traces - traces of mc. Jacobians, 



validateattributes(T, {'numeric'},{'scalar','nonnegative','finite'} );
validateattributes(J, {'cell'}, {});
cellfun( @(x)validateattributes(x, {'numeric'}, {}), J)

% inputs are Jacobian matrices
[Cp, Mp] = matpolys(J); % extract characteristic and minimal polynomials

% 3d polynomial
% P = x^2 - x Trace + Det

if iscell(Cp)
    Traces = cellfun( @(x)(-x(2)), Cp );
    Dets = cellfun( @(x)(x(3)), Cp );
else
    Traces = -Cp(:,2);
    Dets = Cp(:,3);
end

if ~iscell(Mp)
    Mp = mat2cell(Mp, ones(1,size(Mp,1)), size(Mp,2));
end

%% Compute quantifiers

% Finite-Time Lyapunov Exponents
quants.FTLE = cellfun(  @(M)ftle(T,M), J );

quants.Hyp = (T^2 * Dets - 4) .* Dets; % mesohyperbolic when positive

% non-normality: Frobenius norm of the commutator of Jacobian
quants.NonNml = cellfun( @(A)norm( ctranspose(A)*A - A*ctranspose(A), 'fro' ), J  );

% defect: smallest distance between roots of the minimal polynomial
quants.NonDefect = cellfun( @(p) mindist( roots(p) ), Mp ); % mindist.m distributed with package

% compressibility
quants.Compr = T * Dets + Traces;

%% computation of mesohyperbolicity

% figure out which indices of initial conditions correspond to each class
classes_ind.hyp = quants.Hyp >= 0;

classes_ind.mh_flipping = (Dets >= 4/(T^2)) & classes_ind.hyp;
classes_ind.mh_pure = (Dets <= 0) & classes_ind.hyp;

% assign pseudocolors
classes = nan(size(Dets));
classes(classes_ind.mh_flipping) = 1;
classes(classes_ind.mh_pure) = -1;
classes(~classes_ind.hyp) = 0;

spectral.Dets = Dets;
spectral.Traces = Traces;
