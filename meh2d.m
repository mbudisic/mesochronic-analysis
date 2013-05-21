function [classes, quants, spectral] = meh2d( Ts, Js )
% meh2d( T, tol , ... )
%
% inputs:
%
% Ts - vector of integration periods (length K)
% Js - mesochronic Jacobian, 2 x 2 x K
%
%
% all the return fields are vectors of size of T
%
% -- classes
%    Mesohyperbolicity class
%    -1 - hyperbolicity, orientation preserving
%    0  - ellipticity
%    1  - hyperbolicity, orientation reversing
%
% -- quants.Compr
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



validateattributes(Ts, {'numeric'},{'nonnegative','finite'} );
validateattributes(Js, {'numeric'}, {});

D = 2;
assert( size(Js, 1) == D && size(Js, 2) == D, 'Jacobians are not 2x2xK matrices')
assert( numel(Ts) == size(Js,3), 'Number of Jacobians and integration periods has to be the same');

K = numel(Ts);

% storage
classes = nan(size(Ts));
quants.Compr = zeros(size(Ts));
quants.FTLE = zeros(size(Ts));
quants.Hyp = zeros(size(Ts));
quants.NonNml = zeros(size(Ts));
quants.NonDefect = zeros(size(Ts));

spectral.Traces = zeros(size(Ts));
spectral.Dets = zeros(size(Ts));

for k = 1:K
    
    % select appropriate averaging interval and Jacobian matrix
    T = Ts(k);
    J = Js(:,:,k);
    
    %% Compute quantifiers
    
    % Finite-Time Lyapunov Exponents
    quants.FTLE(k) = ftle(T,J);
    
    % non-normality: Frobenius norm of the commutator of Jacobian
    quants.NonNml(k) = norm( ctranspose(J)*J - J*ctranspose(J), 'fro' );
    
    % defect: smallest distance between roots of the minimal polynomial
    quants.NonDefect(k) = mindist( roots( minpoly(J)  ) ); % mindist.m distributed with package
    
    % characteristic polynomial
    Cp = poly(J);
    Det = Cp(3);
    Trace = -Cp(2);
    spectral.Traces(k) = Trace;
    spectral.Dets(k) = Det;
    
    % compressibility
    quants.Compr(k) = T * Det + Trace;
    
    % mesohyperbolicity
    quants.Hyp(k) = (T^2 * Det - 4) .* Det;
    
    %% computation of mesohyperbolicity
    
    if quants.Hyp(k) > 0 % mesohyperbolicity
        
        if Det >= 4/(T^2) % orientation reversing
            classes(k) = 1;
        elseif Det <= 0     % orientation preserving
            classes(k) = -1;
        else
            error('Detected mesohyperbolicity, but cannot detect the class');
        end
    else
        classes(k) = 0; % mesoellipticity
    end
end