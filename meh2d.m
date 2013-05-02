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
% classes - classification into mesohyperbolic types
% compr - numerical compressibility
% spectral.dets - determinants (if passed as input, this is just a copy)
% spectral.traces - traces (if passed as input, this is just a copy)

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
quants.Hyp = T^2 * (Dets - 4) .* Dets; % mesohyperbolic when positive

% non-normality: Frobenius norm of the commutator of Jacobian
quants.NonNml = cellfun( @(A)norm( ctranspose(A)*A - A*ctranspose(A), 'fro' ), J  );

% defect: smallest distance between roots of the minimal polynomial
quants.Defect = cellfun( @(p) mindist( roots(p) ), Mp ); % mindist.m distributed with package

% compressibility
quants.Compr = T * Dets + Traces;

%% computation of mesohyperbolicity
classes_ind.hyp = quants.Hyp > 0;

classes_ind.mh_flipping = (Dets > 4/(T^2)) & classes_ind.hyp;
classes_ind.mh_pure = (Dets < 0) & classes_ind.hyp;

% assign pseudocolors
classes = nan(size(Dets));
classes(classes_ind.mh_flipping) = 1;
classes(classes_ind.mh_pure) = -1;
classes(~classes_ind.hyp) = 0;

spectral.Dets = Dets;
spectral.Traces = Traces;
