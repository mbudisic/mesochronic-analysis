function retval = meh_analysis(T, Jacobians, Ndim, tol)
% meh_analysis(T, Jacobians)
%
% Mesochronic Analysis of a sequence of mesochronic Jacobians.
%
% T - vector of lengths of integration periods (length K)
% Jacobians - cell array of Ndim x Ndim x K or 3x3xK
% Ndim - dimension of state space, 2 or 3

% set up output storage structures
Npoints = numel(Jacobians);
Dets = zeros([Npoints, length(T)]);
Traces = zeros(size(Dets));
Meh = zeros(size(Dets));
Compr = zeros(size(Dets));
NonNml = zeros(size(Dets));
NonDefect = zeros( size(Dets));
TrCof = zeros(size(Dets)); % Trace of Cofactor - used for 3d analysis
FTLE = zeros( size(Dets));
Hyp =  zeros( size(Dets));

if Ndim == 3
    disp('Using 3d mesohyperbolicity')
    assert(exist('meh3d','file') == 2, '3d analysis not yet implemented');
    meh = @(T, J)meh3d( T, J, tol);
else
    disp('Using 2d mesohyperbolicity')
    meh = @(T,J)meh2d(T,J);
end

% evaluate quantifiers for Jacobians
parfor m = 1:Npoints
    
    % select the appropriate jacobian
    myJacobians = Jacobians{m};
    [classes, quants, spectral] = meh( T, myJacobians );
    
    % spectral quantities
    Dets(m,:) = spectral.Dets;
    Traces(m,:) = spectral.Traces;
    if isfield(spectral,'TrCofs')
        TrCof(m,:) = spectral.TrCofs;
    end
    
    % quantifiers of mesochronic analysis criteria
    Compr(m,:) = quants.Compr;
    NonNml(m,:) = quants.NonNml;
    NonDefect(m,:) = quants.NonDefect;
    FTLE(m,:) = quants.FTLE;
    Hyp(m,:) = quants.Hyp;
    
    % mesochronic classes
    Meh(m,:) = classes;
    
end

% create and store output structures
retval.Dets = Dets;
retval.Traces = Traces;
if Ndim == 3
    retval.TrCof = TrCof;
end

retval.Meh = Meh;
retval.Compr = Compr;
retval.NonNml = NonNml;
retval.FTLE = FTLE;
retval.Hyp = Hyp;
retval.NonDefect = NonDefect;
