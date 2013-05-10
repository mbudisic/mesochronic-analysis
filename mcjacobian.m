function [retval_jacs, retval_steps] = mcjacobian(h, Ji, retstep, order)
% MCJACOBIAN
%
% Compute the mesochronic Jacobian evolved over the span t.
%
% h - timestep (double)
% Ji - array of instantaneous Jacobian matrices evaluated at times
%      uniformly spaced by h (each (:,:, k) is a Jacobian )
% steps - vector indices of returned steps, 0-started, mod(, Nsteps)
%       - 0 - first step
%       - -1 - last step
%       - 
% order - order of integration, currently 1 or 2
%
% returns:
%
% retval_jacs - mesochronic jacobians taken at requested steps
% retval_steps - requested steps, corresponding to jacobians in retval_jacs

if h <= 0
    error('Timestep has to be positive')
end

if isempty(find( order == [1,2], 1))
    error('Order needs to be either 1 or 2')
end

Jisize = size(Ji);
Nt = Jisize(3);

% storage for mesochronic jacobian steps; index 1 is current, index 2 is previous step, index 3 two steps ago etc.
F = zeros( [Jisize(1:2), order+1] ); 

% zeroth step
F(:,:,1) = Ji(:,:,1);

if isempty(retstep)
    retstep = -1;
end

% take the mod of retstep to rewrap
assert(isvector(retstep), 'retstep has to be a vector');
retstep = mod(retstep, Nt);
retstep = unique(retstep);

retlength = length(retstep);

% storage for return values
retval_steps = zeros(1, retlength);
retval_jacs = zeros( Jisize(1), Jisize(2), retlength);

maxStep = Nt-1;
% step loop - zero indexed, 0 is the initial time (not computed)
%             maxStep is the final time
for n = 1:maxStep
    
    % move F a step further, dropping oldest record, and inserting empty to the
    % front
        F = stepme1(F);
    
    % Degrade order if we have not computed enough steps to use a higher order method.
    % This way, lower order methods automatically initialize higher ones.
    effectiveorder = min( [n, order] );
    
    if effectiveorder == 1
        F(:,:,1) = F(:,:,2) + h*dJ_dt(n-1, F(:,:,2), Ji, h);
    elseif effectiveorder == 2
        F(:,:,1) = F(:,:,2) + (3*h/2) * dJ_dt(n-1,F(:,:,2), Ji, h) - (h/2) * dJ_dt(n-2, F(:,:,3), Ji, h);
    else
        error('Effective order reported as %d', effectiveorder)
    end
    
    % detect a saving step
    savestep = find(n == retstep, 1, 'first');
    if ~isempty(savestep)
        savestep = min(savestep);
    else
        savestep = NaN;
    end
    
    if ~isnan(savestep)
        retval_steps(savestep) = n;
        retval_jacs(:,:,savestep) = F(:,:,1);
    end
end



function retval = stepme1(mf)
% drop last layer of mf, and append a zero matrix to the first layer
% [a,b,...,m, n] ---> [0, a, b, ..., m]

retval = cat(3,...
        zeros(size(mf,1),size(mf,2)),...
        mf(:,:,1:end-1) );

function retval = dJ_dt(k, F, Ji, h)
% this is numerical derivative of the mesochronic jacobian
Ji_step = Ji(:,:,k+1);
if k == 0
    retval = Ji_step * F;
elseif k > 0
    retval = (Ji_step-F)/(k*h) + Ji_step * F;
else
    error('k has to be non-negative')
end


