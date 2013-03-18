function [retval_jacs, retval_steps] = mcjacobian(h, Ji, retstep, order)
% MCJACOBIAN
%
% Compute the mesochronic Jacobian evolved over the span t.
%
% h - timestep (double)
% Ji - cell-array of instantaneous Jacobian matrices evaluated at times
%      uniformly spaced by h
% retstep - spacing between returned steps
%            - if 0, only the last step is returned
%            - if > 0, every retstep is returned, with last step always
%                      returned
% order - order of integration, currently 1 or 2
%
% returns:
%
% retval_jacs - mesochronic jacobians taken at requested steps
% retval_steps - requested steps, corresponding to jacobians in retval_jacs

if isempty(find( order == [1,2], 1))
    error('Order needs to be either 1 or 2')
end

Jisize = size(Ji);
Nt = Jisize(3);

% storage for mesochronic jacobian steps; index 1 is current, index 2 is previous step, index 3 two steps ago etc.
F = zeros( [Jisize(1:2), order+1] ); 

% zeroth step
F(:,:,1) = Ji(:,:,1);


assert(retstep >= 0, 'Number of return value downsampling (retstep) needs to be positive')
if retstep > 0
    retlength = floor( Nt / retstep);
elseif retstep == 0
    retlength = 1;
else
    error('Retstep not properly formatted')
end

% storage for return values
retval_steps = zeros(1, retlength);
retval_jacs = zeros( Jisize(1), Jisize(2), retlength);

maxStep = Nt-1;
savecounter = 1;
% step loop - this is not an index
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
    if retstep > 0 && mod( maxStep - n, retstep ) == 0
        retval_steps(savecounter) = n;
        retval_jacs(:,:,savecounter) = F(:,:,1);
        savecounter = savecounter+1;
    end
end

% if only the last step was requested for saving
if retstep == 0
    retval_steps = n;
    retval_jacs = F(:,:,1);
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


