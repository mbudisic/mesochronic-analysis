function retval = mcjacobian(h, Ji, order)
% MCJACOBIAN
%
% Compute the mesochronic Jacobian evolved over the span t.
%
% h - timestep (double)
% Ji - cell-array of instantaneous Jacobian matrices evaluated at times
%      uniformly spaced by h
% order - order of integration, currently 1 or 2

if isempty(find( order == [1,2], 1))
    error('Order needs to be either 1 or 2')
end

Jisize = size(Ji);
Nt = Jisize(3);

% storage for mesochronic jacobian steps; index 1 is current, index 2 is previous step, index 3 two steps ago etc.
F = zeros( [Jisize(1:2), order+1] ); 

% zeroth step
F(:,:,1) = Ji(:,:,1);

% step loop - this is not an index
for n = 1:Nt-1
    
    % move F a step further, dropping oldest record, and inserting empty to the
    % front
        F = stepme1(F);
    
    % degrade order if we are not far enough to use higher order method
    % this way, lower order methods automatically initialize higher ones
    effectiveorder = min( [n, order] );
    
    if effectiveorder == 1
        F(:,:,1) = F(:,:,2) + h*updateme(n-1, F(:,:,2), Ji, h);
    elseif effectiveorder == 2
        F(:,:,1) = F(:,:,2) + (3*h/2) * updateme(n-1,F(:,:,2), Ji, h) - (h/2) * updateme(n-2, F(:,:,3), Ji, h);
    else
        error('Effective order reported as %d', effectiveorder)
    end
    
end

% return the most recent F
retval = F(:,:,1);


%%%%%%%%%%%%%%%% Three types of steppers - stepme1 seems to be the fastest
function retval = stepme1(mf)
retval = cat(3,zeros(size(mf,1),size(mf,2)),mf(:,:,1:end-1));

% function mf = stepme2(mf)
% % first type of stepping
% 
% for k = length(mf):-1:2
%     mf{k} = mf{k-1};
% end
% mf{1} = [];
% 
% function retval = stepme3(mf)
% % second type of stepping
% 
% retval = cell(size(mf));
% retval(2:end) = mf(1:end-1);
%%%%%%%%%%%%%%%%%

function retval = updateme(k, F, Ji, h)
% this is basically derivative of the mesochronic jacobian
Ji_step = Ji(:,:,k+1);
if k == 0
    retval = Ji_step * F;
elseif k > 0
    retval = (Ji_step-F)/(k*h) + Ji_step * F;
else
    error('k has to be non-negative')
end


