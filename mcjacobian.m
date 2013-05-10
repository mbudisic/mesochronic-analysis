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
% order - order of integration, currently 1,2,3
%       - if negative, highest allowed order is selected
%
% returns:
%
% retval_jacs - mesochronic jacobians taken at requested steps
% retval_steps - requested steps, corresponding to jacobians in retval_jacs

if h <= 0
    error('Timestep has to be positive')
end

%% choose order of the method
[abrhs, ablhs] = adamsbashforth;

if order < 0
    order = length(ablhs);
elseif order > length(ablhs)
    error('Maximum available Adams-Bashforth order exceeded. Set to -1 for automatically setting to maximal order implemented')
end

%% storage for intermediate results
Jisize = size(Ji);
Nt = Jisize(3);

% storage for mesochronic jacobian steps; index 1 is current, index 2 is previous step, index 3 two steps ago etc.
F = zeros( [Jisize(1:2), order+1] );

% zeroth step
F(:,:,1) = Ji(:,:,1);

%% determine which steps to return back to caller
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

%% main update loop
maxStep = Nt-1;
% step loop - zero indexed, 0 is the initial time (not computed)
%             maxStep is the final time
for n = 1:maxStep
    
    % move F a step further, dropping oldest record, and inserting empty to the
    % front
    F = shiftstep(F);
    
    % Degrade order if we have not computed enough steps to use a higher order method.
    % This way, lower order methods automatically initialize higher ones.
    eord = min( [n, order] ); % effective order
    
    % form multiplication coefficient
    coeff = h * abrhs(eord, 1:eord)/ablhs(eord);
    
    % approximate right hand side of the evolution
    deltaJ = coeff(1) * dJ_dt(n-1, F(:,:,2), Ji, h);
    for bstep = 2:eord
        deltaJ = deltaJ + coeff(bstep) * dJ_dt(n-bstep, F(:,:,1+bstep), Ji, h);
    end
    
    % update jacobian
    F(:,:,1) = F(:,:,2) + deltaJ;
    
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


%% auxiliaries

function retval = shiftstep(mf)
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

function [rhs, lhs] = adamsbashforth
% Adams-Bashforth coefficients (Petzold, Ascher, Linear Multistep,
%                               table 5.1)

maxorder = 6;

rhs = nan(maxorder,maxorder);
lhs = zeros(1,maxorder);

order=1;
rhs(order, 1:order) = 1;
lhs(order) = 1;

order=2;
rhs(order, 1:order) = [3, -1];
lhs(order) = 2;

order=3;
rhs(order, 1:order) = [23, -16, 5];
lhs(order) = 12;

order=4;
rhs(order, 1:order) = [55, -59, 37, -9];
lhs(order) = 24;

order=5;
rhs(order, 1:order) = [1901, -2774, 2616, -1274, 251];
lhs(order) = 720;

order=6;
rhs(order, 1:order) = [4277, -7923, 9982, -7298, 2877, -475];
lhs(order) = 1440;
