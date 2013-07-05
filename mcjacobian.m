function [mcJacobians, mcSteps, orderUsed] = mcjacobian(h, Ji, mcStepsRequested, orderRequested)
% MCJACOBIAN
%
% [mcJacobians, mcSteps, orderUsed] = mcjacobian(h, Ji, mcStepsRequested, orderRequested)
%
% Compute the mesochronic Jacobian using Adams-Bashforth method.
% For reference, see chapter on Linear Multistep methods in
%
% U.M. Ascher and L.R. Petzold, 
% "Computer Methods for Ordinary Differential Equations and Differential-Algebraic Equations"
% (Society for Industrial and Applied Mathematics (SIAM), 
% Philadelphia, PA, 1998), pp. xviii?314.
%
% Adams-Bashforth (A-B) is a fixed timestep method, therefore duration of 
% its solution is given as T = Nt * h
% where Nt is the number of steps, and h is the duration of each timestep.
%
%
% This function applies A-B to evolution ODE of the mesochronic Jacobian.
% The inputs from which the increments are computed depend only on
% the evaluation of instantaneous Jacobians (Ji) along a trajectory.
% Typically, only the final-time mesochronic Jacobian is needed,
% or the values at a relatively-sparse set of intermediate timesteps.
% Therefore, to save memory, the function asks the user to specify
% steps of interest (mcStepsRequested) at which 
% evaluations of the mesochronic Jacobian are required.
%
% It is important to note that this function does not distinguish between
% time-varying and time-invariant systems, nor between forward-time and
% backward-time evolutions. The Ji are assumed to be evaluated along the
% equally spaced (value: h) intervals between 
% initial and terminal time instances. Whether this means
% t0, t0+h, t0+2h, ..., t0 + T
% or
% t0, t0-h, t0-2h, ..., t0 - T
% is irrelevant from the standpoint of this method.
%
% A-B is not a single method, but rather a family of methods, parametrized
% by the order of the method (orderRequested). Currently, orders 
% implemented range from 1 to 6.
%
% *** inputs: ***
% 
% h - length of a timestep (positive)
% Ji - array of instantaneous Jacobian matrices evaluated at times
%      uniformly spaced by h (each (:,:,k) is a Jacobian )
%    D x D x Nt matrix, where 
%    D is the dimension of the state space, 
%    Nt is the length of time evolution (maximum number of timesteps)
%
% mcStepsRequested - 1 x N vector of ordinal numbers of timesteps
%                    for which mesochronic Jacobian will be returned
%                  each element is a value between
%                  0 and Nt-1, where Nt is determined from length of Ji
%
%                  0 -- indicates the inital step
%                  Nt -- indicates the final step
%                  Additionally, user can pass any integer as a value
%                  for elements of mcStepsRequested.
%                  Values outside [0, Nt-1] range will be taken as
%                  modulo-Nt, to allow, 
%                  e.g., last step to be requested as -1,
%                        next-to-last step to be requested as -2, etc.
%                  the actual values of steps for which mc. Jacobians are
%                  returned is passed as an output
%                  
% orderRequested  - preferred order of the A-B method. 
%                   Allowed values:
%                   < 0 -- maximum order advisable to use
%                   1 to 6 -- currently implemented orders
%                   Not allowed values:
%                   0  -- reserved for finite-difference (not A-B) method
%                   > 6 -- not yet implemented
%
% *** outputs: ***
%
% mcJacobians - mesochronic Jacobian
%               D x D x N - where 
%               D is the state dimension,
%               N is the number of iteration steps returned
% mcSteps     - 1 x N -- orders of iteration steps returned
% orderUsed   - the actual order A-B method used; may or may not be equal 
%               to orderRequested (input) 
%

if h <= 0
    error('Timestep has to be positive')
end

%% choose order of the method

% return coefficients of all available orders
[abrhs, ablhs, maxOrder] = adamsbashforth;

assert( orderRequested ~= 0,... 
    ['Zeroth order of Adams-Bashforth requested. ',...
    'By internal convention, we use 0 to denote finite difference',...
    'method (not Adams-Bashforth). Check calling function for errors'])

% determine order that is used based on the request
%
if orderRequested < 0
    % < 0 -> orderUsed to maximum available
    
    % current implementation has roundoff (or other numerical) errors
    % in orders higher than 3, therefore we artificially set 3
    % as the maximum order unless user explicitly requests 4,5,6
    %
    % after debugging, change the following line to read
    % orderUsed = maxOrder;
    orderUsed = 3;
elseif orderRequested > maxOrder
    % > maxOrder - report implementation limitation
    
    error('Maximum available Adams-Bashforth order exceeded. Set to -1 for automatically setting to maximal order implemented.')
else
    % everything is fine, use requested order
    orderUsed = orderRequested;
end
    

%% storage for intermediate results
sizeOfJacobians = size(Ji); % NxN
Nt = sizeOfJacobians(3);    % N of time steps

% storage for mesochronic jacobian steps; index 1 is current, index 2 is previous step, index 3 two steps ago etc.
J_steps = zeros( [sizeOfJacobians(1:2), orderUsed+1] );

% zeroth step
J_steps(:,:,1) = Ji(:,:,1);

%% determine which steps to return back to caller
if isempty(mcStepsRequested)
    mcStepsRequested = -1;
end

% take the mod of retstep to rewrap
assert(isvector(mcStepsRequested), 'mcStepsRequested has to be a vector');
mcStepsRequested = mod(mcStepsRequested, Nt);
mcStepsRequested = unique(mcStepsRequested)+1; % +1 for zero-indexed 

numberOfStepsRequested = length(mcStepsRequested);

% storage for return values
mcSteps = zeros(1, numberOfStepsRequested);
mcJacobians = zeros( sizeOfJacobians(1), ...
                     sizeOfJacobians(2), ...
                     numberOfStepsRequested);

%% main update loop
maxStep = Nt-1;
% step loop - zero indexed, 0 is the initial time (not computed)
%             maxStep is the final time

dJ_storage = zeros( sizeOfJacobians(1), ...
                    sizeOfJacobians(2), ...
                    orderUsed);

for n = 1:maxStep
    
    % move F a step further, dropping oldest record, and inserting empty to the
    % front
    J_steps = shiftstep(J_steps);
    
    % Degrade order if we have not computed enough steps to use a higher order method.
    % This way, lower order methods automatically initialize higher ones.
    effectiveOrder = min( [n, orderUsed] ); % effective order
    
    % select coefficients based on the effective order
    myrhs = abrhs(effectiveOrder, 1:effectiveOrder);
    mylhs = ablhs(effectiveOrder);
    
    
    coeff = myrhs / mylhs;
    
    % index 2 is the base step (first in history) based on which we determine value of next step        
    % approximate right hand side of the evolution
    for bstep = 1:effectiveOrder
         dJ_storage(:,:,bstep) =  coeff(bstep) * dJ_dt(n-bstep, J_steps(:,:,1+bstep), Ji, h);
    end
            
    % update jacobian
    S_matlab = sum(h*dJ_storage(:,:, 1:effectiveOrder),3); % plain summation
    
%   kahan
%   summation - for debugging purposes only    
%   S_kahan = matrixKahanSum(h*dJ_storage(:,:, 1:effectiveOrder)); 

    J_steps(:,:,1) = J_steps(:,:,2) + S_matlab;
    
    % detect a saving step
    [~,savestep] = ismember(n, mcStepsRequested);
    
    if savestep > 0
        mcSteps(savestep) = n;
        mcJacobians(:,:,savestep) = J_steps(:,:,1);
    end
end


%% auxiliaries

function retval = shiftstep(mf)
% drop last layer of mf, and append a zero matrix to the first layer
% [a,b,...,m, n] ---> [0, a, b, ..., m]

retval = cat(3,...
    zeros(size(mf,1),size(mf,2)),...
    mf(:,:,1:end-1) );

function dJdt = dJ_dt(k, J, Jis, h)
% Function approximates the numerical derivative of the mesochronic
% Jacobian.
%
% input:
% k - current step of the iteration
% J - mesochronic Jacobian evaluated at this step
% Jis - entire evolution of instantaneous Jacobians
% h - magnitude of the timestep
%
% returns:
% dJdt - approximate derivative of the mesochronic Jacobian

% instantaneous Jacobian at current step
Ji = Jis(:,:,k+1);

if k == 0
    dJdt = Ji * J;
else
    dJdt = (Ji-J)/(k*h) + Ji * J;
end

function S = matrixKahanSum( M )
% S = matrixKahanSum( M )
%
% Kahan summation for matrices, over 3rd dimension
%
% This summation compensates for finite-precision roundoff effects (but not
% for overflow).
%

S = zeros( size(M,1), size(M,2) ); % result
c = zeros( size(M,1), size(M,2) ); % compensation

N = size(M,3);

for i = 1:N
   y = M(:,:,i) - c;
   T = S + y;
   c = (T - S) - y;
   S = T;
end


function [rhs, lhs, maxOrder] = adamsbashforth
% Adams-Bashforth coefficients (Petzold, Ascher, Linear Multistep,
%                               table 5.1)
%
% This is a "convenience" function, just making it easier
% to locate where coefficients can be found in case one wants to
% code up higher order methods.
%
% Currently, returns Adams-Bashforth coefficients up to maxorder=6.
%
% rhs is a maxorder x maxorder matrix of right-hand-side coefficients
% lhs is a 1 x maxorder vector of left-hand-side coefficients
%
% The reason why the function returns all coefficients, and not only the
% desired order, is because of initialization:
% to initialize Nth order, all coefficients of orders of 1 to N-1 are
% needed, and are used to compute first N-1 steps in the iteration
% (N-th order needs N initial conditions).
%
% To add coefficients of additional order, change maxorder line
% and then append three lines, e.g.,
%
% order=7;
% rhs(order, 1:order) = ...
% lhs(order) = ...
%
% This is clearly not the most efficient way to define the coefficients,
% but I found it to be more readable than just having a big matrix.

% maxorder is implementation dependent, and therefore not a parameter
maxOrder = 6; rhs = nan(maxOrder,maxOrder);

lhs = zeros(1,maxOrder);

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
