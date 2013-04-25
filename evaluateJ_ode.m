function [mJ, sol] = evaluateJ_ode( order, ic, f, t0, T, h, dp )
% EVALUATEJ_ODE( order, ic, f, t0, T, h, dp )
%
% Evaluate mesochronic Jacobian using direct method.
%
% order - order of the method used
% ic - a single initial condition (column vector)
% f - flow field (vectorized function handle)
% t0 - initial time (scalar)
% T - vector of integration lengths (positive)
% h - resampling time
% dp - finite difference used for evaluation of the instantaneous Jacobian

validateattributes(ic, {'numeric'}, {'column'})
validateattributes(t0, {'numeric'}, {'scalar'})

assert( min(diff(T)) > h, 'return steps should differ by more than integration step')
assert( h <= min(T), 'time step should be smaller than integration length' )
assert( dp < 1, 'dp should be small')

Tmax = max(T);

t = t0 + ( 0:h:Tmax );

% simulate
opts = odeset('vectorized','on','RelTol',1e-4,'AbsTol',1e-6);
[tout, yout] = ode23t(f, t, ic, opts);

sol.t = tout;
sol.x = yout;

assert( norm( t - tout.', Inf) < 1e-12, 'Output times do not match input times');


% evaluate Jacobian using finite differences
Ji = jacobian_fd(f, tout, yout, dp);

% return steps are multiples of resampling interval
retstep = fix(T/h);

% solve the Jacobian equation
mJ = mcjacobian_mex(h, Ji, retstep, order);

