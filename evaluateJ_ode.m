function [mJ, sol] = evaluateJ_ode( order, ic, f, t0, T, direction, h, dp )
% EVALUATEJ_ODE( order, ic, f, t0, T, h, dp )
%
% Evaluate mesochronic Jacobian using ODE evolution.
%
% order - order of the method used
% ic - a single initial condition (column vector)
% f - flow field (vectorized function handle)
% t0 - initial time (scalar)
% T - vector of integration lengths (positive)
% direction (nonzero) - positive for forward in time, negative for backward
% h - resampling time
% dp - finite difference used for evaluation of the instantaneous Jacobian

validateattributes(ic, {'numeric'}, {'column'})
validateattributes(t0, {'numeric'}, {'scalar','real'})
validateattributes(direction, {'numeric'}, {'scalar', 'real', 'nonzero'});
validateattributes(T, {'numeric'}, {'positive'});

T = sort(T);
assert( min(diff([t0;T(:)])) > h, 'return steps should differ by more than integration step')
assert( h <= min(T), 'time step should be smaller than integration length' )
assert( dp < 1, 'dp should be small')

Tmax = max(T);

t = t0 + ( 0:h:Tmax );

% simulate
opts = odeset('vectorized','on','RelTol',1e-4,'AbsTol',1e-6);

if direction < 0
    fsim = @(t,x)(-f(-t,x));
else
    fsim = f;
end

[tout, yout] = ode23t(fsim, t, ic, opts);

sol.t = tout;
sol.x = yout;

assert( norm( t - tout.', Inf) < 1e-12, 'Output times do not match input times');


% evaluate Jacobian using finite differences
Ji = jacobian_fd(fsim, tout, yout, dp);

% return steps are multiples of resampling interval
retstep = fix(T/h);

% solve the Jacobian equation
if exist('mcjacobian_mex') ~= 3 && exist('deploytool') == 2
    mypath = mfilename('fullpath');
    mypath = mypath(1:find('/' == mypath, 1, 'last'));
    mcpath = [mypath 'mcjacobian.prj'];
    worker = getCurrentTask();
    if ~isempty(worker)
        error('deploying cannot be done automatically in multithreaded mode. manually compile mcjacobian.prj or re-run command in a single thread mode')
    else
    warning('mcjacobian can be compiled into MEX file. Attempting to build mcjacobian.prj');
    eval(['deploytool -build ' mcpath]);
    error('deploying completed. re-run');
    end
end

try
    mJ = mcjacobian_mex(h, Ji, retstep, order);
catch
    warning(['Runnin slow, non-compiled mcjacobian. mcjacobian.prj did not build on its own. It can be compiled by going to ' mypath])
    mJ = mcjacobian(h, Ji, retstep, order);
end
