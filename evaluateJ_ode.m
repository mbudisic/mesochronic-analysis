function [mJ, retstep, sol, order] = evaluateJ_ode( order, ic, f, t0, T, direction, h, dp )
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
[Ji, fi] = jacobian_fd(fsim, tout, yout, dp);

sol.Ji = Ji;
sol.fi = fi;

% return steps are multiples of resampling interval
retstep = ceil(T/h);

% solve the Jacobian equation
try
    [mJ, ~, order] = mcjacobian_mex(h, Ji, retstep, order);    
catch
    
    if verLessThan('matlab', '8')
        warning('MATLAB:MesochronicAnalysis:Compile', ...
            sprintf(['Cannot run compiled version of code on releases older than 2012b.\n',...
            'Running uncompiled code.\n',...
            'Turn this warning off by issuing: warning(''off'', ''MATLAB:MesochronicAnalysis:Compile'')']));
    else
        warning('MATLAB:MesochronicAnalysis:Compile', ...
            sprintf(['Compiled version does not exist.\nPlease run ''deploytool -build mcproject.prj'' prior to running for faster computations.\n',...
            'Running uncompiled code.\n',...
            'Turn this warning off by issuing: warning(''off'', ''MATLAB:MesochronicAnalysis:Compile'')']));
    end
    
    
    [mJ, ~, order] = mcjacobian(h, Ji, retstep, order);
end
