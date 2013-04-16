function [mJ, sol] = evaluateJ_ode( order, ic, f, T, h, dp )
% EVALUATEJ_ODE
%
% Evaluate mesochronic Jacobian using direct method.
%
% order - order of the method used
% ic - vector of initial condition
% f, Jf - flow field and its jacobian (function handles)
% T - vector of integration lengths
% h - resampling time

validateattributes(ic, {'numeric'}, {'column'})

assert( h < min(T), 'time step should be smaller than integration length' )
assert( dp < 1, 'dp should be small')

Tmax = max(T);

t = 0:h:Tmax;
tc = num2cell(t);
% simulate
opts = odeset('vectorized','on','RelTol',1e-4,'AbsTol',1e-6);
[tout, yout] = ode23t(f, t, ic, opts);

sol.t = tout;
sol.x = yout;

assert( norm( t - tout.', Inf) < 1e-12, 'Output times do not match input times');

%fprintf(1, 'Min dt:%e, Max dt:%e \n', min(diff(S.x)), max(diff(S.x)));

% convert to cell
y = num2cell(yout.',1);
% 
% 
% % jacobians
% Ji = cellfun(Jf, tc, y , 'UniformOutput', false);
% Ji = cat(3, Ji{:});

% evaluate using finite differences
Ji = jacobian_fd(f, tout, yout, dp);

% return steps are multiples of resampling interval
retstep = fix(T/h);

%    [mJ, ti] = mcjacobian_mex(dt, cat(3,Ji{:}), 1000, 2);
mJ = mcjacobian_mex(h, Ji, retstep, order);


