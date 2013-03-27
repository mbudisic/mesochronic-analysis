function mJ = evaluateJ_ode( order, ic, f, Jf, T, h )
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

Tmax = max(T);

t = 0:h:Tmax;
tc = num2cell(t);
% simulate
opts = odeset('vectorized','on');
[tout, yout] = ode23t(f, t, ic, opts);

assert( norm( t - tout.', Inf) < 1e-12, 'Output times do not match input times');

%fprintf(1, 'Min dt:%e, Max dt:%e \n', min(diff(S.x)), max(diff(S.x)));

% uniform resampling
y = num2cell(yout.',1);


% jacobians
Ji = cellfun(Jf, tc, y , 'UniformOutput', false);
Ji = cat(3, Ji{:});

% return steps are multiples of resampling interval
retstep = fix(T/h);

%    [mJ, ti] = mcjacobian_mex(dt, cat(3,Ji{:}), 1000, 2);
mJ = mcjacobian_mex(h, Ji, retstep, order);
