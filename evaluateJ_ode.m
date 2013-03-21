function mJ = evaluateJ_ode( order, ic, f, Jf, T, h )
% EVALUATEJ_ODE
%
% Evaluate mesochronic Jacobian using direct method.
%
% order - order of the method used
% ic - vector of initial condition
% f, Jf - flow field and its jacobian (function handles)
% T - integration length
% h - resampling time

validateattributes(ic, {'numeric'}, {'column'})

tc = 0:h:T;

% simulate
opts = odeset('vectorized','on');
S = ode23t(f, [0, T], ic, opts);

%fprintf(1, 'Min dt:%e, Max dt:%e \n', min(diff(S.x)), max(diff(S.x)));

% uniform resampling
y = num2cell( deval(tc, S), 1);

% jacobians
Ji = cellfun(Jf, num2cell(tc), y , 'UniformOutput', false);

%    [mJ, ti] = mcjacobian_mex(dt, cat(3,Ji{:}), 1000, 2);
mJ = mcjacobian_mex(h, cat(3,Ji{:}), 0, order);
