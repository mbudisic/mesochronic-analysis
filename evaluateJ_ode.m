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

t = 0:h:T;
tc = num2cell(t);
% simulate
opts = odeset('vectorized','on');
[~, yout] = ode23t(f, t, ic, opts);

%fprintf(1, 'Min dt:%e, Max dt:%e \n', min(diff(S.x)), max(diff(S.x)));

% uniform resampling
y = num2cell(yout.',1);


% jacobians
Ji = cellfun(Jf, tc, y , 'UniformOutput', false);
Ji = cat(3, Ji{:});

%    [mJ, ti] = mcjacobian_mex(dt, cat(3,Ji{:}), 1000, 2);
mJ = mcjacobian_mex(h, Ji, 0, order);
