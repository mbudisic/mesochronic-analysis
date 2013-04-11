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
% 
% 
% % jacobians
% Ji = cellfun(Jf, tc, y , 'UniformOutput', false);
% Ji = cat(3, Ji{:});

Ji = jacobian_fd(f, tout, yout.', 1e-9);

% return steps are multiples of resampling interval
retstep = fix(T/h);

%    [mJ, ti] = mcjacobian_mex(dt, cat(3,Ji{:}), 1000, 2);
mJ = mcjacobian_mex(h, Ji, retstep, order);

function Ji = jacobian_fd(f, tout, yout, delta)
% jacobian_fd
%
% evaluate the instantaneous jacobian by finite difference
%

Np = size(yout,2);
assert(size(yout,1) == 2)
Ji = zeros(2,2,Np);

for k = 1:Np
    
    t = tout(k);
    point = yout(:,k);
    
    assert(iscolumn(point))
    
    dvar = [ 1, 0; -1, 0; 0, 1; 0, -1].';
    
    for n = 1:size(dvar,2)
        f_var(:,n) = f( t, point + dvar(:,n) );
    end
    
    Ji(:, :, k) = [f_var(:, 1) - f_var(:,2), f_var(:, 3) - f_var(:,4)]/(2*delta);
end





