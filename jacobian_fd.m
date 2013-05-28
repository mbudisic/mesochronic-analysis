function Ji = jacobian_fd(f, t, y, dp)
% JACOBIAN_FD(f, t, y, dp)
%
% Evaluate the instantaneous jacobian by central finite difference
% along trajectory (t_i, y_i)
%
% f - @(t,x) vector field
% t - Npoints x 1 column vector of times
% y - Npoints x Ndim matrix of states
% dp - finite variation
%
% Returns:
% Ji - Ndim x Ndim x Npoints matrix
%
%
% optimal delta is 
% dp = ( 3 * eps * norm(f) / 2 norm( d^3 f )  )^(1/3)
% where norm(f) is infty-bound on values of f
%       norm(d^3 f) is infty-bound on third derivative of f
%       eps - machine precision (eps command in matlab)

validateattributes(t, {'double'}, {'real','column'});

Npoints = size(t,1);
assert( size(y,1) == Npoints, 'Number of states has to match length of time vector');
Ndim = size(y,2);

% output matrix
Ji = zeros(Ndim,Ndim,Npoints);

% variations
p_var = kron( eye(Ndim), [.5,-.5]) * dp; % Ndim x 2Ndim

% evaluate Jacobian point by point - this can probably be vectorized
for k = 1:Npoints
    
    t_p = t(k);
    point = y(k,:).';
    
    f_var = zeros(size(p_var));
    
    for n = 1:2*Ndim
        f_var(:,n) = f( t_p, point + p_var(:,n) );
    end
    
    Ji(:, :, k) = ( f_var(:,1:2:end) - f_var(:,2:2:end) )/dp;
end



