function Ji = jacobian_fd(f, t, y, delta)
% JACOBIAN_FD(f, t, y, delta)
%
% Evaluate the instantaneous jacobian by finite difference.
%
% f - @(t,x) vector field
% t - Npoints x 1 column vector of times
% y - Npoints x Ndim matrix of states
% delta - finite variation

validateattributes(t, {'double'}, {'real','column'});

Npoints = size(t,1);
assert( size(y,1) == Npoints, 'Number of states has to match length of time vector');
Ndim = size(y,2);

Ji = zeros(Ndim,Ndim,Npoints);

p_var = kron( eye(Ndim), [1,-1]) * delta; % Ndim x 2Ndim

for k = 1:Npoints
    
    t_p = t(k);
    point = y(k,:).';
    
    f_var = zeros(size(p_var));
    
    for n = 1:2*Ndim
        f_var(:,n) = f( t_p, point + p_var(:,n) );
    end
    
    Ji(:, :, k) = ( f_var(:,1:2:end) - f_var(:,2:2:end) )/(2*delta);
end



