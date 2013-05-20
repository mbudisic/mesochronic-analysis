function ftle_val = ftle( T, J )
% Compute FTLE using the mesochronic Jacobian J.
%
% 

validateattributes(J, {'numeric'},  {});
validateattributes(T, {'numeric'}, {'positive'});

J_Phi = eye(size(J)) + T*J; % convert mesochronic Jacobian to flow map Jacobian

ftle_val = log(norm(J_Phi,2)) / abs(T);
