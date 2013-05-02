function ftle_val = ftle( T, J )
% Compute FTLE using the mesochronic Jacobian J.
%
% 

validateattributes(J, {'numeric'},  {});
validateattributes(T, {'numeric'}, {'positive'});

J_Phi = eye(size(J)) + T*J; % convert mesochronic Jacobian to flow map Jacobian

% compute FTLE - (factor 2 accounts for eigenvalue of CG tensor being
%                 square of SV of Phi)
ftle_val = 2 * log10( max( svd( J_Phi ) ) ) / T;