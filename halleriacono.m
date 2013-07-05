function [stretch, shear] = halleriacono(h, Ji, fi)
% HALLERIACONO
%
% (steady flow only)
% 
% [stretch, shear] = halleriacono(h, Ji, stepsRequested)
%
% Computes 2D stretch and shear indicators described in 
% G. Haller and R. Iacono, Phys Rev E 68, 056304 (2003).
%           http://dx.doi.org/10.1103/PhysRevE.68.056304
%
% *** inputs: ***
% 
% h - length of a timestep (positive)
% Ji - array of instantaneous Jacobian matrices evaluated at times
%      uniformly spaced by h (each (:,:,k) is a Jacobian )
%
%    2 x 2 x N matrix, where 
%    N is the length of time evolution (maximum number of timesteps)
%
% 
%%                  
%
% *** outputs: ***
%
% stretch - 1 X N
%               "lambda" indicator, determining amount of material
%               stretching
%               N is the number of iteration steps returned
% shear - 1 X N
%               "mu" indicator, determining amount of material
%               stretching
%

assert( h > 0, 'Timestep has to be positive')
    
%% storage for intermediate results
sizeOfJacobians = size(Ji); % NxN
Nt = sizeOfJacobians(3);    % N of time steps

assert( all(sizeOfJacobians(1:2) == [2,2]), 'Jacobians should be 2x2 matrices');

Spar = zeros(1,Nt);
Scir = zeros(1,Nt);

Ids = [0 -1; 1 0]; % symplectic unity

fi_abs2 = sum( fi.^2,1); % magnitude of velocity squared

for k = 1:Nt
    f = fi(:,k);
    J = Ji(:,:,k);
    Jv = Ids*J;
    Spar(k) = f.' * J * f;
    Scir(k) = f.' * Jv * f;
end
Spar = Spar ./ fi_abs2;
Scir = Scir ./ fi_abs2;

stretch = - cumsum( Spar )*h;

lambda_bwd = - cumsum( Spar(end:-1:1) ) * h;
lambda_bwd = lambda_bwd(end:-1:1);

shear = sum( exp(-2*lambda_bwd) .* Scir )*h;

    