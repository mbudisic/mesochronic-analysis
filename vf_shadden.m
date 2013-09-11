function f = vf_shadden(t,p)
% VF_SHADDEN(t,p)
%
% Running the function without arguments will run a showcase
% demonstration of the flow, showing vector field arrows overlaying the
% stream function.
%
% Vectorized ODE model of the double-gyre with moving boundary, used by
% Shawn Shadden as an example in his LCS tutorial
% t - 1xN time vector
% p - 2xN matrix of states (each column is a state)
%
%
% decent range for p is [0,2] x [0,1] box
% decent range for t is  [0,15] days
%
%
% For detailed parameter set see subfunction vf in source code, and
% website
% http://mmae.iit.edu/shadden/LCS-tutorial
%
% or just google for "Shadden LCS tutorial" and you'll find it.

if nargin < 1
    Nq = 30;
    [Xq,Yq] = meshgrid(linspace(0,2,Nq),linspace(0,1,Nq)); 
    Nls = 200;
    [Xls,Yls] = meshgrid(linspace(0,2,Nls),linspace(0,1,Nls)); 
    
    for t = linspace(0,10,100)
        [~,dPsidx, dPsidy] = vf(Xq(:).',Yq(:).',t);
        [Psi,~,~] = vf(Xls(:).',Yls(:).',t);
        contourf(Xls, Yls, reshape(Psi, Nls,Nls));
        hold all
        quiver(Xq, Yq, -reshape(dPsidy, Nq,Nq), reshape(dPsidx, Nq,Nq),'w')
        xlabel('x')
        ylabel('y')
        title(sprintf('Shadden double-gyre; Color is the stream function; T = %.2f ',t))
        hold off;
        caxis([-0.1,0.1]);
        pause(0.1);
    end
else
    [~,dPsidx, dPsidy] = vf(p(1,:),p(2,:),t);
    f = [-dPsidy; dPsidx];
end


function [Psi,dPsidx, dPsidy] = vf(x,y,t)

A = 0.10;
omega = 2*pi/10;
epsilon = 0.25;  % .10; % magnitude of perturbation.


a = epsilon * sin( omega * t );
b = 1 - 2 * a;

f = a .* x.^2 + b .* x;
df = 2 * a .* x + b;


Psi = A*sin(pi*f).*sin(pi*y);
dPsidx = pi * A * cos( pi * f ) .* sin(pi*y) .* df;
dPsidy = pi*A*sin( pi * f ) .* cos(pi*y);
