function f = vf_polarjet(t,p)
% VF_POLARJET(t,p)
%
% Vectorized ODE model of the atmospheric zonal jet around south pole
% t - 1xN time vector
% p - 2xN matrix of states (each column is a state)
%
% timescale unit - days
% length unit - 10^6 m
%
% decent range for p is [0,20] x [-4,4] box (in Megameters)
% decent range for t is  [0,10] days
%
% Running the function without passing t, p will run a showcase
% demonstration of the flow, showing vector field arrows overlaying the
% stream function.
%
% For detailed parameter set see subfunction vf in source code, and
% paper
% I. Rypina, M.G. Brown, F.J. Beron-Vera, H. Koçak, M.J. Olascoaga, and I.A. Udovydchenkov, J Atmos Sci 64, 3595 (2007).
% "On the Lagrangian Dynamics of Atmospheric Zonal Jets and the
% Permeability of the Stratospheric Polar Vortex"
% doi: 10.1175/JAS4036.1


if nargin < 1
    Nq = 30;
    [Xq,Yq] = meshgrid(linspace(0,20,Nq),linspace(-4,4,Nq)); % unit - 10^6 m
    Nls = 200;
    [Xls,Yls] = meshgrid(linspace(0,20,Nls),linspace(-4,4,Nls)); % unit - 10^6 m
    
    for t = linspace(0,10,100)
        [~,dPsidx, dPsidy] = vf(Xq(:).',Yq(:).',t);
        [Psi,~,~] = vf(Xls(:).',Yls(:).',t);
        contourf(Xls, Yls, reshape(Psi, Nls,Nls));
        hold all
        quiver(Xq, Yq, -reshape(dPsidy, Nq,Nq), reshape(dPsidx, Nq,Nq))
        xlabel('Mm')
        ylabel('Mm')
        title(sprintf('T = %.2f days',t))
        hold off;
        caxis([-60,60]);
        pause(0.1);
    end
else
    [~,dPsidx, dPsidy] = vf(p(1,:),p(2,:),t);
    f = [-dPsidy; dPsidx];
end


function [Psi,dPsidx, dPsidy] = vf(x,y,t)

tscale = 24*3600/1e6;

U0 = 62.66;
re = 6.371; % unit - 10^6 m
L = 1.770; % unit - 10^6 m
k1 = 2/re;
k2 = 4/re;
k3 = 6/re;

c1 = 0;
c3 = 0.461*U0;
c2 = 0.205*U0;

A3 = 0.3;
A2 = 0.1;
A1 = 0.0;


Psi =  c3*y + U0*L*( -tanh(y/L) + (sech(y/L)).^2 .* (A3*cos(k3*(x-c3*t*tscale)) + A2*cos(k2*(x - c2*t*tscale)) + A1*cos(k1*(x - c1*t*tscale))) );
dPsidy = c3 + L.*U0.*((tanh(y/L).^2 - 1)/L - (2.*sinh(y/L).*(A1.*cos(k1.*(x - c1.*t*tscale)) + A2.*cos(k2.*(x - c2.*t*tscale)) + A3.*cos(k3.*(x - c3.*t*tscale))))/(L.*cosh(y/L).^3));
dPsidx = -(L.*U0.*(A1.*k1.*sin(k1.*(x - c1.*t*tscale)) + A2.*k2.*sin(k2.*(x - c2.*t*tscale)) + A3.*k3.*sin(k3.*(x - c3.*t*tscale))))./(cosh(y/L).^2);
