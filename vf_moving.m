function f = vf_moving(t,p)
% VF_MOVING(t,p)
%
% Running the function without arguments will run a showcase
% demonstration of the flow, showing vector field arrows overlaying the
% stream function.
%
% Vectorized ODE model of a nonlinear traveling oscillator
% t - 1xN time vector
% p - 2xN matrix of states (each column is a state)
%
%
% decent range for p is [-1,1] x [-1,1] box
% decent range for t is  [0,5] 
%
%

if nargin < 1
    Nq = 30;
    [Xq,Yq] = meshgrid(linspace(-1,1,Nq)); 
    Nls = 200;
    [Xls,Yls] = meshgrid(linspace(-1,1,Nls)); 
    
    for t = linspace(0,10,100)
        [~,dPsidx, dPsidy] = vf(Xq(:).',Yq(:).',t);
        [Psi,~,~] = vf(Xls(:).',Yls(:).',t);
        contourf(Xls, Yls, reshape(Psi, Nls,Nls));
        hold all
        quiver(Xq, Yq, -reshape(dPsidy, Nq,Nq), reshape(dPsidx, Nq,Nq),'w')
        xlabel('x')
        ylabel('y')
        title(sprintf('Traveling gyre; Color is the stream function; T = %.2f ',t))
        hold off;
%        caxis([-0.1,0.1]);
        pause(0.1);
    end
else
    [~,dPsidx, dPsidy] = vf(p(1,:),p(2,:),t);
    f = [-dPsidy; dPsidx];
end


function [Psi,dPsidx, dPsidy] = vf(x,y,t)

Psi = sin(pi*(x - 0.2*t)).*sin(pi*y);
dPsidx = pi*cos(pi*(x - 0.2*t)).*sin(pi*y);
dPsidy = pi*sin(pi*(x - 0.2*t)).*cos(pi*y);
