function f = vf_mezic(t,p, epsilon)
% VF_MEZIC(t,p)
%
% Running the function without arguments will run a showcase
% demonstration of the flow, showing vector field arrows overlaying the
% stream function.
%
% Vectorized ODE model of the four-gyre flow used by Mezic et al in the
% Science paper
% t - 1xN time vector
% p - 2xN matrix of states (each column is a state)
%
% decent range for p is [0,1] x [0,1] box
% decent range for t is  [0,15] days
% epsilon was 0.1 in the Science paper
%
% For detailed parameter set see subfunction vf in source code, and
%
% Mezi?, I., Loire, S., Fonoberov, V. A., & Hogan, P. J. (2010). 
% A New Mixing Diagnostic and Gulf Oil Spill Movement. Science Magazine, 330(6003), 
% 486?489. doi:10.1126/science.1194607

if nargin < 1
    Nq = 30;
    [Xq,Yq] = meshgrid(linspace(0,1,Nq),linspace(0,1,Nq)); 
    Nls = 200;
    [Xls,Yls] = meshgrid(linspace(0,1,Nls),linspace(0,1,Nls)); 
    epsilon = 1.1;
    for t = linspace(0,10,100)
        [~,dPsidx, dPsidy] = vf(Xq(:).',Yq(:).',t, epsilon);
        [Psi,~,~] = vf(Xls(:).',Yls(:).',t, epsilon);
        contourf(Xls, Yls, reshape(Psi, Nls,Nls));
        hold all
        quiver(Xq, Yq, -reshape(dPsidy, Nq,Nq), reshape(dPsidx, Nq,Nq),'w')
        xlabel('x')
        ylabel('y')
        title(sprintf('Mezic four-gyre; Color is the stream function; T = %.2f ',t))
        hold off;
        caxis([-0.15,0.15]);
        pause(0.1);
    end
else
    [~,dPsidx, dPsidy] = vf(p(1,:),p(2,:),t, epsilon);
    f = [-dPsidy; dPsidx];
end


function [Psi,dPsidx, dPsidy] = vf(x,y,t, epsilon)


PsiU = sin(2*pi*x).*sin(2*pi*y)/2/pi;
dPsiUdx = cos(2*pi*x).*sin(2*pi*y);
dPsiUdy = sin(2*pi*x).*cos(2*pi*y);

PsiP = cos(2*pi*x).*cos(2*pi*y)/2/pi;
dPsiPdx = -sin(2*pi*x).*cos(2*pi*y);
dPsiPdy = -cos(2*pi*x).*sin(2*pi*y);

Psi = PsiU + epsilon*cos(2*pi*t).*PsiP;
dPsidx = dPsiUdx + epsilon*cos(2*pi*t)*dPsiPdx;
dPsidy = dPsiUdy + epsilon*cos(2*pi*t)*dPsiPdy;
