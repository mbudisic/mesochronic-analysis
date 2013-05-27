function uv = vf_fourgyre(t, p, epsilon )
% VF_FOURGYRE
%
% Vector field used by Mezic et al, Science 2010.
%
%
% epsilon - magnitude of perturbation (Science paper set it to 0.1)

x = p(1,:);
y = p(2,:);

uv = [...
    -sin(2*pi*x) .* cos(2*pi*y) + epsilon*cos(2*pi*t) .* cos(2*pi*x) .* sin(2*pi*y);...
    cos(2*pi*x) .* sin(2*pi*y) - epsilon*cos(2*pi*t) .* sin(2*pi*x) .* cos(2*pi*y) ];
