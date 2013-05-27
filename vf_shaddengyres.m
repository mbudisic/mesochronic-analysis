function uv = vf_shaddengyres( t, p )
% SHADDENGYRES
%
% Vector fields for the flow used by Shawn Shadden in his
% tutorial on Lagrangian Coherent Structures.
%

A = 0.10;
omega = 2*pi/10;
epsilon = 0.25;  % .10; % magnitude of perturbation.


x = p(1,:);
y = p(2,:);

a = epsilon * sin( omega * t );
b = 1 - 2 * a;

f = a .* x.^2 + b .* x;
df = 2 * a .* x + b;


u = -pi*A*sin( pi * f ) .* cos(pi*y);
v = pi * A * cos( pi * f ) .* sin(pi*y) .* df;

uv = [u; v];