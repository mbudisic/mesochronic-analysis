function uv = vf_bickley(t, p )
% VF_BICKLEY
%
% Time-invariant linear jet.

x = p(1,:);
y = p(2,:);

uv = [ sech(2*pi*x) .^ 2;...
       zeros(size(y)) ];
