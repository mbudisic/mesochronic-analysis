f = @(t,x)...
    [ -sin(x(1) - pi/2).*cos(x(2)-pi/2) ; ...
    cos(x(1) - pi/2).*sin(x(2)-pi/2)];

Jf = @(t,x)...
    [ -cos(x(1) - pi/2).*cos(x(2)-pi/2), sin(x(1) - pi/2).*sin(x(2)-pi/2); ...
    -sin(x(1) - pi/2).*sin(x(2)-pi/2), cos(x(1) - pi/2).*cos(x(2)-pi/2)];


N = 10;
grid1d = linspace(0,2*pi,N+1);
grid1d = grid1d(1:end-1);

[X1, X2] = meshgrid(grid1d, grid1d);
ics = [X1(:), X2(:)];
Nic = size(ics,1);

T = 10;
dt = 0.01;
t = 0:dt:T;
tc = num2cell(t);

Dets = zeros(size(X1));
Traces = zeros(size(X1));
for kr = 1:N
    parfor kc = 1:N
        ic = [ X1(1,kc), X2(kr,1)] ;
        
    % simulate
    S = ode23t(f, [0, T], ic.');
    
    % uniform resampling
    y = num2cell( deval(t, S), 1);
    
    % jacobians
    Ji = cellfun(Jf, tc, y , 'UniformOutput', false);

    %    [mJ, ti] = mcjacobian_mex(dt, cat(3,Ji{:}), 1000, 2);
    [mJ, ti] = mcjacobian_mex(dt, cat(3,Ji{:}), 0, 2);
    
    ts = t(ti+1);
    cP = squeeze(charpoly_sequence( mJ ));
    Dets(kr, kc) = cP(end);
    Traces(kr, kc) = -cP(end-1);
    fprintf(1,'Row %03d Col %03d completed\n', kr, kc );
    end
end

[mhclasses, compr] = meh2d( T, 0, Dets, Traces );