offset = 0.5;
eps = 0.1;


f = @(t,x)...
    [ -sin(2*pi*(x(1) - offset)).*cos(2*pi*(x(2)-offset)) ;...
    cos(2*pi*(x(1) - offset)).*sin(2*pi*(x(2)-offset))];

Jf = @(t,x)2*pi*...
    [ -cos(2*pi*(x(1) - offset)).*cos(2*pi*(x(2)-offset)), sin(2*pi*(x(1) - offset)).*sin(2*pi*(x(2)-offset));...
    -sin(2*pi*(x(1) - offset)).*sin(2*pi*(x(2)-offset)), cos(2*pi*(x(1) - offset)).*cos(2*pi*(x(2)-offset))];


fp = @(t,x)eps*cos(2*pi*t).*...
    [ -sin(2*pi*(x(1) - offset)).*cos(2*pi*(x(2)-offset)) ;...
    cos(2*pi*(x(1) - offset)).*sin(2*pi*(x(2)-offset))];

Jfp = @(t,x)eps*cos(2*pi*t).*2*pi*...
    [ -cos(2*pi*(x(1) - offset)).*cos(2*pi*(x(2)-offset)), sin(2*pi*(x(1) - offset)).*sin(2*pi*(x(2)-offset));...
    -sin(2*pi*(x(1) - offset)).*sin(2*pi*(x(2)-offset)), cos(2*pi*(x(1) - offset)).*cos(2*pi*(x(2)-offset))];



N = 100;
grid1d = linspace(1/N,1,N)-1/(2*N);

[X1, X2] = meshgrid(grid1d, grid1d);
ics = [X1(:), X2(:)];
Nic = size(ics,1);

%T = pi*10;
T = 1;
dt = 1e-3;
t = 0:dt:T;
tc = num2cell(t);

steps = 10;

Dets = zeros([size(X1), steps]);
Traces = zeros([size(X1), steps]);
Times = zeros(steps,1);
parfor kr = 1:N
    for kc = 1:N
        ic = [ X1(1,kc), X2(kr,1)] ;
        
        % simulate
        S = ode23t(f, [0, T], ic.');
        
        % uniform resampling
        y = num2cell( deval(t, S), 1);
        
        % jacobians
        Ji = cellfun(Jf, tc, y , 'UniformOutput', false);
        
        %    [mJ, ti] = mcjacobian_mex(dt, cat(3,Ji{:}), 1000, 2);
        [mJ, ti] = mcjacobian_mex(dt, cat(3,Ji{:}), fix(length(t)/steps), 2);
        
        ts = t(ti+1);
        cP = squeeze(charpoly_sequence( mJ ));
        for step = 1:steps
            Dets(kr, kc, step) = cP(step,end);
            Traces(kr, kc, step) = -cP(step, end-1);
        end
    end
    fprintf(1,'Row %03d completed\n', kr);
end

tind = length(t):-fix(length(t)/steps):0;
Times = t(tind(end:-1:1));
Times = Times( end-steps+1:end );

mh = zeros(size(Dets));

% for step = 1:steps
% [mhclasses, compr] = meh2d( T, 0, Dets, Traces );
% mh(:,:,step) = mhclasses;
% end

save 'meh2d_rotation' Dets Traces Times steps N