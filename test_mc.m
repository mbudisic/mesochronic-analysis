eps = 0.0;
if eps > 1e-10
    name = 'meh2d_perturbed';
else
    name = 'meh2d_unperturbed';
end

u = @(x,y)[ -sin(x).*cos(y); cos(x) .* sin(y) ];
up = @(x,y)[ cos(x).*sin(y); -sin(x).*cos(y) ];

J = @(x,y)[ -cos(x).*cos(y), sin(x).*sin(y);...
     -sin(x).*sin(y), cos(x).*cos(y) ];

Jp = @(x,y)[-sin(x).*sin(y), cos(x).*cos(y);...
     -cos(x).*cos(y), sin(x).*sin(y)];

f = @(t,x)u(x(1), x(2)) + eps*cos(2*pi*t)*up(x(1), x(2));
Jf = @(t,x)J(x(1), x(2)) + eps*cos(2*pi*t)*Jp(x(1),x(2));

N = 20;
grid1d = 2*pi*linspace(1/N,1,N)-1/(2*N);

[X1, X2] = meshgrid(grid1d, grid1d);

figure('Name','Fields')
subplot(1,2,1);
F = u(X1, X2);
quiver(X1, X2, F(1:N,:), F(N+1:end,:));
title('Unperturbed field');
axis([0,2*pi,0,2*pi]);

subplot(1,2,2);
G = up(X1, X2);
quiver(X1, X2, G(1:N,:), G(N+1:end,:));
title('Perturbation')
axis([0,2*pi,0,2*pi]);

ics = [X1(:), X2(:)];
Nic = size(ics,1);

%T = pi*10;
T = 10;
dt = 1e-2;
t = 0:dt:T;
tc = num2cell(t);

steps = 20;

Dets = zeros([size(X1), steps]);
Traces = zeros([size(X1), steps]);
for kr = 1:N
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
    fprintf(1,'Row %03d/%03d completed\n', kr,N);
end
disp('All done');
tind = length(t):-fix(length(t)/steps):0;
Times = t(tind(end-1:-1:1));
Times = Times( end-steps+1:end );

mh = zeros(size(Dets));
comprs = zeros(size(Dets));

for step = 1:steps
    [mhclasses, compr] = meh2d( Times(step), 0, ...
        Dets(:,:,step), ...
        Traces(:,:,step) );
    mh(:,:,step) = mhclasses;
    comprs(:,:,step) = compr;
end

save(name,'Dets','Traces','Times','steps','N','mh','comprs','X1','X2','eps')
