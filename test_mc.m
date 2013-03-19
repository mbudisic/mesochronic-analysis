if ~exist('eps','var')
    error('need to define eps perturbation')
else
    fprintf(1, 'Running flow with perturbation %.1f', eps);
end
name = sprintf('meh2d_eps_%.1f',eps);

u = @(x,y)[ -sin(x).*cos(y); cos(x) .* sin(y) ];
up = @(x,y)[ cos(x).*sin(y); -sin(x).*cos(y) ];

J = @(x,y)[ -cos(x).*cos(y), sin(x).*sin(y);...
     -sin(x).*sin(y), cos(x).*cos(y) ];

Jp = @(x,y)[-sin(x).*sin(y), cos(x).*cos(y);...
     -cos(x).*cos(y), sin(x).*sin(y)];

f = @(t,x)u(2*pi*x(1), 2*pi*x(2)) + eps*cos(2*pi*t)*up(2*pi*x(1), 2*pi*x(2));
Jf = @(t,x)2*pi*J(2*pi*x(1), 2*pi*x(2)) + 2*pi*eps*cos(2*pi*t)*Jp(2*pi*x(1), 2*pi*x(2));

N = 50;
grid1d = linspace(1/N,1,N)-1/(2*N);

[X1, X2] = meshgrid(grid1d, grid1d);

% %% plotting fields
% figure('Name','Fields')
% subplot(1,2,1);
% F = u(2*pi*X1, 2*pi*X2);
% quiver(X1, X2, F(1:N,:), F(N+1:end,:));
% title('Unperturbed field');
% axis([0,1,0,1]);
% 
% subplot(1,2,2);
% G = up(2*pi*X1, 2*pi*X2);
% quiver(X1, X2, G(1:N,:), G(N+1:end,:));
% title('Perturbation')
% axis([0,1,0,1]);

%% computation


ics = [X1(:), X2(:)];
Nic = size(ics,1);

%T = pi*10;
T = 5;
dt = 5e-3;
t = 0:dt:T;
tc = num2cell(t);

steps = 20;

Dets = zeros([size(X1), steps]);
Traces = zeros([size(X1), steps]);
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

save([name '.mat'],'Dets','Traces','Times','steps','N','mh','comprs','X1','X2','eps')
