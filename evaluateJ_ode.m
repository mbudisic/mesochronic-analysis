function [mJ, retstep, sol, order] = evaluateJ_ode( frame, order, ic, f, t0, T, direction, h, dp )
% EVALUATEJ_ODE( order, ic, f, t0, T, h, dp )
%
% Evaluate mesochronic Jacobian using ODE evolution.
%
% frame - 'reference' or 'frenet'
% order - order of the method used
% ic - a single initial condition (column vector)
% f - flow field (vectorized function handle)
% t0 - initial time (scalar)
% T - vector of integration lengths (positive)
% direction (nonzero) - positive for forward in time, negative for backward
% h - resampling time
% dp - finite difference used for evaluation of the instantaneous Jacobian

validateattributes(ic, {'numeric'}, {'column'})
validateattributes(t0, {'numeric'}, {'scalar','real'})
validateattributes(direction, {'numeric'}, {'scalar', 'real', 'nonzero'});
validateattributes(T, {'numeric'}, {'positive'});

T = sort(T);
assert( min(diff([t0;t0+T(:)])) > h, 'return steps should differ by more than integration step')
assert( h <= min(T), 'time step should be smaller than integration length' )
assert( dp < 1, 'dp should be small')

Tmax = max(T);
t = t0 + ( 0:h:Tmax );

%% simulate in reference frame

opts = odeset('vectorized','on','RelTol',1e-4,'AbsTol',1e-6); % tolerances

% direction of time arrow
if direction < 0 
    fsim = @(t,x)(-f(-t,x));
else
    fsim = f;
end

% simulation
[tout, yout] = ode23t(fsim, t, ic, opts);
sol.t = tout;
sol.x = yout;

% make sure we obtained results at required times
assert( norm( t - tout.', Inf) < 1e-12, 'Output times do not match input times');

% evaluate instantaneous Jacobian using finite differences
[Ji, fi] = jacobian_fd(fsim, tout, yout, dp);
sol.Ji = Ji;
sol.fi = fi;

% return steps are multiples of resampling interval
retstep = fix(T/h);

% solve the Jacobian ODE
try
    [mJ, steps, order] = mcjacobian_mex(h, Ji, retstep, order); 
catch
    
    if verLessThan('matlab', '8')
        warning('MATLAB:MesochronicAnalysis:Compile', ...
            sprintf(['Cannot run compiled version of code on releases older than 2012b.\n',...
            'Running uncompiled code.\n',...
            'Turn this warning off by issuing: warning(''off'', ''MATLAB:MesochronicAnalysis:Compile'')']));
    else
        warning('MATLAB:MesochronicAnalysis:Compile', ...
            sprintf(['Compiled version does not exist.\nPlease run ''deploytool -build mcproject.prj'' prior to running for faster computations.\n',...
            'Running uncompiled code.\n',...
            'Turn this warning off by issuing: warning(''off'', ''MATLAB:MesochronicAnalysis:Compile'')']));
    end
    
    
    [mJ, steps, order] = mcjacobian(h, Ji, retstep, order);
end

% note: steps(k) indexes elements in tout, Ji, fi vectors (instantaneous time steps)
%       k        indexes elements in mJ, T (steps of averaging times requested)

%% Change to Frenet frame if needed
switch lower(frame)
    case 'frenet' % modify jacobians by frenet frame if needed
        
        if size(mJ,1) ~= 2
            error('Frenet frame is currently available only for 2d flows')
        end
        
        cframe = @(f)[f, [0 -1; 1 0]*f]/norm(f); % construct frenet frame from vf f
        f0 = fi(:,1);
        if norm(f0) > 1e-16 % if it's a zero-vector field, skip
            for k = 1:numel(steps)
                ft = sol.fi(:,steps(k)); % select the vector field
                Lfrenet_t = cframe(ft); % form the frenet frames
                Lfrenet_0 = cframe(f0);
                Jf = mJ(:,:,k);     % select the mesochronic jacobian
                
                % change coordinate systems:
                Tdelta = T(k);
                if Tdelta ~= 0
                    mJ(:,:,k) = (Lfrenet_t.' * Lfrenet_0 - eye(2))/Tdelta ...
                        + Lfrenet_t.' * Jf * Lfrenet_0;
                else
                    mJ(:,:,k) =   Lfrenet_t.' * Jf * Lfrenet_0;
                end
            end
        end
end

