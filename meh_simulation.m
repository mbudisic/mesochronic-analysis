function retval = meh_simulation(f, t0, Ts, direction, method, ics, h, dp, order, tol)
% retval = meh_simulation(f, t0, Ts, direction, method, ics, h, dp, order, tol)
%
% Compute mesochronic analysis of a dynamical system given by
% the vector field f under assumption of incompressibility.
%
% f - vector field (vectorized function handle)
%     Please use vectorized flow fields. In other words, if your
%     vector field is
%     f = @(t,x)[ -x(2); x(1) ]
%     its vectorized analogue would be
%     f = @(t,x)[ -x(2,:); x(1,:) ]
% t0 - initial time
% Ts - vector of *positive* integration periods
% direction - direction of integration; positive for forward, negative for
%             backward
% method - method of evaluating mesochronic Jacobian
%          'ode' (Adams-Bashforth evolution) or 'fd' (finite difference)
% ics - Npoints x D list of initial conditions, where D is dimension of state
%       space (2 or 3)
% h - discretization of time
% dp - finite difference variation (for instantaneous Jacobian evaluation, 1e-6 is a good guess)
% order - order of Adams-Bashforth
%         -1 for highest possible Adams-Bashforth method (currently 6),
%          1-6 for appropriate order of Adams-Bashforth
% tol - tolerance on evaluating zero-matching criteria (currently only for
%       3d)
%
% returns:
% retval - structure with parameters of the method
%          and a field Jacobians containing mesochronic Jacobians
%          computed at intervals contained in input Ts
%
%
% Open matlabpool before running if parallel computation is desired.

validateattributes(f, {'function_handle'},{})
validateattributes(tol, {'numeric'},{'positive'})
validateattributes(direction, {'numeric'}, {'scalar', 'real', 'nonzero'});
direction = sign(direction);

fprintf(1, 'Running vector field %s\n in %+d direction.\n', func2str(f), direction);

allowed.ode_ref_frame = {1, 'ODE evolution, reference (static) frame'};
allowed.ode_frenet_frame = {2, 'ODE evolution, Frenet (moving) frame'};
allowed.fd_ref_frame = {3, 'Finite Difference, reference (static) frame'};

try
    sel = allowed.(method);
    disp(sel{2});
catch
    error('Method: %s not allowed. Allowed methods: %s', method,...
        [char(fields(allowed)), repmat(' ',numel(fields(allowed)),1)].');
end

fprintf(1,'h = %.2e\n', h);

%% computation

% determine dimension of the state space
try
    validateattributes(ics, {'numeric'}, {'ncols',3})
    Ndim = 3;
catch
    validateattributes(ics, {'numeric'}, {'ncols',2})
    Ndim = 2;
end

% number of initial conditions
Npoints = size(ics,1);
Jacobians = cell(Npoints,1);

% storage for quantifiers other than mesochronic Jacobian
if Ndim == 2
    hi_stretch = zeros(Npoints,length(Ts));
    hi_shear = zeros(Npoints,length(Ts));
end

% output structure
retval = struct;
retval.ics = ics;
retval.T  = Ts;
retval.t0 = t0;
retval.h = h;
retval.dp = dp;
retval.method = method;
retval.order = order;
retval.f = f;
retval.tol = tol;
retval.direction = direction;
retval.Ndim = Ndim;

if matlabpool('size') < 2 && ~verLessThan('matlab', '8')
    warning('Matlab running in serial mode. On multicore computers you can reduce computation time by opening parallel matlab jobs, e.g., run "matlabpool open"')
end

% list of orders for individually computed points
orderlist = nan(1,Npoints);

% for every initial condition, evaluate mesochronic jacobian
for k = 1:Npoints
    ic = ics(k, :).';
    
    switch sel{1}
        
        case allowed.ode_ref_frame{1}
            [mJ, retstep, sol, myorder] = evaluateJ_ode( 'reference', order, ic, f, t0, Ts, direction, h, dp );
            orderlist(k) = myorder;
            
        case allowed.ode_frenet_frame{1}
            [mJ, retstep, sol, myorder] = evaluateJ_ode( 'frenet', order, ic, f, t0, Ts, direction, h, dp );
            orderlist(k) = myorder;
            
        case allowed.fd_ref_frame{1}
            orderlist(k) = 0;
            mJ = evaluateJ_fd( ic, f, Ts, h, dp );
            
    end
    
    % store the output
    Jacobians{k} = mJ;
    
    if Ndim == 2
        myJi = sol.Ji;
        myfi = sol.fi;
        [stretch, shear] = halleriacono_mex(h, myJi, myfi);
        
        hi_stretch(k,:) = stretch(retstep);
        hi_shear(k,:) = shear(retstep);
    end
    
end

% set the overall order to the lowest order used by computation
retval.order = min(orderlist);
retval.Jacobians = Jacobians;
if Ndim == 2
    retval.hi_stretch = hi_stretch;
    retval.hi_shear = hi_shear;
end
