function retval = meh_simulation(f, t0, T, direction, method, ics, h, dp, order, tol, filetag)
% meh_simulation(f, t0, T, direction, method, ics, h, dp, order, tol,  name)
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
% T - vector of *positive* integration periods
% direction - direction of integration; positive for forward, negative for
%             backward
% method - method of evaluating mesochronic Jacobian
%          'ode' (Adams-Bashforth evolution) or 'fd' (finite difference)
% ics - Npoints x D list of initial conditions, where D is dimension of state
%       space (2 or 3)
% h - discretization of time
% dp - finite difference variation (for instantaneous Jacobian evaluation)
% order - order of Adams-Bashforth
%         -1 for highest possible Adams-Bashforth method (currently 6), 
%          1-6 for appropriate order of Adams-Bashforth
% tol - tolerance on evaluating zero-matching criteria (currently only for
%       3d)
% filetag (optional) - tag of file to be saved to
%
% returns: structure saved to the file.
%
% The output filename is generated from the filetag and basic
% information about the initial conditions simulated.
% If the filename of the output is found in the current directory
% the code will attempt to load the Jacobian data from it,
% therefore shortening the time required to re-analyze the system.
%
% Open matlabpool before running if parallel computation is desired.

validateattributes(f, {'function_handle'},{})
validateattributes(filetag, {'char'},{})
validateattributes(tol, {'numeric'},{'positive'})
validateattributes(direction, {'numeric'}, {'scalar', 'real', 'nonzero'});

direction = sign(direction);

fprintf(1, 'Running vector field %s\n in %+d direction.\n', func2str(f), direction);

if strcmpi(method,'ode')
    disp('Using ODE evolution');
elseif strcmpi(method, 'fd')
    disp('Using finite difference')
else
    error('mytype has to be either ODE or FD')
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

Npoints = size(ics,1);

if direction > 0
    dirlabel = 'fwd';
else
    dirlabel = 'bwd';
end

filename = sprintf('mcan_%s_%s_%sT_%.1f_N_%05d.mat',filetag,method, dirlabel, max(T), Npoints);

fileexists = exist(filename, 'file');

if fileexists
    disp([filename ' exists. Loading from it.']);
    retval = load(filename);
else
    retval = struct;
end

t1 = tic;

% Jacobian data is found in the current directory
if isfield(retval, 'Jacobians') && isfield(retval,'ics')
    
    if (size(retval.Jacobians{1}, 1) == Ndim) && all( ics(:) == retval.ics(:) )
        disp([filename ' has Jacobian data. Skipping computation']);
        Jacobians = retval.Jacobians;
    else
        error('Jacobian data does not match input data. Exiting');
    end
    
% no Jacobian data found - run the simulations    
else
    
    retval.ics = ics;
    retval.T  = T;
    retval.t0 = t0;
    retval.h = h;
    retval.dp = dp;
    retval.method = method;
    retval.order = order;
    retval.f = f;
    retval.tol = tol;
    retval.direction = direction;
    
    Jacobians = cell(Npoints,1);
    
    if matlabpool('size') < 2
        warning('Matlab running in serial mode. On multicore computers you can reduce computation time by opening parallel matlab jobs, e.g., run "matlabpool open"')
    end
    % first just compute Jacobians
    disp([filename ': Started Jacobian computation.']); pause(0.5);
    
    orderlist = nan(1,Npoints);
    
    parfor k = 1:Npoints
        ic = ics(k, :).';
        
        if strcmpi(method,'ode')
            [mJ,~,myorder] = evaluateJ_ode( order, ic, f, t0, T, direction, h, dp );
            orderlist(k) = myorder;
        else
            mJ = evaluateJ_fd( ic, f, T, h, dp );
        end
        
        Jacobians{k} = mJ;
        
    end
    
    % set the overall order to the lowest order used by computation
    retval.order = min(orderlist);
    fprintf(1, '%s : Jacobians computed in %.2f sec.\n', filename, toc(t1));
    pause(0.5);
    
    retval.Jacobians = Jacobians;

end

save(filename, '-struct','retval')
fprintf(1, '%s : Saved Jacobian data. Starting classification\n',filename);

% set up output storage structures
Dets = zeros([Npoints, length(T)]);
Traces = zeros(size(Dets));
Meh = zeros(size(Dets));
Compr = zeros(size(Dets));
NonNml = zeros(size(Dets));
NonDefect = zeros( size(Dets));
TrCof = zeros(size(Dets)); % Trace of Cofactor - used for 3d analysis
FTLE = zeros( size(Dets));
Hyp =  zeros( size(Dets));

if Ndim == 3
    disp('Using 3d mesohyperbolicity')
    assert(exist('meh3d','file'), '3d analysis not yet implemented');
    meh = @(T, J)meh3d( T, {J}, tol);
else
    disp('Using 2d mesohyperbolicity')
    meh = @(T,J)meh2d(T,{J});
end

% evaluate quantifiers for Jacobians
parfor m = 1:Npoints
    
    % temporary parallel job storage
    myDets = zeros(1,length(T));
    myTraces = zeros(size(myDets));
    myMeh = zeros(size(myDets));
    myCompr = zeros(size(myDets));
    myNonNml = zeros(size(myDets));
    myNonDefect = zeros(size(myDets));
    myFTLE = zeros(size(myDets));
    myHyp = zeros(size(myDets));
    
    myTrCof = zeros(size(myDets));
    myJacobians = Jacobians{m};
    
    % for each integration period - this could likely be automatized
    % within meh2d/meh3d functions
    % The faculty within meh2d/3d already exists, but we have to make
    % sure that it can handle a vector of Ts, appropriate slice
    % Jacobians, etc.
    for n = 1:length(T)
        
        [classes, quants, spectral] = meh( T(n), myJacobians(:,:,n) );
        
        % spectral quantities
        myDets(n) = spectral.Dets;
        myTraces(n) = spectral.Traces;
        if isfield(spectral,'TrCofs')
            myTrCof(n) = spectral.TrCofs;
        end
        
        % quantifiers of mesochronic analysis criteria
        myCompr(n) = quants.Compr;
        myNonNml(n) = quants.NonNml;
        myNonDefect(n) = quants.NonDefect;
        myFTLE(n) = quants.FTLE;
        myHyp(n) = quants.Hyp;
        
        % mesochronic classes
        myMeh(n) = classes;
        
    end
    
    % copy temporary to output storage
    Dets(m,:) = myDets(:);
    Traces(m,:) = myTraces(:);
    Meh(m,:) = myMeh(:);
    Compr(m,:) = myCompr(:);
    NonNml(m,:) = myNonNml(:);
    FTLE(m,:) = myFTLE(:);
    Hyp(m,:) = myHyp(:);
    
    NonDefect(m,:) = myNonDefect(:);
    TrCof(m,:) = myTrCof(:);
    
end
fprintf(1, '%s : done in %.2f sec.\n', filename, toc(t1));
pause(0.5);

% create and store output structures
retval.Dets = Dets;
retval.Traces = Traces;
if Ndim == 3
    retval.TrCof = TrCof;
end

retval.Meh = Meh;

retval.Compr = Compr;
retval.NonNml = NonNml;
retval.FTLE = FTLE;
retval.Hyp = Hyp;
retval.NonDefect = NonDefect;


save(filename, '-struct','retval')
fprintf(1, '%s : Saved classification data. All done.\n',filename);

