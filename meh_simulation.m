function retval = meh_simulation(f, t0, T, direction, method, ics, h, dp, order, tol,  name)
% meh_simulation(f, t0, T, method, N, h, dp, order, name)
%
% Compute mesohyperbolicity for the vector field f. (2d on [-0.5,0.5]^2
% grid for now)
%
% f - vector field
% t0 - initial time
% T - vector of periods of averaging
% direction - direction of integration; positive for forward, negative for
%             backward
% mytype - 'ode' or 'fd' method of computation
% ics - Npoints x 3 list of initial conditions
% h - discretization of time
% dp - finite difference variation
% order - order of A-B method
% tol - tolerance on evaluating zero-matching criteria (currently only for
%       3d)
% direction
% name (optional) - name of file to be saved to
%
% returns: structure saved to the file.
%
% Open matlabpool if parallel computation is desired.

%warning('Computing only for 2d fields for now')

validateattributes(f, {'function_handle'},{})
validateattributes(name, {'char'},{})
validateattributes(tol, {'numeric'},{'positive'})
validateattributes(direction, {'numeric'}, {'scalar', 'real', 'nonzero'});

fprintf(1, 'Running vector field %s.\n', func2str(f));


if strcmpi(method,'ode')
    disp('Using ODE evolution');
elseif strcmpi(method, 'fd')
    disp('Using finite difference')
else
    error('mytype has to be either ODE or FD')
end

fprintf(1,'h = %.2e\n', h);

%% computation
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

filename = sprintf('mcan_%s_%s_%sT_%.1f_N_%05d.mat',name,method, dirlabel, max(T), Npoints);

fileexists = exist(filename, 'file');

if fileexists
    disp([filename ' exists. Loading from it.']);
    retval = load(filename);
else
    retval = struct;
end

t1 = tic;
if isfield(retval, 'Jacobians') && isfield(retval,'ics')
    
    if (size(retval.Jacobians{1}, 1) == Ndim) && all( ics(:) == retval.ics(:) )
        disp([filename ' has Jacobian data. Skipping computation']);
        Jacobians = retval.Jacobians;
    else
        error('Jacobian data does not match input data. Exiting');
    end
    
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
    
    parfor k = 1:Npoints
        ic = ics(k, :).';
        
        if strcmpi(method,'ode')
            mJ = evaluateJ_ode( order, ic, f, t0, T, direction, h, dp );
        else
            mJ = evaluateJ_fd( ic, f, T, h, dp );
        end
        
        Jacobians{k} = mJ;
        
    end
    
    fprintf(1, '%s : Jacobians computed in %.2f sec.\n', filename, toc(t1));
    pause(0.5);
    
    retval.Jacobians = Jacobians;

end

save(filename, '-struct','retval')
fprintf(1, '%s : Saved Jacobian data. Starting classification\n',filename);

Dets = zeros([Npoints, length(T)]);
Traces = zeros(size(Dets));
Meh = zeros(size(Dets));
Compr = zeros(size(Dets));
NonNml = zeros(size(Dets));
NonDefect = zeros( size(Dets));
TrCof = zeros(size(Dets));
FTLE = zeros( size(Dets));
Hyp =  zeros( size(Dets));

if Ndim == 3
    disp('Using 3d mesohyperbolicity')
    meh = @(T, J)meh3d( T, {J}, tol);
else
    disp('Using 2d mesohyperbolicity')
    meh = @(T,J)meh2d(T,{J});
end

% evaluate quantifiers for Jacobians
parfor m = 1:Npoints
    
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
    
    for n = 1:length(T)
        
        [classes, quants, spectral] = meh( T(n), myJacobians(:,:,n) );
        
        % spectral
        myDets(n) = spectral.Dets;
        myTraces(n) = spectral.Traces;
        if isfield(spectral,'TrCofs')
            myTrCof(n) = spectral.TrCofs;
        end
        
        % quants
        myCompr(n) = quants.Compr;
        myNonNml(n) = quants.NonNml;
        myNonDefect(n) = quants.NonDefect;
        myFTLE(n) = quants.FTLE;
        myHyp(n) = quants.Hyp;
        
        % classes
        myMeh(n) = classes;
        
    end
    
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

