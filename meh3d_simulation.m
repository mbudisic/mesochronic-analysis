function retval = meh3d_simulation(f, t0, T, method, ics, h, dp, order, name)
% meh3d_simulation(f, t0, T, method, N, h, dp, order, name)
%
% Compute mesohyperbolicity for the vector field f. (2d on [-0.5,0.5]^2
% grid for now)
%
% f - vector field
% t0 - initial time
% T - vector of periods of averaging
% mytype - 'ode' or 'fd' method of computation
% ics - Npoints x 3 list of initial conditions
% h - discretization of time
% dp - finite difference variation
% order - order of A-B method
% name (optional) - name of file to be saved to
%
% returns: structure saved to the file.
%
% Open matlabpool if parallel computation is desired.

%warning('Computing only for 2d fields for now')

validateattributes(f, {'function_handle'},{})
validateattributes(name, {'char'},{})

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
validateattributes(ics, {'numeric'}, {'ncols',3})
Npoints = size(ics,1);
Ndim = size(ics,2);

filename = sprintf('system%s_%s_T_%.1f_N_%05d.mat',name,method, max(T), Npoints);

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
    
    
    Jacobians = cell(Npoints,1);
    
    % first just compute Jacobians
    disp([filename ': Started Jacobian computation.']); pause(0.5);
    
    parfor k = 1:Npoints
        ic = ics(k, :).';
        
        if strcmpi(method,'ode')
            mJ = evaluateJ_ode( order, ic, f, t0, T, h, dp );
        else
            mJ = evaluateJ_fd( ic, f, T, h, dp );
        end
        
        Jacobians{k} = mJ;
        
    end
    
    fprintf(1, '%s : Jacobians computed in %.2f sec. Started classification.', filename, toc(t1));
    pause(0.5);
    
end

Dets = zeros([Npoints, length(T)]);
Traces = zeros(size(Dets));
Meh = zeros(size(Dets));
Compr = zeros(size(Dets));
Nml = zeros(size(Dets));
Defect = zeros( size(Dets));
TrCof = zeros(size(Dets));

% evaluate quantifiers for Jacobians
parfor k = 1:Npoints
    
    myDets = zeros(1,length(T));
    myTraces = zeros(size(myDets));
    myMeh = zeros(size(myDets));
    myCompr = zeros(size(myDets));
    myNml = zeros(size(myDets));
    myDefect = zeros(size(myDets));
    
    myTrCof = zeros(size(myDets));
    
    
    for n = 1:length(T)
        
        [classes, quants, spectral] = meh3d( T(n), {Jacobians{k}(:,:,n)}, 1e-3 );
        
        % spectral
        myDets(n) = spectral.Dets;
        myTraces(n) = spectral.Traces;
        if Ndim == 3
            myTrCof(n) = spectral.TrCofs;
        end
        
        % quants
        myCompr(n) = quants.Compr;
        myNml(n) = quants.Nml;
        myDefect(n) = quants.Defect;
        
        % classes
        myMeh(n) = classes;
        
    end
    
    Dets(k,:) = myDets(:);
    Traces(k,:) = myTraces(:);
    Meh(k,:) = myMeh(:);
    Compr(k,:) = myCompr(:);
    Nml(k,:) = myNml(:);
    Defect(k,:) = myDefect(:);
    TrCof(k,:) = myTrCof(:);
    
    
end
fprintf(1, '%s : done in %.2f sec.', filename, toc(t1));
pause(0.5);

retval.Jacobians = Jacobians;

retval.Dets = Dets;
retval.Traces = Traces;
retval.Meh = Meh;
retval.Compr = Compr;
retval.Nml = Nml;
retval.Defect = Defect;
retval.ics = ics;
retval.T  = T;
retval.t0 = t0;
retval.h = h;
retval.dp = dp;
retval.method = method;
retval.order = order;
retval.f = f;

if Ndim == 3
    retval.TrCof = TrCof;
end

save(filename, '-struct','retval')

