function retval = meh3d_simulation(f, T, method, ics, h, dp, order, name)
% meh3d_simulation(f, T, method, N, h, dp, order, name)
%
% Compute mesohyperbolicity for the vector field f. (2d on [-0.5,0.5]^2
% grid for now)
%
% f - vector field
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

filename = sprintf('system%s_%s_T_%.1f_N_%03d',name,method, max(T), N);

Dets = zeros([Npoints, length(T)]);
Traces = zeros(size(Dets));
Meh = zeros(size(Dets));
Compr = zeros(size(Dets));
Nml = zeros(size(Dets));
Defect = zeros( size(Dets));
Jacobians = cell(Npoints,1);

if Ndim == 3
    TrCof = zeros(size(Dets));
end

for k = 1:Npoints
    ic = ics(k, :).';
    
    if strcmpi(method,'ode')
        mJ = evaluateJ_ode( order, ic, f, T, h, dp );
    else
        mJ = evaluateJ_fd( ic, f, T, h, dp );
    end
    
    myDets = zeros(1,length(T));
    myTraces = zeros(size(myDets));
    myMeh = zeros(size(myDets));
    myCompr = zeros(size(myDets));
    myNml = zeros(size(myDets));
    myDefect = zeros(size(myDets));
    
    if Ndim == 3
        myTrCof = zeros(size(myDets));
    end
    
    for n = 1:length(T)
        
        [classes, compr, spectral] = meh3d( T(n), {mJ(:,:,n)} );
        
        myDets(n) = spectral.dets;
        myTraces(n) = spectral.traces;
        myMeh(n) = classes;
        myCompr(n) = compr;
        myNml(n) = spectral.nml;
        myDefect(n) = spectral.defect;
        
        if Ndim == 3
            myTrCof(n) = spectral.trcofs;
        end
        
    end
    
    Jacobians{k} = mJ;
    Dets(k,:) = myDets(:);
    Traces(k,:) = myTraces(:);
    Meh(k,:) = myMeh(:);
    Compr(k,:) = myCompr(:);
    Nml(k,:) = myNml(:);
    Defect(k,:) = myDefect(:);
    
    if Ndim == 3
        TrCof(k,:) = myTrCof(:);
    end
    
    
end
disp('All done');

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



save([filename '.mat'], '-struct','retval')

