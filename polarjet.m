function mydata = polarjet(mydata, T, N, direction)
% retval = polarjet(mydata, T, N, direction)
%
% A demo run for linear Bickley jet flow - simulation and visualization.
%
% If 'mydata' was passed, then the quantities in it are plotted for
% time period T (scalar).
%
% If 'mydata' is [] (empty), then a simulation is started
% using 'meh_simulation' file, results stored in a file
% whose name is output to Matlab window.
%
% input:
% T   - vector of integration times (positive values)
% N   - dimension of the 2d grid: number of initial conditions
%       per axis (total simulated is N^2)
% direction - direction of time (+1 or -1)
%
%
% output fields (D is dimension of state space):
%
%       ics - N^2 x dim initial conditions
%         T - integration times used (vector of length K)
%        t0 - inital time (scalar)
%         h - integration step used for trajectories (positive scalar)
%        dp - spatial step for evaluating inst. Jacobians using finite-difference
%    method - method used for evaluating mesohyperbolic Jacobian (ode or fd)
%     order - order of method used for evaluating mesohyperbolic Jacobian (if method = ode)
%         f - vector field simulated
%       tol - tolerance for zero matching criteria used
% direction - direction of time flow (1 for forward, -1 for backward)
% Jacobians - N^2-long cell array of mesochronic Jacobians,
%             each element is D x D x K or D x D x K, where K is length of vector T
%      Dets - determinants of mc. Jacobians
%             N^2 x K (columns correspond to elements of T)
%    Traces - traces of mc. Jacobians,
%             N^2 x K (columns correspond to elements of T)
%       Meh - mesohyperbolicity/mesoellipticity classes
%             N^2 x K (columns correspond to elements of T)
%     Compr - numerical compressibility (see below)
%             N^2 x K (columns correspond to elements of T)
%    NonNml - non-normality (see below) of mc. Jacobians
%             N^2 x K (columns correspond to elements of T)
%      FTLE - Finite Time Lyapunov Exponent (see meh2d.m and ftle.m)
%             N^2 x K (columns correspond to elements of T)
%       Hyp - hyperbolicity (see below) of mc. Jacobians
%             N^2 x K (columns correspond to elements of T)
% NonDefect - non-defectiveness (see below) of mc. Jacobians
%             N^2 x K (columns correspond to elements of T)
%
% Mesohyperbolicity class:
%    -1 - hyperbolicity, orientation preserving
%    0  - ellipticity
%    1  - hyperbolicity, orientation reversing
%
% Compressibility of mesochronic Jacobian J (theoretically identical to 0):
%   T * det J + tr J
%   This is a coarse measure of numerical error.
%
% Non-normality of mesochronic Jacobian J: Frobenius norm of the commutator
%   || J* . J - J . J*||
%   When non-normality is zero, J has complete orthogonal basis of eigenvectors
%   (it is *unitarily* diagonalizable)p.
%
% Hyperbolicity of mesochronic Jacobian J:
%   (T^2 * det J - 4) .* det J
%   When positive, J has a pair of real eigenvalues (flow map is hyperbolic for the integration time).
%
% Non-defectiveness of mesochronic Jacobian J:
%   smallest distance between roots of the minimal polynomial of J
%   When zero, J is defective, i.e., it has an
%   incomplete basis of eigenvectors (it is not diagonalizable)
%
%


%% SIMULATION
if  isempty(mydata)
    
    order = 2; % use 2nd order - safest before I check whether there are finite-precision errors in higher orders
    h = 1e-2; % uniform timestep
    dp = 1e-6; % spatial step for finite difference evaluation of inst. Jacobian
    tol = 1e-3; % tolerance on zero-matching criteria (irrelevant for 2d analysis)
    t0 = 0; % initial time
        
    % determine the grid of inital conditions
    icgridx = linspace(0,20,N);
    icgridy = linspace(-4,4,N);
    [X,Y] = meshgrid( icgridx, icgridy);
    ics = [X(:), Y(:)];
            
    % form the filename for saving the a Jacobians
    if direction > 0
        dirlab = 'fwd';
    else
        dirlab = 'bwd';
    end
    
    commonname = 'polarjet';
    
    filename = sprintf('%s_jac_o%d_N%d_%sT_%.1f.mat', commonname, order, N, dirlab, max(T));
    
    % if file exist, load Jacobian data
    if exist(filename,'file')
        disp(['Load ' filename]);
        Jdata = load(filename);
        
    % if file does not exist, simulate the system
    else
        disp(['Simulating ' filename]);        
        Jdata = meh_simulation(@(t,x)vf_polarjet(t,x), t0, T, direction, 'ode', ics, h, dp, order, tol);
        save(filename,'-struct', 'Jdata');
    end
    
    % analyze Jacobian data using mesohyperbolic analysis
    filename = sprintf('%s_meh_o%d_N%d_%sT_%.1f.mat', commonname, order, N, dirlab, max(T));
    MCdata = meh_analysis(T, Jdata.Jacobians, Jdata.Ndim, tol); %
    
    save(filename,'-struct', 'MCdata');
    
    % group evaluations of the Jacobian and mesochronic analysis into the
    % same data set
    for fname = fieldnames(Jdata).'
        mydata.(fname{1}) = Jdata.(fname{1});
    end
    
    for fname = fieldnames(MCdata).'
        mydata.(fname{1}) = MCdata.(fname{1});
    end
    
    % save joint data set
    filename = sprintf('%s_o%d_N%d_%sT_%.1f.mat', commonname, order, N, dirlab, max(T));
    save(filename,'-struct', 'mydata');
    
else

%% PLOTTING
    disp('Plotting the output.')
    
    % determine the initial condition grid from passed data
    Nic = size(mydata.ics, 1);
    N = fix(sqrt(Nic));
    assert(N == sqrt(Nic), 'Number of initial conditions is not a square of an integer. This is unsuitable for demo plotting');
    icgrid = linspace(0,1,N); % we use [0,2] grid for plotting only
    [X,Y] = meshgrid( icgrid, icgrid);
    
    % use value T if present in data, otherwise use maximum available T
    try
        validateattributes(T, {'numeric'}, {'scalar'})
    catch
        T = max(T);
    end
    if ~any(ismember(T, mydata.T))
        T = max(mydata.T);
        fprintf(1, 'Example of plotting for the max time T = %.1f \n', T);
    else
        fprintf(1, 'Plotting for the time T = %.1f \n',T);
    end
    
    % determine index of plotting T in the vector of available averaging
    % intervals mydata.T
    ind = find(T == mydata.T, 1, 'first');
    
    % all quantities are stored in matrices where
    % rows correspond to initial conditions
    % columns correspond to the averaging interval
    
    % once we select the column corresponding to requested T
    % we use the commands
    %    reshape( V, [N,N] )
    % to make it into a NxN matrix corresponding to NxN grid of initial
    % conditions, suitable for plotting
    invedges = [];
    % -- plotting different quantities --
    if mydata.direction > 0
        dirlabel = 'fwd';
    else
        dirlabel = 'bwd';
    end
    tstampline = sprintf(' for T = %.1f (%s)', T, dirlabel);
    % Mesochronic Classes
    n = 1;
    figure(n)
    pcolor(X,Y, reshape( signedlog10(mydata.Dets(:,ind)), [N,N]));
    setaxes;
    [cm, crange] = mehcolor(T, 64);
    
    colormap(cm); caxis([-crange, crange]);
    titleline = ['Mesochronic classes' tstampline];
    title(titleline)
    set(gcf,'name',titleline);
    
    set(gca, 'Color', 'black');
    if ~isempty(invedges)
        alpha(1-invedges)
    end
    
    cb = findobj(gcf,'tag','Colorbar');
    set(cb, 'YTick',[0, 4/(T^2)])
    set(cb, 'YTickLabel',{'0.0', '4/T^2'});
    %set(cb, 'YTickLabel',{'0.0', sprintf('%.f',4/T^2)});
    
    % Finite-Time Lyapunov Exponent
    if isfield(mydata,'FTLE')
        n = n+1; figure(n);
        pcolor(X,Y, reshape( mydata.FTLE(:,ind), [N,N]));
        setaxes(mydata.FTLE(:,ind));
        set(gca, 'Color', 'green');
        if ~isempty(invedges)
            alpha(1-invedges)
        end
        
        titleline =['FTLE' tstampline];
        title(titleline)
        set(gcf,'name',titleline);
        
        
    else
        disp('No FTLE field (Finite-Time Lyapunov Exponent) available')
    end
    
    % Deviation from a normal jacobian
    if isfield(mydata,'NonNml')
        n = n+1; figure(n);
        pcolor(X,Y, reshape( log10(mydata.NonNml(:,ind)), [N,N]));
        setaxes(log10(mydata.NonNml(:,ind)));
        cb = findobj(gcf,'tag','Colorbar');title(cb,'log_{10}')
        titleline = ['Non-normality' tstampline];
        title(titleline)
        set(gcf,'name',titleline);
        
    else
        disp('No NonNml field (deviation from normal Jacobian) available')
    end
    
    % Deviation from a defective jacobian
    if isfield(mydata,'NonDefect')
        n = n+1; figure(n);
        pcolor(X,Y, reshape( log10(mydata.NonDefect(:,ind)), [N,N]));
        setaxes(log10(mydata.NonDefect(:,ind)));
        map = colormap;
        colormap( map(end:-1:1, :) );
        cb = findobj(gcf,'tag','Colorbar');title(cb,'log_{10}')
        titleline= ['Non-defectiveness' tstampline];
        title(titleline)
        set(gcf,'name',titleline);
        
    else
        disp('No NonDefect field (deviation from defective Jacobian) available')
    end
    

    
    
    % Haller-Iacono shear
    if isfield(mydata,'hi_shear')
        n = n+1; figure(n);
        pcolor(X,Y, reshape( signedlog10(mydata.hi_shear(:,ind)), [N,N]));
        setaxes(signedlog10(mydata.hi_shear(:,ind)));
        colormap(diverging_map(linspace(0,1,64), [0.7,0,0],[0,0,0.7]))
        cb = findobj(gcf,'tag','Colorbar');title(cb,'sign(x)  log_{10}(1+|x|)')
        titleline=['Haller-Iacono shear' tstampline];
        title(titleline)
        set(gcf,'name',titleline);
        
    else
        disp('No hi_shear field (Haller-Iacono shear) available')
    end
    
    % Haller-Iacono stretch
    if isfield(mydata,'hi_stretch')
        n = n+1; figure(n);
        pcolor(X,Y, reshape( mydata.hi_stretch(:,ind), [N,N]));
        setaxes(mydata.hi_stretch(:,ind));
        colormap(diverging_map(linspace(0,1,64), [0.7,0,0],[0,0,0.7]))
        cb = findobj(gcf,'tag','Colorbar');
        titleline=['Haller-Iacono stretch' tstampline];
        title(titleline)
        set(gcf,'name',titleline);
        
    else
        disp('No hi_stretch field (Haller-Iacono stretch) available')
    end
    
end

function v = signedlog10( u )

v = log10( 1 + abs(u) ) .* sign(u);

%%
function retval = setaxes(fulldata)
% Helper function for setting axes appropriately
shading flat; axis([0,1, 0, 1]);
colorbar
try
    colormap morgenstemning
catch
    disp('Default color scheme "morgenstemning" missing. Download it from MATLAB Central (google: colormap morgenstemning). Using "hot" instead.');
    colormap hot
end
set(gca, 'XTick', (0:0.25:1))
set(gca, 'YTick', (0:0.25:1))
axis square
xlabel('x [\pi]')
ylabel('y [\pi]')

if exist('fulldata','var')
    caxis( prctile(fulldata(:), [1, 99]) );
    retval = prctile(fulldata(:), [1, 99]);
end

function retval = ridge( M, pct, sgn )
% Compute ridges/throughs using thresholding
%
% pct - percentage (e.g., 90 for ridge, 10 for through)
% sgn >= 0 -- ridge
% sgn < 0 -- through

if sgn >= 0
    retval = double(M > prctile( M(:), pct ));
else
    retval = double(M < prctile( M(:), pct ));
end
