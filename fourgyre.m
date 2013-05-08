function retval = fourgyre(mydata, T, N, tol)
% retval = fourgyre(mydata, T, N, tol)
%
% A demo run for four gyre flow - simulation and visualization.
%
% If 'mydata' was passed, then the quantities in it are plotted for
% time period T (scalar).
%
% If 'mydata' is empty, then a simulation is started
% using 'meh_simulation' file and stored in a file (and returned).
%
% T   - vector of integration times
% N   - number of initial points per axis (total simulated is N^2)
% tol - tolerance in determining mesohyperbolicity
%
%

%% SIMULATION
if  isempty(mydata)
    
    % determine the grid of inital conditions
    icgrid = linspace(0,2*pi,N);
    [X,Y] = meshgrid( icgrid, icgrid);
    ics = [X(:), Y(:)];
    
    % specification of the system
    epsilon = 0.1;
    f = @(t,x)[...
        -sin(x(1,:)) .* cos(x(2,:)) + 2*pi*epsilon*cos(2*pi*t) .* cos(x(1,:)) .* sin(x(2,:));...
         cos(x(1,:)) .* sin(x(2,:)) - 2*pi*epsilon*cos(2*pi*t) .* sin(x(1,:)) .* cos(x(2,:)) ];
    % run the simulation
    mydata = meh_simulation(f, 0, T, 'ode', ics, 1e-2, 1e-8, 2, tol, 'fourgyre');
    
    %% PLOTTING
end
retval = mydata;

if nargout == 0 
    disp('Plotting the output.')

    % determine the initial condition grid from passed data
    Nic = size(mydata.ics, 1);
    N = fix(sqrt(Nic));
    assert(N == sqrt(Nic), 'Number of initial conditions is not a square of an integer. This is unsuitable for demo plotting');
    icgrid = linspace(0,2,N); % we use [0,2] grid for plotting only
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
    tstampline = sprintf(' for T = %.1f', T);
    % Mesochronic Classes
    n = 1;
    figure(n)
    pcolor(X,Y, reshape( mydata.Dets(:,ind), [N,N]));
    setaxes
    [cm, crange] = mehcolor(T, 64);
    
    colormap(cm); caxis([-crange, crange]);
    title(['Mesochronic classes' tstampline])
    
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
        setaxes
        set(gca, 'Color', 'green');
        if ~isempty(invedges)
        alpha(1-invedges)
        end
        
        title(['FTLE' tstampline])
    else
        disp('No FTLE field (Finite-Time Lyapunov Exponent) available')
    end
    
    % Deviation from a normal jacobian
    if isfield(mydata,'NonNml')
        n = n+1; figure(n);
        pcolor(X,Y, reshape( log10(mydata.NonNml(:,ind)), [N,N]));
        setaxes(log10(mydata.NonNml))
        cb = findobj(gcf,'tag','Colorbar');title(cb,'log_{10}')
        title(['Non-normality' tstampline])
    else
        disp('No NonNml field (deviation from normal Jacobian) available')
    end
    
    % Deviation from a defective jacobian
    if isfield(mydata,'NonDefect')
        n = n+1; figure(n);
        pcolor(X,Y, reshape( log10(mydata.NonDefect(:,ind)), [N,N]));
        setaxes(log10(mydata.NonDefect))
        map = colormap;
        colormap( map(end:-1:1, :) );
        cb = findobj(gcf,'tag','Colorbar');title(cb,'log_{10}')
        title(['Non-defectiveness' tstampline])
    else
        disp('No NonDefect field (deviation from defective Jacobian) available')
    end
    
    
    % Hyperbolicity
    if isfield(mydata,'Hyp')
        
        % State space coloring
        n = n+1; figure(n);
        
        hypfield = sign(mydata.Hyp) .* log10( 1 + abs(mydata.Hyp) );
        
        
        
        pcolor(X,Y, reshape( hypfield(:,ind), [N,N]));
        setaxes
        cax = [-1,1]*max(abs(hypfield(:)));
        caxis(cax);
        sigmoid = @(nsteps)(1 - exp(-5*linspace(0,1,nsteps)));
        vars = sigmoid(32);
        vars = [-vars(end:-1:1).'; vars(:)];
        vars = ( vars - min(vars) ) / (max(vars) - min(vars));

        colormap(diverging_map(vars, [0,1,0], [1,0,0]))
        title(['Hyperbolicity' tstampline])
        
        % Histogram
        n = n+1; figure(n);
        boxes = linspace(min(hypfield(:)), max(hypfield(:)),50);
        hist(hypfield(:,ind), boxes); % histogram with persistent bounds
        
        pct = [10,50,90];
        pct_v = prctile(hypfield(:,ind),pct);
        
        ax = axis;
        hold on;
        for k = 1:length(pct)
            set( line( [pct_v(k),pct_v(k)], ax(3:4)), 'color','red');
            text( pct_v(k), ax(4) - 0.05*diff(ax(3:4)), sprintf('%d %%',pct(k)));
        end
        hold off;
        title(['Distribution of hyperbolicity' tstampline])
        xlabel('Hyperbolicity (log_{10})')
        ylabel('Number of initial conditions')
        set(gca, 'XTick', [ax(1), pct_v, ax(2)].')
        set(gca, 'XTickLabel', num2str([ax(1), pct_v, ax(2)].','%.1f'))
    else
        disp('No Hyp field (hyperbolicity quantifier) available')
    end
    
    % Numerical Compressibility (quantifies error in computation of
    % Jacobian)
    if isfield(mydata,'Compr')
        n = n+1; figure(n);
        pcolor(X,Y, reshape( log10(abs(mydata.Compr(:,ind))), [N,N]));
        setaxes(log10(abs(mydata.Compr)))
        cb = findobj(gcf,'tag','Colorbar');title(cb,'log_{10}')
        title(['Numerical compressibility' tstampline])
    else
        disp('No Compr field (numerical compressibility) available')
    end
    
end

%%
function setaxes(fulldata)
% Helper function for setting axes appropriately
shading flat; axis([0,2, 0, 2]);
colorbar
try
colormap morgenstemning
catch
    disp('Download "morgenstemning" color scheme from MATLAB Central (google: colormap morgenstemning) for grayscale-compatible colors ;)');
    colormap hot
end
set(gca, 'XTick', (0:0.25:1)*2)
set(gca, 'YTick', (0:0.25:1)*2)
axis square
xlabel('x [\pi]')
ylabel('y [\pi]')

if exist('fulldata','var')
     disp('Normalizing axes')
%     caxis( [min(fulldata(:)), max(fulldata(:))] );
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

function overlay_ridge(axis_handle, mycolor, datamatrix)
% set background to desired color and alpha
% of the foreground to the datamatrix
set(axis_handle, 'Color', mycolor);
alpha(1-ridge( datamatrix, 90,1))

function E = getinvedges(N, filename)
% Retrieve edges of invariant sets

eqfield = load(filename);
[X,Y] = meshgrid(linspace(0,2*pi,N),linspace(0,2*pi,N));
E = double(edge(eqfield.F(X.',Y.'),'zerocross'));

