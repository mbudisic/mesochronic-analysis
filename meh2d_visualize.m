function handles = meh2d_visualize(Mdata, T)
% MEH2D_VISUALIZE
% handles = meh2d_visualize(Mdata, T)
% 
% Visualize the flow using data computed by meh2d_harness.
%
% Mdata - structure loaded from *.meh.mat file, containing all processed
% data.
% T - time interval requested (if omitted, max T available in Mdata is
%     used)
%
% If no output arguments are requested, the data is visualized in figures
% with handles 1, 2, etc.
%
% If output arguments are requested, figures are left hidden, and the
% handles are copied to higher numbers (so they are not overwritten by
% subsequent calls). 
% handles is a structure whose fields point to (hidden) figures showing
%         data named in field names.
%         Display figures by issuing, e.g.,
%             figure(handles.meh)

%% PLOTTING
disp('Plotting the output.')

% determine the initial condition grid from passed data
params = Mdata.params;

% determine the grid of inital conditions
icgridx = linspace(params.minx,params.maxx,params.Nx);
icgridy = linspace(params.miny,params.maxy,params.Ny);

[X,Y] = meshgrid( icgridx, icgridy);

% use value T if present in data, otherwise use maximum available T
try
    validateattributes(T, {'numeric'}, {'scalar'})
catch
    T = max(Mdata.T);
end
fprintf(1, 'Plotting for the time T = %.1f \n',T);


% determine index of plotting T in the vector of available averaging
% intervals mydata.T
ind = find(T == Mdata.T, 1, 'first');

% all quantities are stored in matrices where
% rows correspond to initial conditions
% columns correspond to the averaging interval

% once we select the column corresponding to requested T
% we use the commands
%    reshape( V, [params.Ny, params.Nx] )
% to make it into a NxN matrix corresponding to NxN grid of initial
% conditions, suitable for plotting
invedges = [];

% -- plotting different quantities --
tstampline = sprintf(' for t0 = %.2f, T = %.2f (%s)', params.t0, T, params.direction);

names = {};

% Mesochronic Classes
n = 1;
newfigure(n); names{n} = 'meh';

pcolor(X,Y,  reshape( Mdata.Dets(:,ind), [params.Ny,params.Nx]));
setaxes(params);
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
if isfield(Mdata,'FTLE')
    n = n+1; newfigure(n); names{n} = 'ftle';
    pcolor(X,Y, reshape( Mdata.FTLE(:,ind), [params.Ny, params.Nx]));
    setaxes(params, Mdata.FTLE(:,ind));
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
if isfield(Mdata,'NonNml')
    n = n+1; newfigure(n); names{n} = 'nonnml';
    pcolor(X,Y, reshape( Mdata.NonNml(:,ind), [params.Ny, params.Nx]));
    setaxes(params, Mdata.NonNml(:,ind));
    cb = findobj(gcf,'tag','Colorbar');
    titleline = ['Non-normality' tstampline];
    title(titleline)
    set(gcf,'name',titleline);
    
else
    disp('No NonNml field (deviation from normal Jacobian) available')
end

% Deviation from a defective jacobian
if isfield(Mdata,'NonDefect')
    n = n+1; newfigure(n); names{n} = 'nondefect';
    pcolor(X,Y, reshape( log10(Mdata.NonDefect(:,ind)), [params.Ny, params.Nx]));
    setaxes(params, log10(Mdata.NonDefect(:,ind)));
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
if isfield(Mdata,'hi_shear')
    n = n+1; newfigure(n);names{n} = 'hi_shear';
    pcolor(X,Y, reshape( signedlog10(Mdata.hi_shear(:,ind)), [params.Ny, params.Nx]));
    setaxes(params, signedlog10(Mdata.hi_shear(:,ind)));
    colormap(diverging_map(linspace(0,1,64), [0.7,0,0],[0,0,0.7]))
    cb = findobj(gcf,'tag','Colorbar');title(cb,'sign(x)  log_{10}(1+|x|)')
    titleline=['Haller-Iacono shear' tstampline];
    title(titleline)
    set(gcf,'name',titleline);
    
else
    disp('No hi_shear field (Haller-Iacono shear) available')
end

% Haller-Iacono stretch
if isfield(Mdata,'hi_stretch')
    n = n+1; newfigure(n);names{n} = 'hi_stretch';
    pcolor(X,Y, reshape( Mdata.hi_stretch(:,ind), [params.Ny, params.Nx]));
    setaxes(params, Mdata.hi_stretch(:,ind));
    colormap(diverging_map(linspace(0,1,64), [0.7,0,0],[0,0,0.7]))
    cb = findobj(gcf,'tag','Colorbar');
    titleline=['Haller-Iacono stretch' tstampline];
    title(titleline)
    set(gcf,'name',titleline);
    
else
    disp('No hi_stretch field (Haller-Iacono stretch) available')
end

% Numerical Compressibility (quantifies error in computation of
% Jacobian)
if isfield(Mdata,'Compr')
    n = n+1; newfigure(n);names{n} = 'compr';
    pcolor(X,Y, reshape( log10(abs(Mdata.Compr(:,ind))), [params.Ny, params.Nx]));
    setaxes(params, log10(abs(Mdata.Compr(:,ind))));
    cb = findobj(gcf,'tag','Colorbar');title(cb,'log_{10}')
    titleline=['Numerical compressibility' tstampline];
    title(titleline)
    set(gcf,'name',titleline);
    
else
    disp('No Compr field (numerical compressibility) available')
end


if nargout == 0
    for h = 1:n
        set(h, 'Visible','on');
    end
else
    for k = 1:numel(names)
        handles.(char(names{k})) = copyfig(k);
    end
end

end

function k = newfigure(k)
% create invisible figure
try
    set(0,'CurrentFigure',k);
    clf(k, 'reset');
catch
    figure(k);
end
set(k, 'visible','off');


end


function v = signedlog10( u )

v = log10( 1 + abs(u) ) .* sign(u);
end

%%
function retval = setaxes(params, fulldata)
% Helper function for setting axes appropriately
shading flat; axis([params.minx, params.maxx, params.miny, params.maxy]);
colorbar
try
    colormap morgenstemning
catch
    disp('Default color scheme "morgenstemning" missing. Download it from MATLAB Central (google: colormap morgenstemning). Using "hot" instead.');
    colormap hot
end
set(gca, 'XTick', linspace( params.minx, params.maxx,5 ) )
set(gca, 'YTick', linspace( params.miny, params.maxy,5 ) )
axis square
xlabel('x')
ylabel('y')

if exist('fulldata','var')
    caxis( prctile(fulldata(:), [1, 99]) );
    retval = prctile(fulldata(:), [1, 99]);
end
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
end

function overlay_ridge(axis_handle, mycolor, datamatrix)
% set background to desired color and alpha
% of the foreground to the datamatrix
set(axis_handle, 'Color', mycolor);
alpha(1-ridge( datamatrix, 90,1))
end
function E = getinvedges(N, filename)
% Retrieve edges of invariant sets

eqfield = load(filename);
[X,Y] = meshgrid(linspace(0,2*pi,N),linspace(0,2*pi,N));
E = double(edge(eqfield.F(X.',Y.'),'zerocross'));
end
