function sims = meh2d_harness( yamlfile, names, dryrun )
% sims = meh2d_harness( yamlfile, names, dryrun )
%
% Reads a YAML configuration file and creates a set of
% mesochronic analysis simulations, by looping over all possible
% combinations of parameters given in yamlfile
%
% Master - struct representing parameters, e.g., ReadYAML('config.yaml')
% names - list of field names in the configuration file
%       - from which to generate simulation label
%       - can be empty - will be populated using any
%       - property that was a sequence

if nargin < 2
    names = {};
end

if ~exist('dryrun','var')
    dryrun = false;
else
    dryrun = true;
end

validateattributes( names, {'cell'}, {});
validateattributes( dryrun, {'logical'}, {})


sims = expandconf(ReadYaml(yamlfile), names );

fprintf(1,'Number of simulations expected: %d\n', numel(sims)) 

if dryrun
    for s = sims
        disp(s.filename)
        s.params
    end
    error('!!!Dry run!!! Pass ''false'' the third argument or skip it.')
else
    for s = sims
        disp(s.filename)
    end
end

for s = sims
    disp(['Running simulation:', s.filename]); tic;
    s.params
    
    %% simulation
            
    % determine the grid of inital conditions
    icgridx = linspace(s.params.minx,s.params.maxx,s.params.Nx);
    icgridy = linspace(s.params.miny,s.params.maxy,s.params.Ny);
    
    [X,Y] = meshgrid( icgridx, icgridy);
    ics = [X(:), Y(:)];
            
    % form the filename for saving the a Jacobians
    if strcmpi(s.params.direction,'fwd')
        direction = 1;
    else
        direction = -1;
    end
    
    % select the model
    switch lower(s.params.model)
        case {'rypina', 'polarjet'}
            vf = @(t,x)vf_rypina(t,x);
        case {'shadden', 'twogyre'}
            vf = @(t,x)vf_shadden(t,x);
        case { 'mezic', 'fourgyre'}
            vf = @(t,x)vf_mezic(t,x, s.params.epsilon);
        case { 'harmonic'}
            vf = @(t,x)vf_harmonic(t,x);
        case { 'nonlinear'}
            vf = @(t,x)vf_nonlin(t,x);
            
        otherwise
            error('Model not recognized')
    end
            
            
    jac_filename = [s.filename '.jac.mat'];
    meh_filename = [s.filename '.meh.mat'];
    yaml_filename = [s.filename '.yaml'];
    WriteYaml( yaml_filename, s.params );
    
    if exist(jac_filename, 'file')
        warning([jac_filename ' jacobian file already exists. Loading.']);
        Jdata = load(jac_filename);
    else
        disp(['Simulating ' jac_filename]);    
% f, t0, Ts, direction, method, ics, h, dp, order, tol        
        Jdata = meh_simulation(vf, s.params.t0, s.params.tdelta, direction, ...
                               s.params.method, ics, ...
                               s.params.tstep, s.params.sstep, ...
                               s.params.order, 1e-3);
        Jdata.params = s.params;
		Jdata.paramyaml = char(WriteYaml( [], s.params));
        save(jac_filename,'-struct', 'Jdata', '-v6');
    end
    
    % analyze Jacobian data using mesohyperbolic analysis
    if exist(meh_filename, 'file')
        warning([meh_filename ' mesochronic analysis file already exists.']);
    else
        MCdata = meh_analysis(Jdata.T, Jdata.Jacobians, Jdata.Ndim, 1e-6); %
        MCdata.params = s.params;

		% copying analyses to MCdata
		MCdata.T = Jdata.T;
		MCdata.hi_stretch = Jdata.hi_stretch;
		MCdata.hi_shear = Jdata.hi_shear;

		MCdata.paramyaml = char(WriteYaml( [], s.params));
        save(meh_filename,'-struct', 'MCdata', '-v6');
    end
    toc
end


end
%% AUXILIARY SECTION


function retval = expandconf( Master, names )
% EXPANDCONF
%
% Recursion (depth first)
%
% Expands array-based configuration to a set of non-arrayed
% configurations, for all combinations of parameters in arrays.

retval = [];

% expansion recursion
for fname = fieldnames(Master).' % pass through all properties
    
    myname = char(fname);
    if iscell( Master.(myname) ) % if property was read-in as array
        
        % skipped
        if strcmpi(myname, 'tdelta') % *** Tdelta is the only one passed as array, it is looped internally ***
            Master.(myname) = [Master.(myname){:}];
            continue
        end
        
        names =  [ names, {myname} ]; % add its name to labels list
        for k = Master.(myname) % for each property-array element
            
            % copy the master, but replace the property-array with the element
            Child = Master;                 
            Child.(myname) = cell2mat(k);
            
            childretval = expandconf(Child, names); % continue descending
            retval = [ retval(:).', childretval(:).' ]; % add all returned configs to higher level
        end
        return
    end
end

% if there were no property-arrays, formulate the filename

% concatenate required labels into one big label
labels = cellfun( @(n)repr(Master, char(n)), names, 'Uniformoutput',false);

% append the label to simulation name and return
retval.filename = [Master.simulationname , labels{:}];
retval.params = Master;

end

function retval = repr(S, name)
% string representation of the field 'name' in structure 'S'
myval = S.(name);
validateattributes(myval, {'numeric', 'char'}, {});
if ~ischar(myval)
    if isfloat(myval)
        myval = sprintf('%.1f', myval);
    else
        myval = sprintf('%d', myval);
    end
end
retval = sprintf('_%s_%s',name,myval);
end
