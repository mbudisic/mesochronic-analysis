function [retval, det_range] = mehcolor( det_range, ellimit, steps)
% MEHCOLOR
%
% [retval, det_range] = mehcolor( det_range, ellimit, steps)
%
%
% det_range - minimum and maximum values of determinant of mesochronic
%              Jacobian visualized
% ellimit   - upper bound of mesoellipticity (4/T^2)
% steps     - number of gradations of color scheme

red = [1,0,0]*0.7;
green = [0,1,0]*0.7;
blue = [0,0,1]*0.7;
    
% ellimit is the upper bound of the mesoelliptic region
% color range has to at least include this limit
if det_range < ellimit
    warning('det_range will be modified')
    det_range = ellimit*1.5;
end

t = linspace(-det_range, det_range, steps);
el_upper = find( t  < ellimit, 1, 'last');
el_lower = find( t  > 0, 1, 'first');

% nonlinear variation of colors
sigmoid = @(nsteps)(1 - exp(-5*linspace(0,1,nsteps)));

% create the color map by stitching three zones of variation 
% (lower bound -> mesohyp/mesoell boundary (det = 0)
% mesohyp/mesoell (det = 0) -> mesoell/mesohyp boundary (det = ellimit)
% mesoell/mesohyp (det = ellimit) -> upper bound
t_blue = diverging_map( sigmoid(el_lower), [0.9,0.9,1], blue );
t_green = diverging_map( sigmoid(fix((el_upper - el_lower)/2)), [0.9,1,0.9], green );
t_red = diverging_map( sigmoid(steps-el_upper), [1,0.9,0.9], red );

retval = [ t_blue(end:-1:1,:); t_green; t_green(end:-1:1,:); t_red ];