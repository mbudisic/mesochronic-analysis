function [retval, det_range] = mehcolor( det_range, ellimit, steps)


red = [1,0,0]*0.7;
green = [0,1,0]*0.7;
blue = [0,0,1]*0.7;
    

if det_range < ellimit
    warning('det_range will be modified')
    det_range = ellimit*1.5;
end

t = linspace(-det_range, det_range, steps);
el_upper = find( t  < ellimit, 1, 'last');
el_lower = find( t  > 0, 1, 'first');

sigmoid = @(nsteps)(1 - exp(-5*linspace(0,1,nsteps)));

t_blue = diverging_map( sigmoid(el_lower), [1,1,1], blue );
t_green = diverging_map( sigmoid(fix((el_upper - el_lower)/2)), [1,1,1], green );
t_red = diverging_map( sigmoid(steps-el_upper), [1,1,1], red );

retval = [ t_blue(end:-1:1,:); t_green; t_green(end:-1:1,:); t_red ];