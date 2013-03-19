function meh2dplot( X1, X2, Dets, T )
levels = 2048;
detrange = max(Dets(:));
fprintf(1, 'Dmax %.2f Ellimit %.2f\n', detrange, 4/(T^2));
[a, drange] = mehcolor(detrange, 4/(T^2), levels); 

pcolor(X1, X2, Dets ); shading flat;
caxis([-drange,drange]); 
colormap(a);
axis(2*pi*[0,1,0,1]);
colorbar

xlabel('x')
ylabel('y')
title(sprintf('T = %.2f', T));