function meh2dplot( X1, X2, Dets, T )
levels = 1024;
detrange = max(Dets(:));
fprintf(1, 'Dmax %.2f Ellimit %.2f\n', detrange, 4/(T^2));
[a, drange] = mehcolor(detrange, 4/(T^2), levels); 

pcolor(X1, X2, Dets ); shading flat;
caxis([-drange,drange]); 
colormap(a);
axis([0,1,0,1]);
colorbar

xlabel('X')
ylabel('Y')
title(sprintf('T = %.2f', T));