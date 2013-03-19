name = 'meh2d_unperturbed';
load(name)
close all; 

h1 = figure('Name','Det');
h2 = figure('Name','Class');

for k = 1:steps; 
    
    figure(h1); 
    meh2dplot(X1, X2, Dets(:,:,k), Times(k)),
    fprintf(1,'Det range: %.2e %.2e\n', min(min(Dets(:,:,k))), max(max(Dets(:,:,k))));
    
    figure(h2); 
    
    c = [0,0,1; 0,1,0; 1, 0, 0]*0.7;
    
    pcolor(X1, X2, mh(:,:,k))
    shading flat
    axis(2*pi*[0,1,0,1]);
    caxis([-1,1]);
    colormap(c)
    colorbar
    xlabel('x')
    ylabel('y')
    title(sprintf('T = %.2f', Times(k)));
    
    saveas(h1,sprintf('img/%s_%02d_det.png',name,k));
    saveas(h2,sprintf('img/%s_%02d_typ.png',name,k));
    
    pause(0.5);
end