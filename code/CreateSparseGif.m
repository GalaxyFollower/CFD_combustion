function CreateSparseGif(T,X,Y,SkippedFrames,dimX,dimY,indexMap,nodeInfo)

figure (10)
filename = 'evolution_gif.gif';
ndof=dimX*dimY;
set(gca, 'FontSize', 18)
set(gcf,'position', [34   164   560   420]);
set(gcf,'paperUnits','centimeters','paperPosition',[0 0 10 8],'paperSize',[10 8])
for i = 1:SkippedFrames:length(T(1,:))
    
    hold off
    T_plot = reshape(mapToFull(T(:,i),indexMap,ndof),[dimY,dimX]);
    T_plot(find(nodeInfo<0))=NaN;
    surf(X,Y,T_plot);
    xlabel('X');
    ylabel('Y');
    title('Time evolution of acoustic waves');
    %campos([58.6464  -30  130.2245]);
    zlim([ min([min(T_plot), -0.001])  max([0.001, max(T_plot)])   ])%for P: -0.01 till 0.01; for u:-10,10
    drawnow
    frame = getframe (gcf);
    im=frame2im(frame);
    [imind,cm] = rgb2ind(im,8);
    if i==1
        imwrite(imind,cm,filename,'gif','loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','Writemode','append','Delaytime',0);
    end  
end