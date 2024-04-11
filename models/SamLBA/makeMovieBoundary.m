function makeMovieBoundary(runs)

fps = 2;
frames = length(runs);
M = moviein(frames);
set(gca,'NextPlot','replacechildren');
for i=1:frames
        curPoly = runs{i}.polys2D{1};
        polyData = cell2mat(curPoly(2:end, :));     
        
        polyData(:,2:3) = 100.*polyData(:, 2:3);
        
        plot(polyData(:, 2), polyData(:, 3), 'o-');
        
        grid on;

    xlabel(curPoly{1,2});
    ylabel(curPoly{1,3});
    hold off;
    legend(['Flux limits ' num2str(runs{i}.bound)]);
    axis([-230 230 -10 10]);
    M(:,i) = getframe;
end
movie(M,3,fps)

end