function [] = animateOrbit(filenames,t_orbit)
    h = figure
    axis tight manual
    filename = 'orbit.gif';
    
    earth_path = filenames(1);
    mars_path = filenames(2);
    [~,x_e,y_e,z_e] = readEphem(earth_path);
    [~,x_m,y_m,z_m] = readEphem(mars_path);
    
    % create sun
    R_sun = 0.004652 * 15; % 15AU
    [X_sun, Y_sun, Z_sun] = sphere(50);
    X_sun = X_sun.*R_sun;
    Y_sun = Y_sun.*R_sun;
    Z_sun = Z_sun.*R_sun;
    sun = surf(X_sun, Y_sun, Z_sun, 'FaceAlpha', 0.8, 'FaceColor', 'y');
    sun.EdgeAlpha = 0.05;
    
    for i = 1:dt
        %draw plots for planets
        plotPlanet(x_e(i),y_e(i),z_e(i));
        plotPlanet(x_m(i),y_m(i),z_m(i));
        drawnow
        
        %capture plot as image
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        %save image to gif
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end
end

