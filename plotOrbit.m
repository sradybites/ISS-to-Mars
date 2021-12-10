function [] = plotOrbit(filenames,t_orbit,planets, points)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    z_scale = 0.05;
    
    hold on
    for n = 1:length(filenames)
        fn = filenames(n);
        raw_ephem = fileread(fn);
        ephem = extractBetween(raw_ephem,'$$SOE','$$EOE');
        ephem = ephem{1}; % turn cell array to character vector

        %pull data from ephemeris readout (csv format)
        expr = '%f %*s %f %f %f %*f %*f %*f %*f %*f %*f';
        C = textscan(ephem,expr, 'Delimiter', ',');
        A = cell2mat(C).';
        t = A(1,:); X = A(2,:); Y = A(3,:); Z = A(4,:);

        % plot orbit of planet
        r0 = A(2:4,1);
        rf = A(2:4,end);
        nf = length(t);
        plot3(X,Y,Z,'Color','black');
        %plot3(r0(1),r0(2),r0(3),'o','MarkerEdgeColor','g','MarkerFaceColor','g')
        %plot3(rf(1),rf(2),rf(3),'o','MarkerEdgeColor','r','MarkerFaceColor','r')
        
        % fade in planets along trajectory!
        if planets
            dindex = floor((nf)/12);
            for i = 0:12
                j = 1 + dindex*i;
                R_p = 0.06;
                [X_p,Y_p,Z_p] = sphere(50);
                X_p = X_p.*R_p + X(j);
                Y_p = Y_p.*R_p + Y(j);
                Z_p = Z_p.*R_p.*z_scale + Z(j);

                p0 = surf(X_p,Y_p,Z_p,'FaceAlpha', ((j/nf)^10+0.05)/1.05);
                p0.EdgeColor= 'none';
            end
        end
    end
    
    % plot transfer orbit
    X_t = t_orbit{2}.';
    Y_t = t_orbit{3}.';
    Z_t = t_orbit{4}.';
    plot3(X_t,Y_t,Z_t,'Color','r','LineWidth',1);
    text(X_t(1),Y_t(1),Z_t(1),'SOM')
    text(X_t(end),Y_t(end),Z_t(end),'EOM')
    
    % create sun
    if planets
        R_sun = 0.004652 * 15; % 10^1 AU
        [X_sun, Y_sun, Z_sun] = sphere(50);
        X_sun = X_sun.*R_sun;
        Y_sun = Y_sun.*R_sun;
        Z_sun = Z_sun.*R_sun.*z_scale;

        sun = surf(X_sun, Y_sun, Z_sun, 'FaceAlpha', 0.8, 'FaceColor', 'y');
        sun.EdgeAlpha = 0.05;
    end
    
    % plot options
    title('ISS, Earth, and Mars during transfer')
    xlabel('x (AU)')
    ylabel('y (AU)')
    zlabel('z (AU)')
    
    grid on
    daspect([1 1 z_scale])
    view(3)
    
    % plot final velocity vector (x10)
    vx = t_orbit{5}; vy = t_orbit{6}; vz = t_orbit{7};
    v_f = 50.*[vx(end);vy(end);vz(end)];
    r_f = [X_t(end);Y_t(end);Z_t(end)];
    P1 = r_f + v_f;
    
    for i=1:length(points)
        r = points{i}
        plot3(r(1), r(2), r(3), "o", color="black");
    end
    
    hold off;
end

