function [ret] = plotPlanet(x,y,z,c)
    R_p = 0.06;
    [X_p,Y_p,Z_p] = sphere(50);
    X_p = X_p.*R_p + x;
    Y_p = Y_p.*R_p + y;
    Z_p = Z_p.*R_p + z;
    ret = surf(X_p,Y_p,Z_p,'FaceAlpha', 0.9, 'FaceColor', c,'EdgeColor','none');
end

