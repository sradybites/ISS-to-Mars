function [t,X,Y,Z,VX,VY,VZ] = newton_2body_proj(a,e,mu,I,w,O,tp,dt)
 % Two-body numerical integration
 % Inputs:
 %   a - semi-major axis
 %   e - eccentricity
 %   mu - gravitational parameter (with same distance units as a)
 %
 % Output:
 %   t - 1000 element time array between 0 and the orbital period
 %   x,y - Orbital radius components for one orbit in the perifocal
 %         frame (e,q) directions. Same size as t.
 %   En - Specific orbital energy over one orbital period. Same
 %       size as t.
 
 h = sqrt(mu*a*(1-e^2)); %angular momentum
 Tp = 2*pi*sqrt(a^3/mu); %orbital per/iod
 t = linspace(0,dt,1e3).'; %time array
 
 %initial conditions (at periapse):
 r0 = a*(1-e);
 x0 = r0;
 y0 = 0;
 vx0 = 0;
 vy0 = mu/h*(e+1);
 z0 = [x0,y0,vx0,vy0].';
 
 [~,z] = ode45(@newtongravity2d,t,z0);
 x = z(:,1);
 y = z(:,2);
 vx = z(:,3);
 vy = z(:,4);
 
 R = [x.';y.';zeros(1,length(x))];
 dcms = DCMs();
 helio2peri = @(Omega,inc,omega) dcms{3}(omega) * dcms{1}(inc) * dcms{3}(Omega);
 mat = helio2peri(O,I,w);
 R = cell2mat( arrayfun(@(i) mat\R(:,i), 1:size(R,2), 'uni',0) );

 DRDT = [vx.';vy.';zeros(1,length(vx))];
 DRDT = cell2mat( arrayfun(@(i) mat\DRDT(:,i), 1:size(DRDT,2), 'uni',0) );

 %force all outputs to be col vectors
 X = R(1,:).';
 Y = R(2,:).';
 Z = R(3,:).';
 VX = DRDT(1,:).';
 VY = DRDT(2,:).';
 VZ = DRDT(3,:).';
 
     %integrator function
     %state vector is z = [x,y,dx,dy]
     function dz = newtongravity2d(~,z)
         rn3 = sqrt(z(1)^2+z(2)^2)^3;
         dz = [z(3);...
             z(4);...
             -mu*z(1)/rn3;...
             -mu*z(2)/rn3];
     end
 end