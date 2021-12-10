function [t,X,Y,Z,VX,VY,VZ,points] = kepler_2body_proj(a,e,mu,I,w,O,tp,dt)
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
 
 tol = 2*eps(2*pi); %use this for the tolerance in your Newton-Raphson
 
 n = sqrt(mu/a^3);       %mean motion
 h = sqrt(mu*a*(1-e^2)); %angular momentum
 Tp = 2*pi*sqrt(a^3/mu); %orbital period
 
 %create time array from 0 to To with 1000 total points
 t = linspace(0,dt,dt+1).';
 
 %calculate the true anomaly at each of the values in t
 M = n*(t-tp); %mean anomaly
 
 %Newton-Raphson to find eccentric anomaly
 %initialize
 counter = 0; %use counter to avoid infinite loops
 del = 1;
 E = M./(1-e);
 inds = E > sqrt(6*(1-e)./e);
 E(inds) = (6*M(inds)./e).^(1/3);
 
 %iterate
 while ((del > tol) && (counter < 1000))
     update = (M - E + e.*sin(E))./(e.*cos(E)-1);
     E = E - update;
     
     %option 1: check tolerance against difference between current estimate
     %of M and true M
     del = max(abs(M - (E - e.*sin(E))));
     
     %option 2: check tolerance against update:
     %del = max(abs(update));
  
     counter = counter+1; %increment counter
 end
 
 %solve for true anomaly and orbital radius magnitude
 nu = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
 r = a*(1-e^2)./(1+e*cos(nu));
 r_p = [r.*cos(nu),r.*sin(nu),zeros(size(r))].';

 
 % calcualte the x and y velocities in the perifocal frame
 vx = -mu/h*sin(nu);
 vy = mu/h*(e+cos(nu));
 
 dcms = DCMs();
 helio2peri = @(Omega,inc,omega) dcms{3}(omega) * dcms{1}(inc) * dcms{3}(Omega);
 iCp = helio2peri(O,I,w).';
 r_helio = iCp * r_p;
 
 DRDT = [vx.';vy.';zeros(1,length(vx))];
 DRDT = cell2mat( arrayfun(@(i) iCp*DRDT(:,i), 1:size(DRDT,2), 'uni',0) );

 %force all outputs to be col vectors
 X = r_helio(1,:).';
 Y = r_helio(2,:).';
 Z = r_helio(3,:).';
 VX = DRDT(1,:).';
 VY = DRDT(2,:).';
 VZ = DRDT(3,:).';
 
 % troubleshooting lambert solver
 nu0 = 0.1324;
 dnu1 = 2.7932;
 dnu2 = 3.3555;
 r = @(nu) a.*(1-e^2)./(1+e.*cos(nu));
 r_p =  @(nu) [r(nu).*cos(nu),r(nu).*sin(nu),zeros(size(r(nu)))].';
 r0_p = r_p(nu0);
 rf1_p = r_p(nu0+dnu1);
 rf2_p = r_p(nu0+dnu2);
 
 r0_helio = iCp * r0_p;
 rf1_helio = iCp * rf1_p;
 rf2_helio = iCp * rf2_p;
 points = {r0_helio, rf1_helio, rf2_helio};
 
 end