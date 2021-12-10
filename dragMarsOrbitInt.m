function [t,r,v,te,ze] = dragMarsOrbitInt(r0,v0,tspan,CdAm) 
 % Numerical integration of Earth orbit with atmospheric drag to 200 kms
 %Inputs:
 %   r0 - 3x1 array: initial orbital position vector (km)
 %   v0 - 3x1 array: initial orbital velocity vector (km/s)
 %   tspan - array: input time array to ode solver; units of seconds
 %   CdAm - scalar: drag coefficient multiplied by spacecraft
 %                  cross-sectional area and scaled by the spacecraft mass
 %                  We assume this remains constant throughout the whole 
 %                  trajectory;  units of m^2/kg
 %
 %Outputs:
 %   t - (nx1) array: output time array from ode solver (s)
 %   r - (nx3) array: position trajectory elements (km)
 %   v - (nx3) array: velocity trajectory elements (km/s)
         
 %constants
 mu_m = 42828.37; %km^3/s^2  gravitational parameter of the Earth
 we = 7.292115e-5;    %rad/s     rotation rate of the Earth
 wm = we/1.026;       %rad/s     rotation rate of Mars
 Rm = 3389.5;      %km        radius of the Mars (assumed constant)
 
 %define the angular velocity of the Mars-fixed rotating frame:
 iWg = [0;0;wm]; 
 
 %initial state:
 z0 = [r0;v0].';
 
 %integrate:
 [t,z,te,ze] = ode113(@dragEarthOrbitInt_eom,tspan,z0,...
     odeset('AbsTol',1e-18,'RelTol',1e-20,...
     'Events',@dragEarthOrbitInt_event));
 
 r = z(:,1:3);
 v = z(:,4:6);
 
     %integrator function
     function dz = dragEarthOrbitInt_eom(~,z)
         %state vector is z = [x,y,z,dx,dy,dz]
         
         r = z(1:3); %position vector (km)
         v = z(4:6); %velocity vector (km/s)
         rmag = norm(r); %km
         h = rmag - Rm; %local altitude (assuming spherical Earth; km)
         
         %calculate the relative velocity assuming no winds:
         vrel = v - cross(iWg,r); %km/s
         vrelmag = norm(vrel);
         
         %calculate local atmospheric density and drag force
         rho = exponential_atmosphere_density(h); %kg/m^3
         CdAmrho = CdAm*rho; %this now has units 1/m
         CdAmrho = CdAmrho*1000; %1/m -> 1/km
         Fdrag = -CdAmrho/2*vrelmag*vrel; %km/s^2
         
         %acceleration is:
         a = -mu_m/rmag^3*r + Fdrag;
         dz = [v;a];
     end
 
     %exponential atmosphere function
     % see https://www.grc.nasa.gov/www/k-12/airplane/atmosmrm.html
     function  rho = exponential_atmosphere_density(h)
         %Inputs:
         %   h - height in km
         %
         %Outputs:
         %   rho - Atmospheric densities in kg/m^3
         alt = h*1000;
         if alt > 7000
             T = -23.4 - 0.00222*alt;
         else
             T = -31 - 0.000998*alt;
         end
         p = 0.699*exp(-0.00009*alt);
         rho = (p/(0.1921*(T+273.1)));
     end
 
     %crossing event function
     function [pos,isterm,dir] = dragEarthOrbitInt_event(~,z)
         r = z(1:3); %position vector (km)
         rmag = norm(r); %km
         h = rmag - Rm; %local altitude (assuming spherical Earth; km)
         pos = 450 - h; %zero crossing
         isterm = 0; %halt
         dir = 0;    %will always be from the top, but doesn't matter
     end
 end