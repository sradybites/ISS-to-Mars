function [r_e,r_m,v_e,v_m] = EarthVenusPosVel(t0,dt)

%INPUTS
%   t0 - Initial time (JD)
%   dt - Cruise duration (days)
%
%OUTPUTS
%   r_e - 3x1 col vector - Heliocentric position of Earth (in AU) at t0
%   r_v - 3x1 col vector - Heliocentric position of Venus (in AU) at t0+dt
%   v_e - 3x1 col vector - Heliocentric-Inertial velocity of Earth 
%                          (in AU/day) at t0
%   v_v - 3x1 col vector - Heliocentric-Inertial velocity of Venus 
%                          (in AU/day) at t0+dt

m2au = 6.6845871226706e-12; % AU/m

%Earth and Venus orbital elements
a_e = 1.000373836656026E+00;        %sma (AU)
a_m = 7.233304294142561E-01;
e_e = 1.712127710968187E-02;        %eccentricity
e_m = 6.801883435239245E-03;
I_e = 2.777040607882003E-03*pi/180; %Inclination (rad)
I_m = 3.394490723804590E+00*pi/180; 
w_e = 3.043573249748720E+02*pi/180; %arg. of periapsis (rad)
w_m = 5.476413547432202E+01*pi/180;
O_e = 1.596967974767415E+02*pi/180; %long. of ascending node (rad)
O_m = 7.662845580893872E+01*pi/180;
t_p_e = 2458853.731945450883;       % time of periapsis passage (JD)
t_p_m = 2458029.748282369226; 
mu_e = 2.9591309705483544E-04;      %G(m_sun + m_planet) (AU^3/day^2)
mu_m = 2.9591293263082414E-04;


%add any subfunctions as needed.  remember to put an 'end' statement at the end of any subfunctions    
function [x,y,vx,vy] = kepler_2body(a,e,mu,tp,tf) 

tol = 2*eps(2*pi); %use this for the tolerance in your Newton-Raphson

n = sqrt(mu/a^3);       %mean motion
h = sqrt(mu*a*(1-e^2));       %angular momentum
Tp = 2*pi / n;       %orbital period

M = mod(n*(tf-tp),2*pi); %mean anomaly

%Newton-Raphson to find eccentric anomaly. Iterate to 
%machine-percision tolerance (tol)
iter = 0;
E_old = 0;
upper = sqrt(6*(1-e)/e);
lower = M/(1-e);
if lower < upper
    E = lower;
else
    E = (6*M/e)^(1/3);
end

while abs(E-E_old) >= tol && iter < 1000
    E_old = E;
    E = E - (M - E + e*sin(E))/(e*cos(E)-1);
    iter = iter + 1;
end

%calculate x and y positions in the perifocal frame using Kepler's
%equations
b = sqrt((1-e^2)*a^2);
x = a*(cos(E)-e);
y = b*sin(E);
r = hypot(x,y);

%calcualte the x and y velocities in the perifocal frame
vx = (a*n/r)*(-a*sin(E));
vy = (a*n/r)*(b*cos(E));

%force all outputs to be col vectors
end

rotMatE3 = @(ang) [cos(ang) sin(ang) 0;-sin(ang) cos(ang) 0;0 0 1];
rotMatE1 = @(ang) [1 0 0;0 cos(ang) sin(ang);0 -sin(ang) cos(ang)];
eRot313 = rotMatE3(w_e) * rotMatE1(I_e) * rotMatE3(O_e);
vRot313 = rotMatE3(w_m) * rotMatE1(I_m) * rotMatE3(O_m);

%calculate position and velocity vectors
[x_e, y_e, vx_e, vy_e] = kepler_2body(a_e, e_e, mu_e, t_p_e, t0);

[x_v, y_v, vx_v, vy_v] = kepler_2body(a_m, e_m, mu_m, t_p_m, t0+dt);


r_e = eRot313.' * [x_e;y_e;0];
v_e = eRot313.' * [vx_e;vy_e;0];

r_m = vRot313.' * [x_v;y_v;0];
v_m = vRot313.' * [vx_v;vy_v;0];

%ensure outputs are column vectors
r_e = r_e(:);
r_m = r_m(:);
v_e = v_e(:);
v_m = v_m(:);
end

