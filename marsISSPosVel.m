function [r_e,r_m,v_e,v_m] = marsISSPosVel(t0,dt,at,et,It,wt,Ot,tpt)

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

%Earth and Venus orbital elements
a_e = at;        %sma (AU)
a_m = 1.52371;
e_e = et;        %eccentricity
e_m = 0.09339;
I_e = It; %Inclination (rad)
I_m = 1.850*pi/180;%
w_e = wt;
w_m = 286.5*pi/180;
O_e = Ot;
O_m = 49.6*pi/180;
t_p_e = tpt;
t_p_m = 2454101.80556;
mu_e = 2.9591309705483544E-04;%G(m_sun + m_planet) (AU^3/day^2)
mu_m = 2.95912E-04;
te = (t0+t_p_e);
tv = (t0 + dt+t_p_m);
DCM_e = [cos(w_e),sin(w_e),0;-sin(w_e),cos(w_e),0;0,0,1]*[1,0,0;0,cos(I_e),sin(I_e);0,-sin(I_e),cos(I_e)]*[cos(O_e),sin(O_e),0;-sin(O_e),cos(O_e),0;0,0,1];
DCM_v =  [cos(w_m),sin(w_m),0;-sin(w_m),cos(w_m),0;0,0,1]*[1,0,0;0,cos(I_m),sin(I_m);0,-sin(I_m),cos(I_m)]*[cos(O_m),sin(O_m),0;-sin(O_m),cos(O_m),0;0,0,1];
%add any subfunctions as needed.  remember to put an 'end' statement at the end of any subfunctions


n_e = sqrt(mu_e/(a_e^3));
n_m = sqrt(mu_m/(a_m^3));
Me = mod(n_e*(t0-t_p_e),2*pi);
Mm = mod(n_m*(t0+dt-t_p_m),2*pi);

tol = 2*eps(2*pi);
if (Me/(1-e_e) < sqrt(6*(1-e_e)/e_e))
        E_e = (Me/(1-e_e));
    else
        E_e = (6*Me/e_e)^(1/3);
end
 check = true;
 j = 1;
    while (check) && j < 100
        E_e=E_e-((Me-E_e+e_e*sin(E_e))/(e_e*cos(E_e)-1));
        j = j + 1;
        err = (Me - (E_e - e_e*sin(E_e)));
        check = tol < abs(err);
    end

if (Mm/(1-e_m) < sqrt(6*(1-e_m)/e_m))
        E_v = (Mm/(1-e_m));
    else
        E_v = (6*Mm/e_m)^(1/3);
end
 check = 1;
 j = 1;
    while (check) && j < 100
        E_v=E_v-((Mm-E_v+e_m*sin(E_v))/(e_m*cos(E_v)-1));
        j = j + 1;
        err = (Mm - (E_v - e_m*sin(E_v)));
        check = tol < abs(err);
    end
E_e = mod(E_e,2*pi);
E_v = mod(E_v,2*pi);
    
nu_e = 2*atan(tan(E_e/2)*sqrt((1+e_e)/(1-e_e)));
if mod(E_e,2*pi) >= pi/2 && mod(E_e,2*pi) <= 3*pi/2
    nu_e = nu_e + pi;
end
nu_v = 2*atan(tan(E_v/2)*sqrt((1+e_m)/(1-e_m)));
if mod(E_v,2*pi) >= pi/2 && mod(E_v,2*pi) <= 3*pi/2
    nu_v = nu_v + pi;
end


nu_e = mod(nu_e,2*pi);
nu_v = mod(nu_v,2*pi);
b_e = sqrt(a_e^2 - (a_e*e_e)^2);
b_v = sqrt(a_m^2 - (a_m*e_m)^2);

r_m_val = a_e*(1-e_e^2)/(1+e_e*cos(nu_e));
%r_e_peri = r_e_val*[cos(nu_e),sin(nu_e),0];
r_e_peri = [a_e*(cos(E_e)-e_e),b_e*sin(E_e),0];
r_v_val = a_m*(1-e_m^2)/(1+e_m*cos(nu_v));
r_v_peri = r_v_val*[cos(nu_v),sin(nu_v),0];
r_e = transpose(DCM_e)*transpose(r_e_peri);
r_m = transpose(DCM_v)*transpose(r_v_peri);

Edot_e = n_e/(1-e_e*cos(E_e));
Edot_v = n_m/(1-e_m*cos(E_v));

v_e_peri = [-a_e*Edot_e*sin(E_e), b_e*Edot_e*cos(E_e),0];
v_v_peri = [-a_m*Edot_v*sin(E_v), b_v*Edot_v*cos(E_v),0];
v_e = transpose(DCM_e)*transpose(v_e_peri);
v_m = transpose(DCM_v)*transpose(v_v_peri);

%ensure outputs are column vectors
r_e = r_e(:);
r_m = r_m(:);
v_e = v_e(:);
v_m = v_m(:);
end