function [r_fin, v_fin] = rv_final_from_rv_initial(r,v)
%{
This program plots the ground track of an earth satellite
for which the orbital elements are specified
mu - gravitational parameter (km^3/s^2)
deg - factor that converts degrees to radians
J2 - second zonal harmonic
Re - earth’s radius (km)
we - earth’s angular velocity (rad/s)
rP - perigee of orbit (km)
rA - apogee of orbit (km)
TA, TAo - true anomaly, initial true anomaly of satellite (rad)
RA, RAo - right ascension, initial right ascension of the node (rad)
incl - orbit inclination (rad)
wp, wpo - argument of perigee, initial argument of perigee (rad)
n_periods - number of periods for which ground track is to be plotted
a - semimajor axis of orbit (km)
T - period of orbit (s)
e - eccentricity of orbit
h - angular momentum of orbit (km^2/s)
E, Eo - eccentric anomaly, initial eccentric anomaly (rad)
M, Mo - mean anomaly, initial mean anomaly (rad)
to, tf - initial and final times for the ground track (s)
fac - common factor in Equations 4.53 and 4.53
RAdot - rate of regression of the node (rad/s)
wpdot - rate of advance of perigee (rad/s)
times - times at which ground track is plotted (s)
ra - vector of right ascensions of the spacecraft (deg)
dec - vector of declinations of the spacecraft (deg)
TA - true anomaly (rad)
r - perifocal position vector of satellite (km)
R - geocentric equatorial position vector (km)
R1 - DCM for rotation about z through RA
R2 - DCM for rotation about x through incl
R3 - DCM for rotation about z through wp
QxX - DCM for rotation from perifocal to geocentric equatorial
Q - DCM for rotation from geocentric equatorial
into earth-fixed frame
r_rel - position vector in earth-fixed frame (km)
alpha - satellite right ascension (deg)
delta - satellite declination (deg)
n_curves - number of curves comprising the ground track plot
RA - cell array containing the right ascensions for each of
the curves comprising the ground track plot
Dec - cell array containing the declinations for each of
the curves comprising the ground track plot
User M-functions required: sv_from_coe, kepler_E, ra_and_dec_from_r
%}
  
%...Constants
deg = pi/180;
mu = 398600;
hours = 96;
Re = 6378;
J2 = 0.00108263;

%...End data declaration
%...Compute the initial time (since perigee) and
% the rates of node regression and perigee advance

coe = coe_from_sv(r,v,mu);
h = coe(1);
e = coe(2);
RAo = coe(3);
incl = coe(4);
w = coe(5);
TAo = coe(6);
a = coe(7);
T = 2*pi/sqrt(mu)*a^1.5;
n = 2*pi/T
Eo = 2*atan(tan(TAo/2)*sqrt((1-e)/(1+e)));
Mo = Eo - e*sin(Eo);
to = Mo*(T/2/pi);
tf = to + 3600*hours;
n_periods = fix(tf/T);
t_peri = T*(tf/T - fix(tf/T));
M_peri = n*t_peri;
E_peri = kepler_E(e, M_peri);
TA_peri = 2*atan(tan(E_peri/2)*sqrt((1+e)/(1-e)));
rp = (h^2/mu) * (1/(1 + e*cos(TA_peri))) * (cos(TA_peri)*[1;0;0] + sin(TA_peri)*[0;1;0]);
vp = (mu/h) * (-sin(TA_peri)*[1;0;0] + (e + cos(TA_peri))*[0;1;0]);
fac = -3/2*sqrt(mu)*J2*Re^2/(1-e^2)^2/a^(7/2);
Wdot = fac*cos(incl);
wpdot = fac*(5/2*sin(incl)^2 - 2);
RA_peri = RAo + Wdot*3600*hours;
w_peri = w + wpdot*3600*hours;
coe_fin = [h e RA_peri incl w_peri TA_peri a];
[r_fin, v_fin] = sv_from_coe(coe_fin,mu);

display(r_fin);


