function Oblateness_Drag_Code_test
% 
%
% This function solves Example 10.2 by using Encke's method together
% with MATLAB's ode45 to integrate Equation 10.2 for a J2 gravitational
% perturbation given by Equation 10.30.
%
% User M-functions required: sv_from_coe, coe_from_sv, rv_from_r0v0
% User subfunction required: rates
% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
%...Preliminaries:
clc, close all, clear all
%...Conversion factors:
hours = 3600; %Hours to seconds
days = 24*hours; %Days to seconds
deg = pi/180; %Degrees to radians

%...Constants:
global mu
mu = 398600; %Gravitational parameter (km^3/s^2)
RE = 6378; %Earth's radius (km)
J2 = 1082.63e-6;

wE = [ 0 0 7.2921159e-5]'; %Earth's angular velocity (rad/s)
%...Satellite data:




t0 = 0; tf = 200*60; %Initial and final time (s)

%...Obtain the initial state vector from Algorithm 4.5 (sv_from_coe):
R0 = [-4487.58 5428.67 12.252];
V0 = [-2.50989 -2.0271 6.78003];
a0=7.0168e+03;
T0 = 2*pi/sqrt(mu)*a0^1.5;
%R0 is the initial position vector

%R0 is the initial position vector
r0 = norm(R0); v0 = norm(V0); %Magnitudes of T0 and V0
del_t = T0/5000; %Time step for Encke procedure
options = odeset('maxstep', del_t);
%...Begin the Encke integration;
t = t0; %Initialize the time scalar
tsave = t0; %Initialize the vector of solution times
y = [R0 V0]; %Initialize the state vector
del_y0 = zeros(6,1); %Initialize the state vector perturbation
t = t + del_t; %First time step
% Loop over the time interval [t0, tf] with equal increments del_t:
while t <= tf + del_t/2
    
% Integrate Equation 12.7 over the time increment del_t:
[dum,z] = ode45(@rates, [t0 t], del_y0, options);
% Compute the osculating state vector at time t:
[Rosc,Vosc] = rv_from_r0v0(R0, V0, t-t0);
% Rectify:
R0 = Rosc + z(end,1:3);
V0 = Vosc + z(end,4:6);
t0 = t;
% Prepare for next time step:
tsave = [tsave;t];
y = [y; [R0 V0]];
t = t + del_t;
del_y0 = zeros(6,1);
end
% End the loop
t = tsave; %t is the vector of equispaced solution times
%...End the Encke integration;
%...At each solution time extract the orbital elements from the state
% vector using Algorithm 4.2:
n_times = length(t); %n_times is the number of solution times
for j = 1:n_times
R = [y(j,1:3)]
V = [y(j,4:6)]
r(j) = norm(R);
v(j) = norm(V);
coe = coe_from_sv(R,V, mu);
h(j) = coe(1);
e(j) = coe(2);
RA(j) = coe(3);
i(j) = coe(4);
w(j) = coe(5);
TA(j) = coe(6);
end

function dfdt = rates(t,f)
% 
%
% This function calculates the time rates of Encke's deviation in position

% del_r and velocity del_v.
% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
del_r = f(1:3)'; %Position deviation
del_v = f(4:6)'; %Velocity deviation
%...Compute the state vector on the osculating orbit at time t
% (Equation 12.5) using Algorithm 3.4:
[Rosc,Vosc] = rv_from_r0v0(R0, V0, t-t0);
%...Calculate the components of the state vector on the perturbed orbit
% and their magnitudes:
Rpp = Rosc + del_r;
Vpp = Vosc + del_v;
rosc = norm(Rosc);
rpp = norm(Rpp);
B=26.748;
alt = rpp - RE; %Altitude (km)
rho = atmosphere(alt); %Air density from US Standard Model (kf/m^3)
Vrel = Vpp - cross(wE,Rpp); %Velocity relative to the atmosphere (km/s)

vrel = norm(Vrel); %Speed relative to the atmosphere (km/s)
uv = Vrel/vrel; %Relative velocity unit vector
ap = -rho/(0.157)*B*(1000*vrel)^2*uv;  %Acceleration due to drag (m/s^2)
%(converting units of vrel from km/s to m/s)


%...Compute the J2 perturbing acceleration from Equation 12.30:
xx = Rpp(1); yy = Rpp(2); zz = Rpp(3);
fac = 3/2*J2*(mu/rpp^2)*(RE/rpp)^2;
ao = -fac*[(1 - 5*(zz/rpp)^2)*(xx/rpp) ...
(1 - 5*(zz/rpp)^2)*(yy/rpp) ...
(3 - 5*(zz/rpp)^2)*(zz/rpp)];
%...Compute the total perturbing ecceleration from Equation 12.7:
F = 1 - (rosc/rpp)^3;
del_a = -mu/rosc^3*(del_r - F*Rpp) + ao + ap/1000;
dfdt = [del_v(1) del_v(2) del_v(3) del_a(1) del_a(2) del_a(3)]';
dfdt = [del_v del_a]'; %Return the deviative velocity and acceleration
%to ode45.
end %rates
end